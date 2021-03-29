
import pickle
import copy
import itertools

from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import MMFFSanitizeMolecule, MMFFOptimizeMolecule
import openff
import MDAnalysis as mda

from .monomer import MonomerUnit
# from .lib import PyPolymer
from . import utils, rdutils
from .base import ArrayList
from .alignment import Alignment


class Polymer(MonomerUnit):
    
    def __init__(self, name: str="UNK"):
        super().__init__(Chem.rdchem.RWMol(), name)
        self.units = ArrayList()

    @property
    def n_units(self):
        return len(self.units)

    def add_unit_only(self, unit):
        self._add_unit_atoms(unit)
        self._add_unit_params(unit)
        self._add_unit_tags(unit)

    def _add_unit_atoms(self, unit):
        # https://github.com/rdkit/rdkit/issues/3386
        # self.rdmol.InsertMol(unit.rdmol)
        inc = self.n_atoms

        for atom in unit.rdmol.GetAtoms():
            i = self.rdmol.AddAtom(atom)
        for bond in unit.rdmol.GetBonds():
            i = bond.GetBeginAtomIdx() + inc
            j = bond.GetEndAtomIdx() + inc
            self.rdmol.AddBond(i, j, bond.GetBondType())
        self.atoms.extend(unit.atoms)

    def _add_unit_params(self, unit):
        for attr in ("bonds", "angles", "dihedrals",
                     "impropers", "pairs", "exclusions"):
            pattr = getattr(self, attr)
            uattr = getattr(unit, attr)
            pattr.extend(uattr)

    def _add_unit_tags(self, unit):
        self.tags.extend(unit.tags)

    def add_unit(self, unit: MonomerUnit, unit_tags: list=[],
                 polymer_tags: list=[],
                 replace_polymer_atoms: bool=True,
                 minimize: bool=True):
        # TODO: slow. Can I replace this with try/except?
        if not unit_tags or not polymer_tags or not self.n_atoms:
            self.add_unit_only(unit)

        monomer_tags = utils.asiterable(unit_tags)
        polymer_tags = utils.asiterable(polymer_tags)

        n_mtags = len(monomer_tags)
        err = "Must provide same number of tags for monomer and polymer"
        assert len(polymer_tags) == n_mtags, err

        replace_polymer_atoms = utils.asiterable(replace_polymer_atoms)
        if (len(replace_polymer_atoms) != n_mtags
            and len(replace_polymer_atoms) == 1):
            replace_polymer_atoms = replace_polymer_atoms * n_mtags
        
        mp_indices = []
        del_rep_atoms = []
        for mtag, ptag, rep in zip(monomer_tags, polymer_tags,
                                   replace_polymer_atoms):
            mp_indices.extend(list(zip(mtag.indices, ptag.indices)))

            rep = utils.asiterable(rep)
            if len(rep) != len(mtag.atoms) and len(rep) == 1:
                rep = rep * len(mtag.atoms)
            for m, p, r in zip(mtag.atoms, ptag.atoms, rep):
                atoms = [p, m]
                if not r:
                    atoms = atoms[::-1]
                del_rep_atoms.append(atoms)

        Chem.rdMolAlign.AlignMol(unit.rdmol, self.rdmol,
                                 atomMap=mp_indices)
        self.add_unit_only(unit)

        del_rep_atoms = sorted(del_rep_atoms, key=lambda x: x[0].index)
        delete_atoms = [x[0] for x in del_rep_atoms]

        self.remove_params_within_atoms(delete_atoms)

        for to_delete, to_replace in del_rep_atoms[::-1]:
            self.update_atom(to_delete, to_replace)
        
        self.remove_atoms(delete_atoms[::-1])
        self.clean()

    def add_monomer_with_tag_names(self, monomer: Monomer,
                                   monomer_tag_names: list=[],
                                   polymer_tag_names: list=[],
                                   monomer_tag_indices=-1,
                                   polymer_tag_indices=-1,
                                   replace_polymer_atoms: bool=True,
                                   minimize: bool=True):
        unit = MonomerUnit(monomer)

        mnames = utils.asiterable(monomer_tag_names)
        pnames = utils.asiterable(polymer_tag_names)
        
        n_mtags = len(mnames)

        mix = utils.asiterable(monomer_tag_indices)
        pix = utils.asiterable(polymer_tag_indices)
        if len(pix) != n_mtags and len(pix) == 1:
            pix = pix * n_mtags
        if len(mix) != n_mtags and len(mix) == 1:
            mix = mix * n_mtags
        mtags = [unit.get_tag(x, i)
                 for x, i in zip(mnames, mix)]
        ptags = [self.get_tag(x, i)
                 for x, i in zip(pnames, pix)]
        
        return self.add_unit(unit, unit_tags=mtags, polymer_tags=ptags,
                             replace_polymer_atoms=replace_polymer_atoms,
                             minimize=minimize)

    def add_monomer(self, monomer: Monomer,
                    monomer_tags: list=[],
                    polymer_tags: list=[],
                    **kwargs):
        """Bit of an annoying workaround"""

        mtags = utils.asiterable(monomer_tags)
        ptags = utils.asiterable(polymer_tags)

        mnames = [x.name for x in mtags]
        pnames = [x.name for x in ptags]

        mindices = [monomer.tags[x.name].index(x) for x in mtags]
        pindices = [self.tags[x.name].index(x) for x in ptags]

        return self.add_monomer_with_tag_names(monomer,
                                               monomer_tag_names=mnames,
                                               polymer_tag_indices=pnames,
                                               monomer_tag_indices=mindices,
                                               polymer_tag_indices=pindices,
                                               **kwargs)


    # def add_monomers(self, monomers: list=[], ratio: list=[],
    #                  monomer_tag_names: list=[],
    #                  polymer_tag_names: list=[],
    #                  replace_polymer_atoms: list=[], total: int=1,
    #                  random: bool=False, minimize: bool=True):
    #     ...

    def get_capped_rdunits(self, n_neighbors: int=3):
        self.clean_unit_atoms()
        self.reindex_atoms()

        capped_units = {}
        smiles = []
        for unit in self.units:
            capped = self.cap_unit(unit, n_neighbors)
            smi = Chem.MolToSmiles(capped)
            capped_units[smi] = capped
            smiles.append(smi)
        
        keys, values = zip(capped_units.items())
        keys = list(keys)
        indices = [keys.index(x) for x in smiles]
        return values, indices
        
        


    def cap_unit(self, unit, n_neighbors=3):
        indices = sorted(unit.atoms.indices)
        capped = Chem.rdchem.RWMol()
        n_indices = len(indices)
        if not n_indices:
            return capped
        neighbor_ints = rdutils.get_neighbor_ints(self.rdmol, indices,
                                                  n_neighbors)
        new_conf = Chem.rdchem.Conformer()
        conf = self.rdmol.GetConformer(0)
        capped.AddConformer(new_conf)
        only_on_atoms = []

        # add unit atoms
        for i, ix in enumerate(indices):
            at = self.rdmol.GetAtomWithIdx(ix)
            capped.AddAtom(at)
            new_conf.SetAtomPosition(i, conf.GetAtomPosition(ix))
        
        # add neighbor atoms with AltLoc "+"
        for i, ix in enumerate(neighbor_ints, n_indices):
            at = self.rdmol.GetAtomWithIdx(ix)
            old_info = at.GetPDBResidueInfo()
            try:
                info = Chem.AtomPDBResidueInfo(old_info)
            except TypeError:  # I think
                info = Chem.AtomPDBResidueInfo()
            info.SetAltLoc("+")
            at.SetPDBResidueInfo(info)
            at.SetNoImplicit(False)

            capped.AddAtom(at)
            new_conf.SetAtomPosition(i, conf.GetAtomPosition(ix))
            only_on_atoms.append(i)
        
        # add bonds
        all_ints = indices + neighbor_ints
        capped_ints = np.arange(len(all_ints))
        for i, j in itertools.combinations(capped_ints, 2):
            bond = self.rdmol.GetBondBetweenAtoms(all_ints[i], all_ints[j])
            if bond is not None:
                capped.AddBond(i, j, bond.GetBondType())

        # add Hs
        hcapped = Chem.AddHs(capped, explicitOnly=False, addCoords=True,
                             onlyOnAtoms=only_on_atoms)
        info_template = hcapped.GetAtomWithIdx(0).GetPDBResidueInfo()
        for i in range(capped_ints[-1], hcapped.GetNumAtoms()):
            at = hcapped.GetAtomWithIdx(i)
            new_info = Chem.AtomPDBResidueInfo(info_template)
            new_info.SetName("H")
            new_info.SetAltLoc("-")
            at.SetPDBResidueInfo(new_info)
        
        Chem.SanitizeMol(hcapped)
        return hcapped




# class Polymer(Monomer):
    
#     def __init__(self, name: str="UNK"):
#         self._cymol = PyPolymer(name)
#         self._n_neighbors = 3

#     @property
#     def n_units(self):
#         return self._cymol.n_units

#     @property
#     def _universe(self):
#         if not self.n_atoms:
#             return mda.Universe.empty(0)
#         u = mda.Universe(self._rdmol, format="RDKIT")
#         return u

#     @property
#     def _rdmol(self):
#         molbin = self._cymol.get_rdmol_binary()
#         return Chem.Mol(molbin)

#     def add_monomer(self, monomer: Monomer, monomer_tag: list=[],
#                     polymer_tag: list=[],
#                     replace_polymer_atoms: bool=True,
#                     minimize: bool=True):
#         if not monomer_tag or not polymer_tag or not self.n_atoms:
#             self._cymol.add_monomer_only(monomer._cymol)
#             return

#         monomer_tag = utils.asiterable(monomer_tag)
#         polymer_tag = utils.asiterable(polymer_tag)
#         if not utils.isiterable(replace_polymer_atoms):
#             replace_polymer_atoms = [replace_polymer_atoms] * len(monomer_tag)
        
#         if len(replace_polymer_atoms) != len(monomer_tag):
#             if len(replace_polymer_atoms) == 1:
#                 replace_polymer_atoms = replace_polymer_atoms * len(monomer_tag)

#         # create alignment structures
#         alignments = []
#         for mon, pol, rep in zip(monomer_tag, polymer_tag,
#                                  replace_polymer_atoms):
#             n_monomer_atoms = len(monomer.get_tag_indices_by_name(mon))
#             if not utils.isiterable(rep):
#                 rep = [rep] * n_monomer_atoms
#             alignments.append([[mon, pol], rep])

#         self._cymol.add_monomer_with_tags(monomer._cymol, alignments)
#         if minimize:
#             self.minimize(50)

#     def add_monomers(self, monomers: list=[], ratio: list=[],
#                      monomer_tags: list=[], polymer_tags: list=[],
#                      replace_polymer_atoms: list=[], total: int=1,
#                      random: bool=False, minimize: bool=True):
#         """
#         Add monomers.

#         Example
#         -------

#             pol.add_monomers([mon1, mon2], monomer_tags=[["left"], ["left"]],
#                              polymer_tags=[["right"], ["right"]],
#                              replace_polymer_atoms=[[True], [True]],
#                              ratio=[0.3, 0.7], total=100)

#             pol.add_monomers([mon1, mon2], monomer_tags="left",
#                              polymer_tags="right", replace_polymer_atoms=True,
#                              ratio=[0.3, 0.7])
#         """

#         # my god there's a lot of checks
#         # can I pass this off to C++...?
#         monomer_tags = utils.asiterable(monomer_tags)
#         polymer_tags = utils.asiterable(polymer_tags)
#         replace_polymer_atoms = utils.asiterable(replace_polymer_atoms)

#         n_monomers = len(monomers)
#         if len(monomer_tags) == 1:
#             monomer_tags = monomer_tags * n_monomers
#         if len(polymer_tags) == 1:
#             polymer_tags = polymer_tags * n_monomers
#         if len(replace_polymer_atoms) == 1:
#             replace_polymer_atoms = replace_polymer_atoms * n_monomers
        
#         ferr = "Must specify ratio for each monomer alignment"
#         assert len(ratio) == n_monomers, ferr
#         merr = "Must specify monomer_tags for each monomer alignment"
#         assert len(monomer_tags) == n_monomers, merr
#         perr = "Must specify polymer_tags for each monomer alignment"
#         assert len(polymer_tags) == n_monomers, merr
#         rerr = "Must specify replace_polymer_atoms for each monomer alignment"
#         assert len(replace_polymer_atoms) == n_monomers, merr

#         monomer_tags = [utils.asiterable(x) for x in monomer_tags]
#         polymer_tags = [utils.asiterable(x) for x in polymer_tags]
#         replace_polymer_atoms = [utils.asiterable(x) for x 
#                                  in replace_polymer_atoms]

#         # create alignment structures
#         # AtomwiseTagNameAlignment:  pair[  pair[string, string], vector[bool]  ]
#         # vector[  pair[  Monomer*, vector[AtomwiseTagNameAlignment]  ]  ]
#         alignments = []
#         for mon, mtags, ptags, repl in zip(monomers, monomer_tags,
#                                            polymer_tags,
#                                            replace_polymer_atoms):
#             # make sure tags are the same length
#             n_mtags = len(mtags)
#             n_ptags = len(ptags)
#             n_tags = max([n_mtags, n_ptags])
#             if n_tags > 1:
#                 if n_mtags == 1:
#                     mtags = mtags * n_tags
#                 if n_ptags == 1:
#                     ptags = ptags * n_tags
#                 if len(repl) == 1:
#                     repl = repl * n_tags
            
#             atomwise_vec = []
#             for mtag, ptag, rep in zip(mtags, ptags, repl):
#                 n_monomer_atoms = len(mon.get_tag_indices_by_name(mtag))
#                 if not utils.isiterable(rep):
#                     rep = [rep] * n_monomer_atoms

#                 atomwise = [[mtag, ptag], rep]
#                 atomwise_vec.append(atomwise)
            
#             alignments.append([mon._cymol, atomwise_vec])
        
#         self._cymol.add_monomers(alignments, ratio, total, random)
#         # if minimize:
#         #     self.minimize(50)
        
#     def get_unit_indices(self):
#         # TODO: err, come up with better fix than this
#         return self._cymol.get_unit_indices(self._n_neighbors)

#     def get_capped_rdunits(self, n_neighbors: int=3):
#         self._n_neighbors = n_neighbors
#         bins = self._cymol.get_capped_rdunits_binary(n_neighbors)
#         mols = [Chem.Mol(x) for x in bins]
#         hmols = []

#         # idk why I can't do this in C++ :(
#         for i, m in enumerate(mols):
#             Chem.Cleanup(m)
#             Chem.SanitizeMol(m)
#             u = mda.Universe(m, format="RDKIT")
#             ix = [int(x) for x in u.select_atoms("altLoc +").indices]
            
#             newmol = Chem.AddHs(m, explicitOnly=False, addCoords=True)#, onlyOnAtoms=ix)

#             info = newmol.GetAtomWithIdx(0).GetMonomerInfo()

#             for i in range(len(u.atoms), newmol.GetNumAtoms()):
#                 atom = newmol.GetAtomWithIdx(i)
#                 info2 = type(info)("H", i+1, "-", info.GetResidueName(),
#                                    info.GetResidueNumber())
                
#                 atom.SetMonomerInfo(info2)

#             Chem.SanitizeMol(newmol)
#             MMFFSanitizeMolecule(newmol)
#             MMFFOptimizeMolecule(newmol)
#             hmols.append(newmol)
#         return hmols


#     def get_capped_offunits(self, n_neighbors: int=3):
#         mols = self.get_capped_rdunits(n_neighbors)
#         offmols = [openff.toolkit.topology.Molecule.from_rdkit(x) for x in mols]
#         return offmols


#     def parametrize_capped_units(self, n_neighbors: int=3):
#         offmols = self.get_capped_offunits(n_neighbors)
#         top = openff.toolkit.topology.Topology.from_molecules(offmols)
#         return top
 