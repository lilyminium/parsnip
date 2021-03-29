
import pickle
import copy
import itertools

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdForceFieldHelpers import MMFFSanitizeMolecule, MMFFOptimizeMolecule
# import openff
import numpy as np
import MDAnalysis as mda

from .monomer import Monomer, MonomerUnit
from . import utils, rdutils
from .base import ArrayList


class Polymer(Monomer):
    
    def __init__(self, name: str="UNK"):
        mol = Chem.rdchem.RWMol()
        super().__init__(mol, name=name)
        self.rdconf = Chem.rdchem.Conformer()
        self.rdmol.AddConformer(self.rdconf)
        self.units = ArrayList()

    @property
    def n_units(self):
        return len(self.units)

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, values):
        self._units = values

    def add_unit_only(self, unit):
        self._add_unit_atoms(unit)
        self._add_unit_params(unit)
        self._add_unit_tags(unit)
        self.units.append(unit)
        unit.resid = len(self.units)

    def _add_unit_atoms(self, unit):
        # https://github.com/rdkit/rdkit/issues/3386
        # self.rdmol.InsertMol(unit.rdmol)
        inc = self.n_atoms
        uconf = unit.rdmol.GetConformer()
        pconf = self.rdconf
        
        self.rdconf = Chem.rdchem.Conformer(inc + unit.rdmol.GetNumAtoms())
        copied = Chem.rdchem.RWMol(unit.rdmol)
        for i, atom in enumerate(copied.GetAtoms()):
            j = self.rdmol.AddAtom(atom)            
            self.rdconf.SetAtomPosition(j, list(uconf.GetAtomPosition(i)))
        for bond in unit.rdmol.GetBonds():
            i = bond.GetBeginAtomIdx() + inc
            j = bond.GetEndAtomIdx() + inc
            self.rdmol.AddBond(i, j, bond.GetBondType())

        for i in range(inc):
            self.rdconf.SetAtomPosition(i, list(pconf.GetAtomPosition(i)))
        self.rdmol.RemoveConformer(0)
        self.rdmol.AddConformer(self.rdconf)
        self.atoms.extend(unit.atoms)

        for rdatom, pyatom in zip(self.rdmol.GetAtoms(), self.atoms):
            pyatom.rdatom = rdatom

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
            return

        monomer_tags = utils.asiterable(unit_tags)
        polymer_tags = utils.asiterable(polymer_tags)

        n_mtags = len(monomer_tags)
        err = "Must provide same number of tags for monomer and polymer"
        assert len(polymer_tags) == n_mtags, err

        replace_polymer_atoms = utils.asiterable(replace_polymer_atoms)
        if (len(replace_polymer_atoms) != n_mtags
            and len(replace_polymer_atoms) == 1):
            replace_polymer_atoms = replace_polymer_atoms * n_mtags
        
        m_indices = []
        p_indices = []
        del_rep_atoms = []
        for mtag, ptag, rep in zip(monomer_tags, polymer_tags,
                                   replace_polymer_atoms):
            m_indices.extend(mtag.indices)
            p_indices.extend(ptag.indices)

            rep = utils.asiterable(rep)
            if len(rep) != len(mtag.atoms) and len(rep) == 1:
                rep = rep * len(mtag.atoms)
            for m, p, r in zip(mtag.atoms, ptag.atoms, rep):
                atoms = [p, m]
                if not r:
                    atoms = atoms[::-1]
                del_rep_atoms.append(atoms)

        # great
        atom_map = tuple(zip(map(int, m_indices), map(int, p_indices)))

        AllChem.AlignMol(unit.rdmol, self.rdmol, 0, 0, atomMap=atom_map)
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
                                               polymer_tag_names=pnames,
                                               monomer_tag_indices=mindices,
                                               polymer_tag_indices=pindices,
                                               **kwargs)

    def clean_unit_atoms(self):
        for unit in self.units:
            mask = np.isin(unit.atoms, self.atoms)
            unit.atoms = unit.atoms[mask]


    # def add_monomers(self, monomers: list=[], ratio: list=[],
    #                  monomer_tag_names: list=[],
    #                  polymer_tag_names: list=[],
    #                  replace_polymer_atoms: list=[], total: int=1,
    #                  random: bool=False, minimize: bool=True):
    #     ...

    def get_capped_rdunits(self, n_neighbors: int=3, minimize: bool=True):
        self.clean_unit_atoms()
        self.reindex_atoms()
        if minimize:
            self.minimize()

        capped_units = {}
        smiles = []
        for unit in self.units:
            capped = self.cap_unit(unit, n_neighbors)
            smi = Chem.MolToSmiles(capped)
            capped_units[smi] = capped
            smiles.append(smi)
        
        keys, values = zip(*capped_units.items())
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
        new_conf = Chem.rdchem.Conformer(n_indices + len(neighbor_ints))
        conf = self.rdmol.GetConformer(0)
        only_on_atoms = []

        # add unit atoms
        indices = list(map(int, indices))
        for ix in indices:
            at = self.rdmol.GetAtomWithIdx(ix)
            new = capped.AddAtom(at)
            new_conf.SetAtomPosition(new, list(conf.GetAtomPosition(ix)))
        
        # add neighbor atoms with AltLoc "+"

        for ix in neighbor_ints:
            at = self.rdmol.GetAtomWithIdx(ix)

            new = capped.AddAtom(at)
            at = capped.GetAtomWithIdx(new)
            old_info = at.GetPDBResidueInfo()
            info = rdutils.copy_pdbinfo(old_info)
            info.SetAltLoc("+")
            at.SetMonomerInfo(info)
            at.SetNoImplicit(False)
            new_conf.SetAtomPosition(new, list(conf.GetAtomPosition(ix)))
            only_on_atoms.append(new)
        
        # add bonds
        all_ints = indices + neighbor_ints
        # things RDKit complains about: np.uint,  np.uint32, np.uintp, np.uintc
        # capped_ints = np.arange(len(all_ints), dtype=np.uintc)
        capped_ints = list(map(int, range(len(all_ints))))  # omg
        for i, j in itertools.combinations(capped_ints, 2):
            bond = self.rdmol.GetBondBetweenAtoms(all_ints[i], all_ints[j])
            if bond is not None:
                capped.AddBond(i, j, bond.GetBondType())

        capped.AddConformer(new_conf)
        # add Hs
        Chem.SanitizeMol(capped)
        hcapped = Chem.AddHs(capped, explicitOnly=False, addCoords=True,
                             onlyOnAtoms=only_on_atoms)
        info_template = hcapped.GetAtomWithIdx(0).GetPDBResidueInfo()
        for i in range(capped_ints[-1], hcapped.GetNumAtoms()):
            at = hcapped.GetAtomWithIdx(i)
            new_info = rdutils.copy_pdbinfo(info_template)
            new_info.SetName("H")
            new_info.SetAltLoc("-")
            at.SetMonomerInfo(new_info)
        
        Chem.SanitizeMol(hcapped)
        return hcapped


