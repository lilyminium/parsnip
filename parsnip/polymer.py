
import pickle
import copy

from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import MMFFSanitizeMolecule, MMFFOptimizeMolecule
import openff
import MDAnalysis as mda

from .monomer import Monomer
from .lib import PyPolymer
from . import utils

class Polymer(Monomer):
    
    def __init__(self, name: str="UNK"):
        self._cymol = PyPolymer(name)
        self._n_neighbors = 3

    @property
    def n_units(self):
        return self._cymol.n_units

    @property
    def _universe(self):
        if not self.n_atoms:
            return mda.Universe.empty(0)
        u = mda.Universe(self._rdmol, format="RDKIT")
        return u

    @property
    def _rdmol(self):
        molbin = self._cymol.get_rdmol_binary()
        return Chem.Mol(molbin)

    def add_monomer(self, monomer: Monomer, monomer_tag: list=[],
                    polymer_tag: list=[],
                    replace_polymer_atoms: bool=True,
                    minimize: bool=True):
        if not monomer_tag or not polymer_tag or not self.n_atoms:
            self._cymol.add_monomer_only(monomer._cymol)
            return

        monomer_tag = utils.asiterable(monomer_tag)
        polymer_tag = utils.asiterable(polymer_tag)
        if not utils.isiterable(replace_polymer_atoms):
            replace_polymer_atoms = [replace_polymer_atoms] * len(monomer_tag)
        
        if len(replace_polymer_atoms) != len(monomer_tag):
            if len(replace_polymer_atoms) == 1:
                replace_polymer_atoms = replace_polymer_atoms * len(monomer_tag)

        # create alignment structures
        alignments = []
        for mon, pol, rep in zip(monomer_tag, polymer_tag,
                                 replace_polymer_atoms):
            n_monomer_atoms = len(monomer.get_tag_indices_by_name(mon))
            if not utils.isiterable(rep):
                rep = [rep] * n_monomer_atoms
            alignments.append([[mon, pol], rep])

        self._cymol.add_monomer_with_tags(monomer._cymol, alignments)
        if minimize:
            self.minimize(50)

    def add_monomers(self, monomers: list=[], ratio: list=[],
                     monomer_tags: list=[], polymer_tags: list=[],
                     replace_polymer_atoms: list=[], total: int=1,
                     random: bool=False, minimize: bool=True):
        """
        Add monomers.

        Example
        -------

            pol.add_monomers([mon1, mon2], monomer_tags=[["left"], ["left"]],
                             polymer_tags=[["right"], ["right"]],
                             replace_polymer_atoms=[[True], [True]],
                             ratio=[0.3, 0.7], total=100)

            pol.add_monomers([mon1, mon2], monomer_tags="left",
                             polymer_tags="right", replace_polymer_atoms=True,
                             ratio=[0.3, 0.7])
        """

        # my god there's a lot of checks
        # can I pass this off to C++...?
        monomer_tags = utils.asiterable(monomer_tags)
        polymer_tags = utils.asiterable(polymer_tags)
        replace_polymer_atoms = utils.asiterable(replace_polymer_atoms)

        n_monomers = len(monomers)
        if len(monomer_tags) == 1:
            monomer_tags = monomer_tags * n_monomers
        if len(polymer_tags) == 1:
            polymer_tags = polymer_tags * n_monomers
        if len(replace_polymer_atoms) == 1:
            replace_polymer_atoms = replace_polymer_atoms * n_monomers
        
        ferr = "Must specify ratio for each monomer alignment"
        assert len(ratio) == n_monomers, ferr
        merr = "Must specify monomer_tags for each monomer alignment"
        assert len(monomer_tags) == n_monomers, merr
        perr = "Must specify polymer_tags for each monomer alignment"
        assert len(polymer_tags) == n_monomers, merr
        rerr = "Must specify replace_polymer_atoms for each monomer alignment"
        assert len(replace_polymer_atoms) == n_monomers, merr

        monomer_tags = [utils.asiterable(x) for x in monomer_tags]
        polymer_tags = [utils.asiterable(x) for x in polymer_tags]
        replace_polymer_atoms = [utils.asiterable(x) for x 
                                 in replace_polymer_atoms]

        # create alignment structures
        # AtomwiseTagNameAlignment:  pair[  pair[string, string], vector[bool]  ]
        # vector[  pair[  Monomer*, vector[AtomwiseTagNameAlignment]  ]  ]
        alignments = []
        for mon, mtags, ptags, repl in zip(monomers, monomer_tags,
                                           polymer_tags,
                                           replace_polymer_atoms):
            # make sure tags are the same length
            n_mtags = len(mtags)
            n_ptags = len(ptags)
            n_tags = max([n_mtags, n_ptags])
            if n_tags > 1:
                if n_mtags == 1:
                    mtags = mtags * n_tags
                if n_ptags == 1:
                    ptags = ptags * n_tags
                if len(repl) == 1:
                    repl = repl * n_tags
            
            atomwise_vec = []
            for mtag, ptag, rep in zip(mtags, ptags, repl):
                n_monomer_atoms = len(mon.get_tag_indices_by_name(mtag))
                if not utils.isiterable(rep):
                    rep = [rep] * n_monomer_atoms

                atomwise = [[mtag, ptag], rep]
                atomwise_vec.append(atomwise)
            
            alignments.append([mon._cymol, atomwise_vec])
        
        self._cymol.add_monomers(alignments, ratio, total, random)
        # if minimize:
        #     self.minimize(50)
        
    def get_unit_indices(self):
        # TODO: err, come up with better fix than this
        return self._cymol.get_unit_indices(self._n_neighbors)

    def get_capped_rdunits(self, n_neighbors: int=3):
        self._n_neighbors = n_neighbors
        bins = self._cymol.get_capped_rdunits_binary(n_neighbors)
        mols = [Chem.Mol(x) for x in bins]
        hmols = []

        # idk why I can't do this in C++ :(
        for i, m in enumerate(mols):
            Chem.Cleanup(m)
            Chem.SanitizeMol(m)
            u = mda.Universe(m, format="RDKIT")
            ix = [int(x) for x in u.select_atoms("altLoc +").indices]
            
            newmol = Chem.AddHs(m, explicitOnly=False, addCoords=True)#, onlyOnAtoms=ix)

            info = newmol.GetAtomWithIdx(0).GetMonomerInfo()

            for i in range(len(u.atoms), newmol.GetNumAtoms()):
                atom = newmol.GetAtomWithIdx(i)
                info2 = type(info)("H", i+1, "-", info.GetResidueName(),
                                   info.GetResidueNumber())
                
                atom.SetMonomerInfo(info2)

            Chem.SanitizeMol(newmol)
            MMFFSanitizeMolecule(newmol)
            MMFFOptimizeMolecule(newmol)
            hmols.append(newmol)
        return hmols


    def get_capped_offunits(self, n_neighbors: int=3):
        mols = self.get_capped_rdunits(n_neighbors)
        offmols = [openff.toolkit.topology.Molecule.from_rdkit(x) for x in mols]
        return offmols


    def parametrize_capped_units(self, n_neighbors: int=3):
        offmols = self.get_capped_offunits(n_neighbors)
        top = openff.toolkit.topology.Topology.from_molecules(offmols)
        return top
 