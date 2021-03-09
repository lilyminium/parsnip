
import pickle

from rdkit import Chem

import MDAnalysis as mda

from .monomer import Monomer
from .lib import PyPolymer
from . import utils

class Polymer(Monomer):
    
    def __init__(self, name: str="UNK"):
        self._cymol = PyPolymer(name)

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
                    replace_polymer_atoms: bool=True):
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

    

    def get_capped_rdunits(self, n_neighbors: int=3):
        bins = self._cymol.get_capped_rdunits_binary(n_neighbors)
        mols = [Chem.Mol(x) for x in bins]
        hmols = []

        # idk why I can't do this in C++ :(
        for m in mols:
            Chem.SanitizeMol(m)
            u = mda.Universe(m, format="RDKIT")
            ix = [int(x) for x in u.select_atoms("icode +").indices]
            newmol = Chem.AddHs(m, explicitOnly=True, addCoords=True, onlyOnAtoms=ix)
            new_u = mda.Universe(newmol)
            new_u.atoms[len(u.atoms):].altLocs = "-"
            hmols.append(new_u.atoms.convert_to("RDKIT"))
            # hmols.append(newmol)
        return hmols
 