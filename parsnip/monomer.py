from collections import defaultdict
import pickle
import itertools

import numpy as np
from rdkit import Chem
import MDAnalysis as mda
import networkx as nx

from .molecule import Molecule, RDMolecule



class Monomer(RDMolecule):

    @classmethod
    def from_file(cls, filename, name="UNK"):
        u = mda.Universe(filename)
        rdmol = u.atoms.convert_to("RDKIT")
        return cls(rdmol, name=name)

    def __init__(self, rdmol, name: str="UNK"):
        super().__init__(rdmol)
        self.name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name
        for atom in self.atoms:
            atom.resname = name

    def update_atom(self, old_atom, new_atom):
        for param_list in self.params:
            # TODO: more efficient way to do this
            for param in param_list:
                param.update_atom(old_atom, new_atom)
        
        old_ix = old_atom.rdatom.GetIdx()
        new_ix = new_atom.rdatom.GetIdx()
        for bond in old_atom.rdatom.GetBonds():
            old_begin = bond.GetBeginAtomIdx()
            old_end = bond.GetEndAtomIdx()
            if old_begin == old_ix:
                begin = new_ix
                end = old_end
            else:
                begin = old_begin
                end = new_ix
            existing = self.rdmol.GetBondBetweenAtoms(begin, end)
            if existing is None:
                self.rdmol.AddBond(begin, end, order=bond.GetBondType())
            self.rdmol.RemoveBond(old_begin, old_end)

    def align_to(self, other, self_atom_indices=[],
                 other_atom_indices=[]):
        atom_map = list(zip(self_atom_indices, other_atom_indices))
        Chem.rdMolAlign.AlignMol(self.rdmol, other.rdmol, atomMap=atom_map)

    

class MonomerUnit(Monomer):

    def __init__(self, monomer):
        super().__init__(Chem.rdchem.RWMol(monomer.rdmol), monomer.name)
        self.copy_tags_from(monomer)
        for p in self.params:
            p.clear()
        self.copy_params_from(monomer)
        self.resid = 1

    @property
    def resid(self):
        return self._resid

    @resid.setter
    def resid(self, value):
        self._resid = int(value)
        for at in self.atoms:
            at.resid = self._resid

    
