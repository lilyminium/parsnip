from collections import defaultdict
import pickle
import itertools

import numpy as np
from rdkit import Chem
import MDAnalysis as mda
import networkx as nx

from .molecule import Molecule, RDMolecule



class Monomer(Molecule):

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
            self.rdmol.AddBond(begin, end, order=bond.GetBondType())
            self.rdmol.RemoveBond(old_begin, old_end)
    

class MonomerUnit(Monomer):

    def __init__(self, monomer):
        self.name = monomer.name  # 
        self.set_rdmol(monomer.rdmol)

    def copy_tags_from(self, other):
        self.tags.extend(other.tags)

    def align_to(self, other, self_atom_indices=[],
                 other_atom_indices=[]):
        atom_map = list(zip(self_atom_indices, other_atom_indices))
        Chem.rdMolAlign.AlignMol(self.rdmol, other.rdmol, atomMap=atom_map)




# class Monomer:

#     @classmethod
#     def from_file(cls, filename, name="UNK"):
#         u = mda.Universe(filename)
#         mol = u.atoms.convert_to("RDKIT")
#         return cls(mol, name)
    
#     def __init__(self, rdmol, name="UNK"):
#         rdmol.UpdatePropertyCache()
#         Chem.GetSymmSSSR(rdmol)
#         self._rdmol = rdmol
#         self._universe = mda.Universe(rdmol, format="RDKIT")
#         pickled = rdmol.ToBinary()
#         self._cymol = PyMonomer(pickled, name)
#         self._atom_names = list(self._universe.atoms.names)

    
        
#     @property
#     def n_atoms(self):
#         return self._cymol.n_atoms

#     @property
#     def n_tags(self):
#         return self._cymol.n_tags

#     @property
#     def name(self):
#         return self._cymol.name

#     @property
#     def n_bonds(self):
#         return self._cymol.n_bonds

#     def view(self):
#         self.minimize()
#         view = nv.show_mdanalysis(self._universe.atoms)
#         view.clear_representations()
#         view.add_representation("ball+stick", selection="all")
#         return view

#     def view_2d(self):
#         return Chem.MolFromSmiles(Chem.MolToSmiles(self._rdmol))

    # def show_tag(self, name, index=-1):
    #     indices = self.get_tag_indices_by_name(name, index)
    #     view = self.view()
    #     view.add_representation("ball+stick", indices,
    #                             opacity=0.3, color="orange",
    #                             radius=0.3)
    #     return view

#     def add_tag(self, tag_name="tag", indices=[]):
#         self._cymol.add_tag(tag_name, indices)

#     def get_tag_indices_by_name(self, name, index=-1):
#         return self._cymol.get_tag_indices_by_name(name, index)


#     def del_atom_by_index(self, index):
#         self._cymol.del_atom_by_index(index)
#         self._atom_names = list(self._universe.atoms.names)

#     def get_atom_indices_by_name(self, names=[]):
#         try:
#             indices = [self._atom_names.index(i) for i in names]
#         except IndexError:
#             raise ValueError("At least one atom name not "
#                              "found. Available names: " + 
#                              ', '.join(self._atom_names))
#         return indices

#     def del_atoms_by_name(self, names=[]):
#         indices = self.get_atom_indices_by_name(names)
#         self.del_atoms_by_indices(indices)
    
#     def del_atoms_by_indices(self, indices=[]):
#         ix = sorted(indices)[::-1]
#         for i in ix:
#             self._cymol.remove_atom_by_index(i, False)
#         self._cymol.reindex_atoms()
        
    
#     def add_tag_by_atom_name(self, tag_name="tag", names=[]):
#         indices = self.get_atom_indices_by_name(names)
#         self.add_tag(tag_name, indices)

#     def write(self, filename):
#         self._universe.atoms.write(filename)

#     def minimize(self, max_iter: int=1000, vdw_threshold: float=10):
#         self._cymol.minimize_rdmol(max_iter, vdw_threshold)