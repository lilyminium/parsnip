
import pickle

import MDAnalysis as mda

from .lib import PyMonomer

class Monomer:
    
    def __init__(self, rdmol, name="UNK"):
        self._universe = mda.Universe(rdmol, format="RDKIT")
        pickled = rdmol.ToBinary()
        self._pymonomer = PyMonomer(pickled, name)
        
    @property
    def n_atoms(self):
        return self._pymonomer.n_atoms

    @property
    def n_tags(self):
        return self._pymonomer.n_tags

    @property
    def name(self):
        return self._pymonomer.name

    @property
    def n_bonds(self):
        return self._pymonomer.n_bonds

    def add_tag(self, tag_name="tag", indices=[]):
        self._pymonomer.addTag(tag_name, indices)
    
    def add_tag_by_atom_name(self, tag_name="tag", names=[]):
        all_names = list(self._universe.atoms.names)
        try:
            indices = [all_names.index(i) for i in names]
        except IndexError:
            raise ValueError("At least one atom name not "
                             "found. Available names: " + 
                             ', '.join(all_names))
        else:
            self.add_tag(tag_name, indices)
