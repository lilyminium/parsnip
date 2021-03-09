
import pickle

import MDAnalysis as mda
import nglview as nv
from rdkit import Chem

from .lib import PyMonomer

class Monomer:

    @classmethod
    def from_file(cls, filename, name="UNK"):
        u = mda.Universe(filename)
        mol = u.atoms.convert_to("RDKIT")
        return cls(mol, name)
    
    def __init__(self, rdmol, name="UNK"):
        rdmol.UpdatePropertyCache()
        Chem.GetSymmSSSR(rdmol)
        self._universe = mda.Universe(rdmol, format="RDKIT")
        pickled = rdmol.ToBinary()
        self._cymol = PyMonomer(pickled, name)
        self._atom_names = list(self._universe.atoms.names)
        
    @property
    def n_atoms(self):
        return self._cymol.n_atoms

    @property
    def n_tags(self):
        return self._cymol.n_tags

    @property
    def name(self):
        return self._cymol.name

    @property
    def n_bonds(self):
        return self._cymol.n_bonds

    def view(self):
        self.minimize()
        view = nv.show_mdanalysis(self._universe.atoms)
        view.clear_representations()
        view.add_representation("ball+stick", selection="all")
        return view

    def show_tag(self, name, index=-1):
        indices = self.get_tag_indices_by_name(name, index)
        view = self.view()
        view.add_representation("ball+stick", indices,
                                opacity=0.3, color="orange",
                                radius=0.3)
        return view

    def add_tag(self, tag_name="tag", indices=[]):
        self._cymol.add_tag(tag_name, indices)

    def get_tag_indices_by_name(self, name, index=-1):
        return self._cymol.get_tag_indices_by_name(name, index)


    def del_atom_by_index(self, index):
        self._cymol.del_atom_by_index(index)
        self._atom_names = list(self._universe.atoms.names)

    def get_atom_indices_by_name(self, names=[]):
        try:
            indices = [self._atom_names.index(i) for i in names]
        except IndexError:
            raise ValueError("At least one atom name not "
                             "found. Available names: " + 
                             ', '.join(self._atom_names))
        return indices

    def del_atoms_by_name(self, names=[]):
        indices = self.get_atom_indices_by_name(names)
        self.del_atoms_by_indices(indices)
    
    def del_atoms_by_indices(self, indices=[]):
        ix = sorted(indices)[::-1]
        for i in ix:
            self._cymol.remove_atom_by_index(i, False)
        self._cymol.reindex_atoms()
        
    
    def add_tag_by_atom_name(self, tag_name="tag", names=[]):
        indices = self.get_atom_indices_by_name(names)
        self.add_tag(tag_name, indices)

    def write(self, filename):
        self._universe.atoms.write(filename)

    def minimize(self, max_iter: int=1000, vdw_threshold: float=10):
        self._cymol.minimize_rdmol(max_iter, vdw_threshold)