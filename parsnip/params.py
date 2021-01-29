
import numpy as np

from .atom import Atom

class Param:
    n_atoms = 0
    def __init__(self, atoms: list=[]):
        self._atoms = np.empty(self.n_atoms, dtype=object)
        self.atoms = atoms

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, values):
        err = f"expected {self.n_atoms} atoms but {len(values)} were given"
        assert len(values) == self.n_atoms, err
        if values[0].index > values[-1].index:
            values = values[::-1]
        self._atoms[:] = values

    @property
    def indices(self):
        return np.array([atom.index for atom in self.atoms])
    
    def update_atom(self, old_atom, new_atom):
        self.atoms = np.where(self.atoms == old_atom, new_atom, self.atoms)

    def contains_atom(self, atom):
        return atom in self.atoms
    
    def within_atoms(self, atoms):
        # TODO: can numpy do sets better??
        atoms = np.fromiter(atoms, Atom)
        return np.all(np.isin(self.atoms, atoms))

    def _within_atom_array(self, atoms):
        return np.all(np.isin(self.atoms, atoms))
    

class Bond(Param):
    n_atoms = 2

class Pair(Param):
    n_atoms = 2

class Exclusion(Param):
    n_atoms = 2

class Angle(Param):
    n_atoms = 3

class Dihedral(Param):
    n_atoms = 4

class Improper(Param):
    n_atoms = 4

