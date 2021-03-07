
import pickle

from rdkit import Chem

import MDAnalysis as mda

from .monomer import Monomer
from .lib import PyPolymer

class Polymer(Monomer):
    
    def __init__(self, name="UNK"):
        self._cymol = PyPolymer(name)

    @property
    def _universe(self):
        if not self.n_atoms:
            return mda.Universe.empty(0)
        molbin = self._cymol.get_rdmol_binary()
        rdmol = Chem.Mol(molbin)
        u = mda.Universe(rdmol, format="RDKIT")
        return u

    def add_monomer(self, monomer, monomer_tag=None,
                    polymer_tag=None, replace_polymer_atoms=True):
        if not monomer_tag or not polymer_tag or not self.n_atoms:
            self._cymol.add_monomer_only(monomer._cymol)
        else:
            self._cymol.add_monomer(monomer._cymol, monomer_tag, polymer_tag,
                                    replace_polymer_atoms)

    def get_capped_rdunits(self, n_neighbors=3):
        bins = self._cymol.get_capped_rdunits_binary(n_neighbors)
        mols = [Chem.Mol(x) for x in bins]
        return mols
 