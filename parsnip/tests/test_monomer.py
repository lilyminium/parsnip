
import pytest
import rdkit

from rdkit.Chem.rdmolfiles import MolFromPDBFile
from parsnip.monomer import Monomer

import MDAnalysis as mda

from parsnip.tests.datafiles import m1c1

def pdb2rdmol(pdb):
    return MolFromPDBFile(pdb, False, False)

class TestMonomer:

    @pytest.mark.parametrize("pdb, n_atoms, n_bonds", [
        (m1c1, 54, 53),
    ])
    def test_create_monomer(self, pdb, n_atoms, n_bonds):
        mol = pdb2rdmol(pdb)
        assert mol.GetNumAtoms() == n_atoms
        monomer = Monomer(mol)
        assert monomer.n_atoms == n_atoms
        assert monomer.n_bonds == n_bonds
        assert monomer.n_tags == 0
        assert monomer.name == "UNK"

    def test_add_tag(self):
        mol = pdb2rdmol(m1c1)
        monomer = Monomer(mol)
        monomer.add_tag("left", [43, 46, 44, 45])
        assert monomer.n_tags == 1


    