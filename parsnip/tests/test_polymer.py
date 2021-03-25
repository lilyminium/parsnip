
import pytest

from parsnip.monomer import Monomer
from parsnip.polymer import Polymer

from parsnip.tests.datafiles import m1c1, m1c2
from .utils import pdb2rdmol

class TestPolymer:

    @pytest.fixture()
    def m1(self):
        mol = pdb2rdmol(m1c1)
        monomer = Monomer(mol, name="M1")
        monomer.add_tag("left", [43, 46, 44, 45])
        monomer.add_tag("right", [48, 51, 53, 52])
        return monomer
    
    @pytest.fixture()
    def m2(self):
        mol = pdb2rdmol(m1c2)
        monomer = Monomer(mol, name="M2")
        monomer.add_tag("left", [43, 46, 44, 45])
        monomer.add_tag("right", [48, 51, 53, 52])
        return monomer


    def test_create_polymer(self):
        polymer = Polymer()
        assert polymer.n_atoms == 0
        assert polymer.n_bonds == 0
        assert polymer.n_tags == 0
        assert polymer.name == "UNK"
        assert len(polymer._universe.atoms) == 0

    def test_add_unit(self, m1):
        polymer = Polymer()
        polymer.add_monomer(m1)

        assert polymer.n_atoms == 54
        assert polymer.n_bonds == 53
        assert polymer.n_tags == 2
        assert len(polymer._universe.atoms) == 54

    def test_add_units_via_tags(self, m1):
        polymer = Polymer()
        polymer.add_monomer(m1)
        polymer.add_monomer(m1, monomer_tag="left", polymer_tag="right",
                            replace_polymer_atoms=[[True, False, True, True]])
        assert polymer.n_atoms == 104
        assert polymer.n_bonds == 104
        assert polymer.n_tags == 2
        assert len(polymer._universe.atoms) == 104

    def test_add_multiple_units(self, m1, m2):
        polymer = Polymer()
        polymer.add_monomers(monomers=[m1, m2], monomer_tags="left",
                             polymer_tags="right", replace_polymer_atoms=True,
                             ratio=[0.3, 0.7], total=100)
        assert polymer.n_units == 100