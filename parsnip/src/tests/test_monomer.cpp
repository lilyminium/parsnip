

#include <vector>

#include "catch2.hpp"
#include "base.hpp"
#include "monomer.hpp"

using namespace polytop;

TEST_CASE( "Monomer-specific functions", "[monomer]" ) {

    RDKit::RWMol *mol = PDBtoRDMol("m1_u4_c001_opt.pdb");
    Monomer mon = Monomer(*mol, "M1");

    REQUIRE( mon.atoms.size() == 54 );
    REQUIRE( mon.name == "M1" );
    REQUIRE( mon.tags->size() == 0 );
    REQUIRE( mon.bonds->size() == 53 );
    REQUIRE( mon.angles->size() == 0 );
    REQUIRE( mon.dihedrals->size() == 0 );
    REQUIRE( mon.impropers->size() == 0 );
    REQUIRE( mon.pairs->size() == 0 );
    REQUIRE( mon.exclusions->size() == 0 );

    REQUIRE( mon.atoms[0]->getResName() == "M1");
    REQUIRE( mon.atoms[0]->getSerial() == 1);


    SECTION("Get atom indices") {
        mon.atoms[0]->index = 3;
        auto indices = mon.getAtomIndices();
        REQUIRE( indices.size() == 54 );
        REQUIRE( indices[0] == 3 );
        REQUIRE( indices[20] == 20 );
    };

    SECTION("Get atoms by indices") {
        auto atoms = mon.getAtomsByIndices({1, 2});
        REQUIRE( atoms.size() == 2 );
        REQUIRE( atoms[0]->index == 1 );
        REQUIRE( atoms[1]->index == 2 );
    };

    SECTION("Delete atom") {
        REQUIRE( mon.atoms.size() == 54 );
        REQUIRE( mon.atoms[0]->name == "C" );
        mon.delAtomFromVector(mon.atoms[0]);
        REQUIRE( mon.atoms.size() == 53 );
        REQUIRE( mon.atoms[0]->name == "H" );
    };

    SECTION("Delete atom by index") {
        REQUIRE( mon.atoms.size() == 54 );
        REQUIRE( mon.rdMol.getNumAtoms() == 54);
        REQUIRE( mon.atoms[0]->name == "C" );
        mon.delAtomByIndex(0);
        REQUIRE( mon.atoms.size() == 53 );
        REQUIRE( mon.rdMol.getNumAtoms() == 54);
        REQUIRE( mon.atoms[0]->name == "H" );
    };

    SECTION("Delete atoms") {
        REQUIRE( mon.atoms.size() == 54 );
        AtomPVec atoms = mon.getAtomsByIndices({1, 2, 3});
        mon.delAtomsFromVector(atoms);
        REQUIRE( mon.atoms.size() == 51 );
    };

    SECTION("Remove atom") {
        REQUIRE( mon.atoms.size() == 54 );
        REQUIRE( mon.rdMol.getNumAtoms() == 54 );
        auto newSize = mon.removeAtom(mon.atoms[0]);

        REQUIRE( newSize == 53 );
        REQUIRE( mon.bonds->size() == 49 );
        REQUIRE( mon.rdMol.getNumAtoms() == 53 );
        REQUIRE( mon.rdMol.getNumBonds() == 49 );
        REQUIRE( mon.atoms[0]->index == 0 ); //reindexes
    };

    SECTION("Remove atom by index") {
        REQUIRE( mon.atoms.size() == 54 );
        REQUIRE( mon.rdMol.getNumAtoms() == 54 );
        auto newSize = mon.removeAtomByIndex(0);

        REQUIRE( newSize == 53 );
        REQUIRE( mon.bonds->size() == 49 );
        REQUIRE( mon.rdMol.getNumAtoms() == 53 );
        REQUIRE( mon.rdMol.getNumBonds() == 49 );
        REQUIRE( mon.atoms[0]->index == 0 ); //reindexes
    };

    SECTION("Contains atom") {
        REQUIRE( mon.containsAtom(mon.atoms[0]) == true );
    };

    SECTION("Does not contain atom") {
        auto atom = mon.atoms[0];
        mon.removeAtom(atom);
        REQUIRE( mon.containsAtom(atom) == false );
    };

    SECTION("Contains all atoms") {
        AtomPVec atoms = mon.getAtomsByIndices({1, 2, 3});
        REQUIRE( mon.containsAllAtoms(atoms) == true );
    };

    SECTION("Does not contain all atoms") {
        AtomPVec atoms = mon.getAtomsByIndices({1, 2, 3});
        mon.removeAtomByIndex(1);
        REQUIRE( mon.containsAllAtoms(atoms) == false );
    };
    
    SECTION("Replacing atoms removes the old atom") {
        auto c = mon.atoms[4];
        auto h = mon.atoms[1];
        mon.removeAtom(c);

        BondVector* bonds = mon.bonds;
        Bond* bond = (*bonds)[0];
        std::array<unsigned int, 2> indices = {1, 0};

        REQUIRE ( mon.atoms.size() == 53 );
        REQUIRE ( mon.bonds->size() == 49 );
        REQUIRE ( bond->getAtomIndices() == indices );

        mon.replaceAtom(h, c);
        // it does remove the old atom and does not
        // insert the new one (as this method is mostly for
        // polymers)
        std::array<unsigned int, 2> indices2 = {4, 0};
        REQUIRE ( mon.atoms.size() == 52 );
        REQUIRE ( mon.bonds->size() == 48 );
        REQUIRE ( bond->getAtomIndices() == indices2 );
    }

    SECTION("Clean parameters") {
        mon.delAtomByIndex(0);
        REQUIRE( mon.bonds->size() == 53 );
        mon.cleanParams();
        REQUIRE( mon.bonds->size() == 49 );
    }

};

TEST_CASE("Tag functions", "[monomer][tags]") {

    RDKit::RWMol *mol = PDBtoRDMol("m1_u4_c001_opt.pdb");
    Monomer mon = Monomer(*mol, "M1");
    REQUIRE( mon.tags->size() == 0 );

    SECTION("Add tag from numbers") {
        REQUIRE( mon.tags->size() == 0 );
        mon.addTag("left", {43, 46, 44, 45});
        REQUIRE( mon.tags->size() == 1 );
    };

    SECTION("Add tag from atoms") {
        REQUIRE( mon.tags->size() == 0 );
        mon.addTag("test", {mon.atoms[0]});
        REQUIRE( mon.tags->size() == 1 );
    };

    SECTION("Add created tag") {
        REQUIRE( mon.tags->size() == 0 );
        Tag* tag = new Tag("test", {mon.atoms[0]});
        REQUIRE( mon.addTag(tag) == 1 );

        delete tag;
    };

    SECTION("Safe add tag only if atoms exist") {
        REQUIRE( mon.tags->size() == 0 );
        Tag* tag = new Tag("test", {mon.atoms[0]});
        mon.removeAtomByIndex(0);
        REQUIRE( mon.safeAddTag(tag) == 0 );

        delete tag;
    };

    SECTION("Clean tags") {
        auto size = mon.addTag("left", {43, 46, 44, 45});
        REQUIRE( size == 1 );
        mon.removeAtomByIndex(43);
        REQUIRE ( mon.cleanTags() == 0 );
    };

    SECTION("Get first/last tag by name") {
        Tag* tag1 = new Tag("test", {mon.atoms[0]});
        Tag* tag2 = new Tag("test", {mon.atoms[2]});

        mon.addTag(tag1);
        mon.addTag(tag2);

        REQUIRE( mon.tags->size() == 2 );
        REQUIRE( mon.getFirstTagByName("test") == tag1 );
        REQUIRE( mon.getLastTagByName("test") == tag2 );

        delete tag1;
        delete tag2;
    }
};

TEST_CASE("Create monomer units", "[monomer][monomerunit]") {
    RDKit::RWMol *mol = PDBtoRDMol("m1_u4_c001_opt.pdb");
    Monomer *mon = new Monomer(*mol, "M1");
    mon->addTag("left", {43, 46, 44, 45});

    MonomerUnit unit = MonomerUnit(mon);

    REQUIRE( unit.name == "M1" );
    REQUIRE( unit.atoms.size() == 54 );
    REQUIRE( unit.bonds->size() == 53 );
    REQUIRE( unit.angles->size() == 0 );
    REQUIRE( unit.atoms[0] != mon->atoms[0] );
    REQUIRE( unit.atoms[0]->name == mon->atoms[0]->name );
    REQUIRE( mon->atoms[0]->getResName() == "M1" );
    REQUIRE( unit.atoms[0]->getResName() == "M1" );
    REQUIRE( unit.tags->size() == 1);

    delete mon;
}

    