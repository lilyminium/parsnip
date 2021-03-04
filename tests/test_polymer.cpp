

#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include "catch2.hpp"
#include "base.hpp"
#include "polymer.hpp"

using namespace polytop;

TEST_CASE( "Creating a polymer", "[polymer]" ) {

    RDKit::RWMol *mol = PDBtoRDMol("m1_u4_c001_opt.pdb");
    Monomer* mon = new Monomer(*mol, "M1");
    mon->addTag("left", {43, 46, 44, 45});
    mon->addTag("right", {48, 51, 53, 52});

    Polymer pol = Polymer("PolyAla");

    REQUIRE( pol.atoms.size() == 0 );
    REQUIRE( pol.units.size() == 0 );
    REQUIRE( pol.name == "PolyAla" );
    REQUIRE( pol.tags->size() == 0 );
    REQUIRE( pol.bonds->size() == 0 );
    REQUIRE( pol.angles->size() == 0 );
    REQUIRE( pol.dihedrals->size() == 0 );
    REQUIRE( pol.impropers->size() == 0 );
    REQUIRE( pol.pairs->size() == 0 );
    REQUIRE( pol.exclusions->size() == 0 );

    SECTION("Add monomer unit to empty") {
        MonomerUnit *unit = new MonomerUnit(mon);
        pol.addMonomerUnit(unit);
        REQUIRE( pol.atoms.size() == 54 );
        REQUIRE( pol.units.size() == 1 );
        REQUIRE( pol.tags->size() == 2 );
        REQUIRE( pol.bonds->size() == 53 );
        REQUIRE( pol.angles->size() == 0 );

        REQUIRE( pol.units[0]->atoms.size() == 54 );
        REQUIRE( pol.units[0]->bonds->size() == 53 );

        delete unit;
        
    };

    SECTION("Add monomer to empty") {
        pol.addMonomer(mon);
        REQUIRE( pol.atoms.size() == 54 );
        REQUIRE( pol.units.size() == 1 );
        REQUIRE( pol.tags->size() == 2 );
        REQUIRE( pol.bonds->size() == 53 );
        REQUIRE( pol.angles->size() == 0 );

        REQUIRE( pol.units[0]->atoms.size() == 54 );
        REQUIRE( pol.units[0]->bonds->size() == 53 );

        SECTION("Add second monomer") {
            pol.addMonomer(mon, "left", "right", true);
            REQUIRE( pol.atoms.size() == 104 );
            REQUIRE( pol.units.size() == 2 );
            REQUIRE( pol.tags->size() == 3 );
            REQUIRE( pol.bonds->size() == 103 );
            REQUIRE( pol.angles->size() == 0 );

            REQUIRE( pol.units[0]->atoms.size() == 50 );
            REQUIRE( pol.units[0]->bonds->size() == 49 );
            REQUIRE( pol.units[1]->atoms.size() == 54 );
            REQUIRE( pol.units[1]->bonds->size() == 53 );

            SECTION("Add third monomer") {
                pol.addMonomer(mon, "left", "right", true);
                REQUIRE( pol.atoms.size() == 154 );
                REQUIRE( pol.units.size() == 3 );
                REQUIRE( pol.tags->size() == 4 );
                REQUIRE( pol.bonds->size() == 153 );
                REQUIRE( pol.angles->size() == 0 );

                REQUIRE( pol.units[0]->atoms.size() == 50 );
                REQUIRE( pol.units[0]->bonds->size() == 49 );
                REQUIRE( pol.units[1]->atoms.size() == 50 );
                REQUIRE( pol.units[1]->bonds->size() == 49 );
                REQUIRE( pol.units[2]->atoms.size() == 54 );
                REQUIRE( pol.units[2]->bonds->size() == 53 );

            }
        };
        
    };

    

    
};

