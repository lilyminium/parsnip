

#include <vector>

#include "catch2.hpp"
#include "base.hpp"
#include "atom.hpp"

using namespace polytop;

TEST_CASE( "Create Atom", "[atom]" ) {

    RDKit::RWMol *mol = PDBtoRDMol("m1_u4_c001_opt.pdb");

    SECTION("atom C") {
        Atom *atom = new Atom(mol->getAtomWithIdx(0));
        
        REQUIRE( atom->index == 0 );
        REQUIRE( atom->name == "C" );
        REQUIRE( atom->charge == 0 );
    };

    SECTION("atom H") {
        Atom *atom = new Atom(mol->getAtomWithIdx(1));
        
        REQUIRE( atom->index == 1 );
        REQUIRE( atom->name == "H" );
        REQUIRE( atom->charge == 0 );
    };

    SECTION("Override rdAtom name") {
        Atom *atom = new Atom(mol->getAtomWithIdx(0), "Y");
        
        REQUIRE( atom->index == 0 );
        REQUIRE( atom->name == "Y" );
        REQUIRE( atom->charge == 0 );
    };

    
    

};
