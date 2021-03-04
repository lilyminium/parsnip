

#include <vector>

#include "catch2.hpp"

#include "base.hpp"
#include "utils.hpp"

using namespace polytop;

TEST_CASE( "Neighbor ints are obtained", "[utils]" ) {
    RDKit::RWMol *mol = PDBtoRDMol("m1_u4_c001_opt.pdb");

    REQUIRE( getNeighborInts(*mol, {0}, 0) == IntVec({}) );
    REQUIRE( getNeighborInts(*mol, {0}, 1) == IntVec({1, 2, 3, 4}) );
    REQUIRE( getNeighborInts(*mol, {0}, 2) == IntVec({1, 2, 3, 4, 5, 8, 47}) );
    REQUIRE( getNeighborInts(*mol, {0}, 3) == IntVec({1, 2, 3, 4, 5, 8, 47,
                                                      6, 7, 9, 10, 43, 48, 49, 50}) );

    REQUIRE( getNeighborInts(*mol, {14, 17, 18}, 0) == IntVec({}) );
    REQUIRE( getNeighborInts(*mol, {14, 17, 18}, 1) == IntVec({11, 15, 16, 19, 20, 21}) );
    REQUIRE( getNeighborInts(*mol, {14, 17, 18}, 2) == IntVec({11, 15, 16, 10, 12, 13,
                                                               19, 20, 21, 22, 23, 24}) );
};

TEST_CASE( "get RDMol subset", "[utils]" ) {
    RDKit::RWMol *mol = PDBtoRDMol("m1_u4_c002_opt.pdb");

    SECTION("subset 0") {
        auto submol = subsetRDMol(*mol, {0});
        REQUIRE( submol->getNumAtoms() == 1 );
        REQUIRE( submol->getAtomWithIdx(0)->getSymbol() == "C" );
    };

    SECTION("subset 4 atoms") {
        auto submol = subsetRDMol(*mol, {14, 17, 18, 11, 15, 16, 19, 20, 21});

        RDKit::MolToPDBFile(*submol, "/Users/lily/Desktop/capped.pdb");
        REQUIRE( submol->getNumAtoms() == 9 );
        REQUIRE( submol->getAtomWithIdx(0)->getSymbol() == "C" );
        REQUIRE( submol->getAtomWithIdx(1)->getSymbol() == "O" );
        REQUIRE( submol->getAtomWithIdx(5)->getSymbol() == "H" );

        REQUIRE( submol->getNumBonds() == 8 );
        REQUIRE( submol->getBondWithIdx(0)->getBeginAtomIdx() == 0 );
        REQUIRE( submol->getBondWithIdx(0)->getEndAtomIdx() == 1 );
    }

};

TEST_CASE( "split integer into ratios", "[utils]") {
    REQUIRE( splitIntegerIntoRatio(5, {0.2, 0.3}) == IntVec({2, 3}) );
    REQUIRE( splitIntegerIntoRatio(10, {4, 6}) == IntVec({4, 6}) );
    REQUIRE( splitIntegerIntoRatio(9, {1, 2, 1}) == IntVec({3, 4, 2}) );
    REQUIRE( splitIntegerIntoRatio(28, {0.21, 0.79}) == IntVec({6, 22}) );
}


