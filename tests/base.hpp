#include <string>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>


#ifndef _TEST_BASE_H
#define _TEST_BASE_H

#include "polymer.hpp"

namespace polytop {

    std::string ptbase = "/Users/lily/pydev/colytop/"; //getenv("POLYBASE");
    std::string testdata = ptbase + "tests/data/";

    typedef std::vector<unsigned int> IntVec;
    typedef std::vector<Atom*> AtomPVec;



    RDKit::RWMol* PDBtoRDMol(std::string filename) {
        std::string fName = testdata + "pdb/" + filename;
        return RDKit::PDBFileToMol(fName, false, false);
    };
};

#endif