#include <algorithm>
#include <set>
#include <string>
#include <iostream>
#include <utility>

#ifndef _PT_MONOMER_C
#define _PT_MONOMER_C

#include <GraphMol/GraphMol.h>
// #include "atom.h"
#include "monomer.h"

namespace polytop {
    

    int Monomer::addAtom(Atom &atom) {
        atoms.emplace_back(&atom);
        return atoms.size();
    }

    Monomer::Monomer(const RDKit::ROMol &mol, std::string name) {
        setRDMol(RDKit::RWMol(mol));
        resName = name;
    }
    


};

#endif