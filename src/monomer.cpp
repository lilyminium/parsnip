#include <algorithm>
#include <set>
#include <string>

#ifndef _PT_MONOMER_C
#define _PT_MONOMER_C

// #include "atom.h"
#include "monomer.h"

namespace polytop {

    // Monomer::Monomer() {};
    Monomer::Monomer(std::string name) {
        resName = name;
    };

    int Monomer::delAtom(Atom &atom) {
        auto found = std::find(atoms.begin(), atoms.end(), &atom);
        if (found != atoms.end()) {
            // int index = std::distance(atoms.begin(), found);
            atoms.erase(found);
        }

        return atoms.size();
    }


    int Monomer::addAtom(Atom &atom) {
        atoms.push_back(&atom);
        return atoms.size();
    }

    MonomerUnit::MonomerUnit(Monomer mol) {
        
    }

    



    // void Monomer::removeParametersInAtomSubset(std::set<int> atomIndices) {
    //     for (auto& bond : bonds) {

    //     }
    // };

    // MonomerUnit Monomer::createUnit() {
    //     MonomerUnit unit;
    //     // atoms
    //     unit.atoms.reserve(atoms.size());
    //     for (size_t i = 0; i < atoms.size(); i++) {
    //         unit.atoms[i] = atoms[i].copy();
    //         unit.atoms[i].owningMol = unit;
    //     };

    //     // parameters
    //     for (auto& el : atoms) {
    //         el.copyParamsTo(unit);
    //     };

    //     return unit;
    // };

};

#endif