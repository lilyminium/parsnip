#include <array>
#include <iostream>

#include <boost/range/iterator_range.hpp>
#ifndef _PT_ATOM_C
#define _PT_ATOM_C


#include "atom.hpp"
#include <GraphMol/MolOps.h>

namespace polytop {

    Atom::Atom(RDKit::Atom *atom, std::string name, double charge) : rdAtom(atom), name(name), charge(charge) {
        setRDAtom(atom);
    };

    void Atom::setRDAtom(RDKit::Atom *atom) {
        index = atom->getIdx();
        if (name == "X") {name = atom->getSymbol();};
    };

};

#endif