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

        atom->setMonomerInfo(resInfo);
        resInfo->setSerialNumber(index + 1);
        resInfo->setName(name);
    };

    void Atom::setResName(std::string resName) {
        resInfo->setResidueName(resName);
    };
    std::string Atom::getResName() {
        return resInfo->getResidueName();
    };

    void Atom::setResNum(unsigned int resNum) {
        resInfo->setResidueNumber(resNum);
    };

    int Atom::getResNum() {
        return resInfo->getResidueNumber();
    }

    void Atom::setIndex(unsigned int atomIndex) {
        index = atomIndex;
        resInfo->setSerialNumber(index + 1);
    }

    unsigned int Atom::getIndex() {
        return index;
    }

    unsigned int Atom::getSerial() {
        return resInfo->getSerialNumber();
    }

    void Atom::updateMonomerInfo() {
        auto newInfo = new RDKit::AtomPDBResidueInfo(*resInfo);
        rdAtom->setMonomerInfo(newInfo);
    }



};

#endif