

#include <array>
#include <set>
#include <algorithm>
#include <iostream>


#ifndef _PT_PARAM_H
#define _PT_PARAM_H

#include "atom.hpp"
namespace polytop {


    enum PARAMS {
        bond,
        angle,
        dihedral,
        improper,
        pair,
        exclusion
    };

    template <std::size_t numAtoms, PARAMS pType> class Param {

        public:
            Param(std::array< Atom*, numAtoms> indices);
            Param(Param &param);
            std::array<Atom*, numAtoms> atoms;
            std::size_t nAtoms = numAtoms;
            std::array<unsigned int, numAtoms> getAtomIndices();
            std::array<Atom*, numAtoms> getOrderedAtomIndices();
            std::array<Atom*, numAtoms> updateAtom(Atom *oldAtom, Atom *newAtom);
            // Param<numAtoms, pType> copyToMol(MonomerUnit &mol);
            bool containsAtom(Atom *atom);
            bool withinAtoms(std::set<Atom*> atomSet);
            bool withinAtoms(std::vector<Atom*> atomVector);
            PARAMS paramType = pType;
    };


    template <std::size_t numAtoms, PARAMS pType>
    Param<numAtoms, pType>::Param(std::array<Atom*, numAtoms> atoms) : atoms(atoms) {}

    template <std::size_t numAtoms, PARAMS pType>
    std::array<unsigned int, numAtoms> Param<numAtoms, pType>::getAtomIndices() {
        std::array<unsigned int, numAtoms> indices;
        for (std::size_t i = 0; i < numAtoms; i++) {
            indices[i] = atoms.at(i)->index;
        }
        return indices;
    };

    template <std::size_t numAtoms, PARAMS pType>
    std::array<Atom*, numAtoms> Param<numAtoms, pType>::getOrderedAtomIndices() {
        std::array<Atom*, numAtoms> indices = getAtomIndices();
        if (indices.back() < indices.front()) {
            std::reverse(indices.begin(), indices.end());
        }
        return atoms;
    };

    template <std::size_t numAtoms, PARAMS pType>
    std::array<Atom*, numAtoms> Param<numAtoms, pType>::updateAtom(Atom *oldAtom, 
                                                                   Atom *newAtom) {
        std::replace(atoms.begin(), atoms.end(), oldAtom, newAtom);
        return atoms;
    };

    // template <std::size_t numAtoms, PARAMS pType>
    // Param<numAtoms, pType> Param<numAtoms, pType>::copyToMol(MonomerUnit &mol) {
    //     std::array<unsigned int, numAtoms> indices = getAtomIndices();
    //     Param<numAtoms, pType> param = Param<numAtoms, pType>(mol, indices);
    //     return param;
    // }

    template <std::size_t numAtoms, PARAMS pType>
    bool Param<numAtoms, pType>::containsAtom(Atom *atom) {
        auto found = std::find(atoms.begin(), atoms.end(), atom);
        return found != atoms.end();
    }

    template <std::size_t numAtoms, PARAMS pType>
    bool Param<numAtoms, pType>::withinAtoms(std::set<Atom*> atomSet) {
        for (auto &el : atoms) {
            if (atomSet.find(el) == atomSet.end()) return false;
        }
        return true;
    }

    template <std::size_t numAtoms, PARAMS pType>
    bool Param<numAtoms, pType>::withinAtoms(std::vector<Atom*> atomVector) {
        for (auto &el : atoms) {
            if (std::find(atomVector.begin(), atomVector.end(), el) == atomVector.end()) {
                return false;
            }
        }
        return true;
    }


    typedef Param<2, bond> Bond;
    typedef Param<2, pair> Pair;
    typedef Param<2, exclusion> Exclusion;
    typedef Param<3, angle> Angle;
    typedef Param<4, dihedral> Dihedral;
    typedef Param<4, improper> Improper;




    class HarmonicBond : Bond {
        public:
            double eqValue = 0;
            double forceConstant = 0;
        
        // private:
        //     int gmxFInt = 1;
    };

}

#endif