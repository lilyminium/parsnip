

#include <array>
#include <set>
#include <algorithm>

#ifndef _PT_PARAM_H
#define _PT_PARAM_H

#include "monomer.h"

namespace polytop {
    class Atom;
    class Monomer;
    class MonomerUnit;

    template <std::size_t numAtoms> class Param {

        public:
            Param(Monomer &mol, std::array<int, numAtoms> indices);
            Param(MonomerUnit &mol, std::array<int, numAtoms> indices);
            std::array<Atom*, numAtoms> atoms;
            std::array<int, numAtoms> getAtomIndices();
            std::array<Atom*, numAtoms> getOrderedAtomIndices();
            std::array<Atom*, numAtoms> updateAtom(Atom *oldAtom, Atom *newAtom);
            Param<numAtoms> copyToMol(MonomerUnit &mol);
            bool containsAtom(Atom *atom);
            bool withinAtomSet(std::set<Atom*> atomSet);
            // std::string getParamName();
    };

    template <std::size_t numAtoms>
    Param<numAtoms>::Param(Monomer &mol, std::array<int, numAtoms> indices) {
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i] = mol.atoms.at(indices.at(i));
        };
    };

    template <std::size_t numAtoms>
    Param<numAtoms>::Param(MonomerUnit &mol, std::array<int, numAtoms> indices) {
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i] = mol.atoms.at(indices.at(i));
        };
    };

    template <std::size_t numAtoms>
    std::array<int, numAtoms> Param<numAtoms>::getAtomIndices() {
        std::array<int, numAtoms> indices;
        for (std::size_t i = 0; i < numAtoms; i++) {
            indices[i] = atoms.at(i)->index;
        }
        return indices;
    };

    template <std::size_t numAtoms>
    std::array<Atom*, numAtoms> Param<numAtoms>::getOrderedAtomIndices() {
        std::array<Atom*, numAtoms> indices = getAtomIndices();
        if (indices.back() < indices.front()) {
            std::reverse(indices.begin(), indices.end());
        }
        return atoms;
    };

    template <std::size_t numAtoms>
    std::array<Atom*, numAtoms> Param<numAtoms>::updateAtom(Atom *oldAtom, 
                                                            Atom *newAtom) {
        std::replace(atoms.begin(), atoms.end(), oldAtom, newAtom);
        return atoms;
    };

    template <std::size_t numAtoms>
    Param<numAtoms> Param<numAtoms>::copyToMol(MonomerUnit &mol) {
        std::array<int, numAtoms> indices = getAtomIndices();
        Param<numAtoms> param = Param<numAtoms>(mol, indices);
        return param;
    }

    template <std::size_t numAtoms>
    bool Param<numAtoms>::containsAtom(Atom *atom) {
        auto found = std::find(atoms.begin(), atoms.end(), atom);
        return found != atoms.end();
    }

    template <std::size_t numAtoms>
    bool Param<numAtoms>::withinAtomSet(std::set<Atom*> atomSet) {
        for (auto &el : atoms) {
            if (atomSet.find(el) == atomSet.end()) return false;
        }
        return true;
    }

    typedef Param<2> TwoAtomParam;
    typedef Param<3> ThreeAtomParam;
    typedef Param<4> FourAtomParam;

    class Bond : public TwoAtomParam {};
    class Pair : public TwoAtomParam {};
    class Exclusion : public TwoAtomParam {};
    class Angle : public ThreeAtomParam {};
    class Dihedral : public FourAtomParam {};
    class Improper : public FourAtomParam {};



    class HarmonicBond : Bond {
        public:
            double eqValue = 0;
            double forceConstant = 0;
        
        private:
            int gmxFInt = 1;
    };

}

#endif