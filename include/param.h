

#include <array>
#include <set>
#include <algorithm>
#include <iostream>


#ifndef _PT_PARAM_H
#define _PT_PARAM_H

#include "monomer.h"
#include "polymer.h"

namespace polytop {
    class Atom;
    // class Monomer;
    // class MonomerUnit;

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
            Param(Monomer &mol, std::array<unsigned int, numAtoms> indices);
            Param(MonomerUnit &mol, std::array<unsigned int, numAtoms> indices);
            Param(Monomer &mol, Param &param);
            Param(MonomerUnit &mol, Param &param);
            std::array<Atom*, numAtoms> atoms;
            std::array<unsigned int, numAtoms> getAtomIndices();
            std::array<Atom*, numAtoms> getOrderedAtomIndices();
            std::array<Atom*, numAtoms> updateAtom(Atom *oldAtom, Atom *newAtom);
            Param<numAtoms, pType> copyToMol(MonomerUnit &mol);
            bool containsAtom(Atom *atom);
            bool withinAtomSet(std::set<Atom*> atomSet);
            PARAMS paramType = pType;
            // std::string getParamName();
    };

    template <std::size_t numAtoms, PARAMS pType>
    Param<numAtoms, pType>::Param(Monomer &mol, std::array<unsigned int, numAtoms> indices) {
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i] = mol.atoms.at(indices.at(i));
        };
    };

    template <std::size_t numAtoms, PARAMS pType>
    Param<numAtoms, pType>::Param(MonomerUnit &mol, std::array<unsigned int, numAtoms> indices) {
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i] = mol.atoms.at(indices.at(i));
        };
    };

    template <std::size_t numAtoms, PARAMS pType>
    Param<numAtoms, pType>::Param(Monomer &mol, Param<numAtoms, pType> &param) {
        auto indices = param.getAtomIndices();
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i] = mol.atoms.at(indices.at(i));
        };
    };

    template <std::size_t numAtoms, PARAMS pType>
    Param<numAtoms, pType>::Param(MonomerUnit &mol, Param<numAtoms, pType> &param) {
        auto indices = param.getAtomIndices();
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i] = mol.atoms.at(indices.at(i));
        };
    };

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

    template <std::size_t numAtoms, PARAMS pType>
    Param<numAtoms, pType> Param<numAtoms, pType>::copyToMol(MonomerUnit &mol) {
        std::array<unsigned int, numAtoms> indices = getAtomIndices();
        Param<numAtoms, pType> param = Param<numAtoms, pType>(mol, indices);
        return param;
    }

    template <std::size_t numAtoms, PARAMS pType>
    bool Param<numAtoms, pType>::containsAtom(Atom *atom) {
        auto found = std::find(atoms.begin(), atoms.end(), atom);
        return found != atoms.end();
    }

    template <std::size_t numAtoms, PARAMS pType>
    bool Param<numAtoms, pType>::withinAtomSet(std::set<Atom*> atomSet) {
        for (auto &el : atoms) {
            if (atomSet.find(el) == atomSet.end()) return false;
        }
        return true;
    }

    // typedef Param<2> TwoAtomParam;
    // typedef Param<3> ThreeAtomParam;
    // typedef Param<4> FourAtomParam;

    typedef Param<2, bond> Bond;
    typedef Param<2, pair> Pair;
    typedef Param<2, exclusion> Exclusion;
    typedef Param<3, angle> Angle;
    typedef Param<4, dihedral> Dihedral;
    typedef Param<4, improper> Improper;


    // class BondVector;
    // class PairVector;
    // class ExclusionVector;
    // class AngleVector;
    // class DihedralVector;
    // class ImproperVector;

    // class Bond : public TwoAtomParam {
    // };
    // class Pair : public TwoAtomParam {
    // };
    // class Exclusion : public TwoAtomParam {
    // };
    // class Angle : public ThreeAtomParam {
    // };
    // class Dihedral : public FourAtomParam {
    // };
    // class Improper : public FourAtomParam {
    // };



    class HarmonicBond : Bond {
        public:
            double eqValue = 0;
            double forceConstant = 0;
        
        // private:
        //     int gmxFInt = 1;
    };

}

#endif