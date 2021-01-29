#include <array>

#ifndef _PT_ATOM_C
#define _PT_ATOM_C


#include "atom.h"

namespace polytop {

    Atom::Atom() {};

    Atom::Atom(const Atom &rhs) {
        index = rhs.index;
        name = rhs.name;
    }

    Atom::Atom(Monomer &mol) { owningMol = &mol; };

    Atom::Atom(MonomerUnit &mol) { owningMol = &mol; };

    void Atom::replaceWithAtom(Atom newAtom) {
        bonds.updateAtom(newAtom);
        pairs.updateAtom(newAtom);
        exclusions.updateAtom(newAtom);
        angles.updateAtom(newAtom);
        dihedrals.updateAtom(newAtom);
        impropers.updateAtom(newAtom);
        newAtom.owningMol->delAtom(*this);
    };

    // // void Atom::copyToMol(MonomerUnit mol) {
    // //     std::shared_ptr<Atom> atom = copy();
    // //     atom->owningMol = &mol;
    // //     mol.atoms.push_back(&atom);
    // // }

    void Atom::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
        bonds.removeParamsWithinAtomSet(atomSet);
        pairs.removeParamsWithinAtomSet(atomSet);
        exclusions.removeParamsWithinAtomSet(atomSet);
        angles.removeParamsWithinAtomSet(atomSet);
        dihedrals.removeParamsWithinAtomSet(atomSet);
        impropers.removeParamsWithinAtomSet(atomSet);
    }

    template<>
    BondVector BondVector::getCorrespondingParamVector(Atom *atom) {return atom->bonds;}
    template<>
    PairVector PairVector::getCorrespondingParamVector(Atom *atom) {return atom->pairs;}
    template<>
    ExclusionVector ExclusionVector::getCorrespondingParamVector(Atom *atom) {return atom->exclusions;}
    template<>
    AngleVector AngleVector::getCorrespondingParamVector(Atom *atom) {return atom->angles;}
    template<>
    DihedralVector DihedralVector::getCorrespondingParamVector(Atom *atom) {return atom->dihedrals;}
    template<>
    ImproperVector ImproperVector::getCorrespondingParamVector(Atom *atom) {return atom->impropers;}
    

    // BondVector Bond::getCorrespondingVector(Atom *atom) {return atom->bonds;}
    // PairVector Pair::getCorrespondingVector(Atom *atom) {return atom->pairs;}
    // ExclusionVector Exclusion::getCorrespondingVector(Atom *atom) {return atom->exclusions;}
    // AngleVector Angle::getCorrespondingVector(Atom *atom) {return atom->angles;}
    // DihedralVector Dihedral::getCorrespondingVector(Atom *atom) {return atom->dihedrals;}
    // ImproperVector Improper::getCorrespondingVector(Atom *atom) {return atom->impropers;}

    // void MonomerUnit::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
    //     for (auto &atom : atoms) { atom->removeParamsWithinAtomSet(atomSet); }
    // }

    int Monomer::delAtom(Atom &atom) {
        auto found = std::find(atoms.begin(), atoms.end(), &atom);
        if (found != atoms.end()) {
            // int index = std::distance(atoms.begin(), found);
            atoms.erase(found);
        }

        return atoms.size();
    };
    
    MonomerUnit::MonomerUnit(Monomer mol) {
        // atoms
        atoms.reserve(mol.atoms.size());
        for (auto &atom : mol.atoms) {
            Atom *newAtom = new Atom(*atom);
            atoms.push_back(newAtom);
            newAtom->owningMol = this;
        }

        // params
        copyParams<Bond, BondVector>(mol);
        copyParams<Pair, PairVector>(mol);
    };



};

#endif