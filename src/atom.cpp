#include <array>

#ifndef _PT_ATOM_C
#define _PT_ATOM_C

#include "atom.h"
#include "monomer.h"

namespace polytop {

    Atom::Atom() { owningMol = new Monomer(); };

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

    std::unique_ptr<Atom> Atom::copy() {
        std::unique_ptr<Atom> atom(new Atom);
        atom->index = index;
        atom->name = name;
        return atom;
    };

    void Atom::copyToMol(MonomerUnit mol) {
        std::unique_ptr<Atom> atom = copy();
        atom->owningMol = &mol;
        mol.atoms.push_back(atom);
    }

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
    
    void MonomerUnit::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
        for (auto &atom : atoms) { atom->removeParamsWithinAtomSet(atomSet); }
    }
    

};

#endif