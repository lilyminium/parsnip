#include <array>
#include <iostream>

#include <boost/range/iterator_range.hpp>
#ifndef _PT_ATOM_C
#define _PT_ATOM_C


#include "atom.h"
#include <GraphMol/MolOps.h>

namespace polytop {

    // Atom::Atom() {};

    // Atom::Atom(const Atom &rhs) {
    //     index = rhs.index;
    //     name = rhs.name;
    // }

    // Atom::Atom(MonomerUnit &mol) {
    //     owningMol = &mol;
    //     int idx = mol.rdMol.addAtom(true);
    //     RDKit::Atom *atom = mol.rdMol.getAtomWithIdx(idx);
    //     setRDAtom(*atom, true);
    // };

    Atom::Atom(Monomer &mol, int atomIndex) {
        owningMol = &mol;
        RDKit::Atom *atom = mol.rdMol.getAtomWithIdx(atomIndex);
        setRDAtom(atom, true);
    };

    Atom::Atom(MonomerUnit &mol, int atomIndex) {
        owningMol = &mol;
        RDKit::Atom *atom = mol.rdMol.getAtomWithIdx(atomIndex);
        setRDAtom(atom, true);
    };

    void Atom::setRDAtom(RDKit::Atom *atom, bool updateCharge) {
        name = atom->getSymbol();
        index = atom->getIdx();
        if (updateCharge) charge = atom->getFormalCharge();
    }


    void Atom::replaceWithAtom(Atom newAtom) {
        bonds.updateAtom(newAtom);
        pairs.updateAtom(newAtom);
        exclusions.updateAtom(newAtom);
        angles.updateAtom(newAtom);
        dihedrals.updateAtom(newAtom);
        impropers.updateAtom(newAtom);
        
        std::cout << index << " atom index" << owningMol->name << std::endl;
        std::cout << owningMol->rdMol.getNumAtoms() << " heavy atoms in rdmol" << std::endl;

        // replace bonds in rdMol
        RDKit::Atom *atom = owningMol->rdMol.getAtomWithIdx(index);
        for (const auto &nbri : boost::make_iterator_range(owningMol->rdMol.getAtomBonds(atom))) {
            RDKit::Bond *bond_ = (owningMol->rdMol)[nbri];
            if (bond_->getBeginAtomIdx() == index) {
                bond_->setBeginAtomIdx(newAtom.index);
            } else {
                bond_->setEndAtomIdx(newAtom.index);
            }
        }
        owningMol->delAtom(*this);
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

    // void MonomerUnit::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
    //     for (auto &atom : atoms) { atom->removeParamsWithinAtomSet(atomSet); }
    // }

    int Monomer::delAtom(Atom &atom) {
        auto found = std::find(atoms.begin(), atoms.end(), &atom);
        rdMol.removeAtom(atom.index);
        if (found != atoms.end()) {
            // int index = std::distance(atoms.begin(), found);
            atoms.erase(found);
        }
        reindexAtoms(false);
        return atoms.size();
    };

    void Monomer::setRDMol(RDKit::RWMol mol) {
        RDKit::MolOps::addHs(mol);
        rdMol = mol;
        numAtoms = mol.getNumAtoms();
        atoms.reserve(numAtoms);

        for (int i = 0; i < numAtoms; i++) {
            Atom *newAtom = new Atom(*this, i);
            atoms.emplace_back(newAtom);
        }

        for (int i = 0; i < mol.getNumBonds(false); i++) {
            auto rdBond = mol.getBondWithIdx(i);
            
            std::array<unsigned int, 2> indices = {rdBond->getBeginAtomIdx(),
                                          rdBond->getEndAtomIdx()};
            Bond *newBond = new Bond(*this, indices);
            for (auto &atom : newBond->atoms) {
                atom->bonds.addParam(newBond);
            }

        }

    }

    void Monomer::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
        for (auto &atom : atoms) { atom->removeParamsWithinAtomSet(atomSet); }
        unsigned int numBonds = rdMol.getNumBonds();
        std::vector< std::pair<int, int> > toDelete;
        for (int i = 0; i < numBonds; i++) {
            RDKit::Bond* rdBond = rdMol.getBondWithIdx(i);
            unsigned int begin = rdBond->getBeginAtomIdx();
            unsigned int end = rdBond->getEndAtomIdx();
            if ((atomSet.count(atoms[begin]) == 1) &&
                (atomSet.count(atoms[end]) == 1)) {
                    toDelete.insert(toDelete.begin(), std::make_pair(begin, end));
                }
        }
        for (auto &el : toDelete) {
            rdMol.removeBond(el.first, el.second);
        };
    }

    void Monomer::reindexAtoms(bool setOwningMol) {
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i]->index = i;
            if (setOwningMol) {atoms[i]->owningMol = this;};
        }
    }

    
    MonomerUnit::MonomerUnit(Monomer mol) {
        // atoms
        rdMol = RDKit::RWMol(mol.rdMol);
        atoms.reserve(mol.atoms.size());
        for (auto &atom : mol.atoms) {
            Atom *newAtom = new Atom(*atom);
            atoms.push_back(newAtom);
            newAtom->owningMol = this;
        }

        // params
        copyParams<Bond, BondVector>(mol);
        copyParams<Pair, PairVector>(mol);
        copyParams<Exclusion, ExclusionVector>(mol);
        copyParams<Angle, AngleVector>(mol);
        copyParams<Dihedral, DihedralVector>(mol);
        copyParams<Improper, ImproperVector>(mol);
    };



};

#endif