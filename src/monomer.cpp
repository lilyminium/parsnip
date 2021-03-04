#include <algorithm>
#include <set>
#include <string>
#include <iostream>
#include <utility>

#ifndef _PT_MONOMER_C
#define _PT_MONOMER_C

#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>

#include "monomer.hpp"

namespace polytop {
    

    int Monomer::addAtom(Atom &atom) {
        atoms.emplace_back(&atom);
        return atoms.size();
    }

    Monomer::Monomer(const RDKit::ROMol &mol, std::string name) : name(name) {
        setRDMol(RDKit::RWMol(mol));
    }

    std::vector<unsigned int> Monomer::getAtomIndices() {
        std::vector<unsigned int> vec;
        for (auto &at : atoms) {
            vec.emplace_back(at->index);
        }
        return vec;
    }
    
    std::size_t Monomer::delAtomFromVector(Atom *atom) {
        auto found = std::find(atoms.begin(), atoms.end(), atom);
        if (found != atoms.end()) {
            atoms.erase(found);
        }
        return atoms.size();
    }

    std::size_t Monomer::delAtomsFromVector(std::vector<Atom*> toDelete) {
        for (auto &at : toDelete) {
            delAtomFromVector(at);
        }
        return atoms.size();
    }

    std::size_t Monomer::delAtomsFromVector(std::set<Atom*> toDelete) {
        for (auto &at : toDelete) {
            delAtomFromVector(at);
        }
        return atoms.size();
    }

    std::size_t Monomer::removeAtom(Atom *atom) {
        rdMol.removeAtom(atom->rdAtom->getIdx());
        delAtomFromVector(atom);
        reindexAtoms();
        cleanParams();
        return atoms.size();
    };

    std::size_t Monomer::cleanTags() {
        for (auto it = tags.begin(); it != tags.end();) {
            if (containsAllAtoms(it->second->atoms)) {
                ++it;
            } else {
                it = tags.erase(it);
            }
        }
        return tags.size();
    }

    void Monomer::setRDMol(RDKit::RWMol mol) {
        RDKit::MolOps::addHs(mol, true, true);
        rdMol = mol;
        auto numAtoms = rdMol.getNumAtoms();
        atoms.reserve(numAtoms);

        for (int i = 0; i < numAtoms; i++) {
            Atom *newAtom = new Atom(rdMol.getAtomWithIdx(i));
            atoms.emplace_back(newAtom);
        }

        for (RDKit::ROMol::BondIterator bndIt = mol.beginBonds(); bndIt != mol.endBonds(); ++bndIt) {
            std::array<unsigned int, 2> indices = {(*bndIt)->getBeginAtomIdx(),
                                                   (*bndIt)->getEndAtomIdx()};
            std::array<Atom*, 2> atoms = getAtomsByIndices(indices);
            Bond *newBond = new Bond(atoms);
            bonds.addParam(newBond);
        }
    };

    std::vector<Atom*> Monomer::getAtomsByIndices(std::vector<unsigned int> indices) {
        std::vector<Atom*> vec;
        for (auto &i : indices) {
            vec.emplace_back(atoms[i]);
        }
        return vec;
    }

    void Monomer::replaceAtom(Atom *oldAtom, Atom *newAtom) {
        bonds.updateAtom(oldAtom, newAtom);
        pairs.updateAtom(oldAtom, newAtom);
        exclusions.updateAtom(oldAtom, newAtom);
        angles.updateAtom(oldAtom, newAtom);
        dihedrals.updateAtom(oldAtom, newAtom);
        impropers.updateAtom(oldAtom, newAtom);
        
        // replace bonds in rdMol
        RDKit::Atom *atom = oldAtom->rdAtom;
        unsigned int oldIndex = oldAtom->rdAtom->getIdx();
        unsigned int newIndex = newAtom->rdAtom->getIdx();
        for (const auto &nbri : boost::make_iterator_range(rdMol.getAtomBonds(atom))) {
            auto bond_ = rdMol[nbri];
            if (bond_->getBeginAtomIdx() == oldIndex) {
                bond_->setBeginAtomIdx(newAtom->rdAtom->getIdx());
                rdMol.addBond(newAtom->rdAtom, bond_->getEndAtom());
            } else {
                bond_->setEndAtomIdx(newAtom->rdAtom->getIdx());
                rdMol.addBond(bond_->getBeginAtom(), newAtom->rdAtom);
            }

        }
        rdMol.updatePropertyCache();
        removeAtom(oldAtom);
    };

    void Monomer::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
        removeParamsWithinAtomSetFromVector(atomSet);
        unsigned int numBonds = rdMol.getNumBonds();
        std::vector< std::pair<unsigned int, unsigned int> > toDelete;
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

    std::size_t Monomer::delAtomByIndex(unsigned int index) {
        auto at = atoms[index];
        return delAtomFromVector(at);
    }

    std::size_t Monomer::removeAtomByIndex(unsigned int index) {
        auto at = atoms[index];
        return removeAtom(at);
    }

    void Monomer::copyParamsFrom(Monomer* other) {
        copyParamsFrom(other->bonds);
        copyParamsFrom(other->angles);
        copyParamsFrom(other->dihedrals);
        copyParamsFrom(other->impropers);
        copyParamsFrom(other->pairs);
        copyParamsFrom(other->exclusions);
    }

    void Monomer::addParamsFromMolecule(Monomer* other) {
        addParamsFromVector(other->bonds);
        addParamsFromVector(other->angles);
        addParamsFromVector(other->dihedrals);
        addParamsFromVector(other->impropers);
        addParamsFromVector(other->pairs);
        addParamsFromVector(other->exclusions);
    }

    void Monomer::reindexAtoms() {
        for (unsigned int i = 0; i < atoms.size(); i++) {
            atoms[i]->index = i;
            // atoms[i]->rdAtom->setIdx(i);
        }
    }

    Tag* Monomer::getFirstTagByName(std::string name) {
        auto it = tags.find(name);
        return it->second;
    }

    Tag* Monomer::getLastTagByName(std::string name) {
        // TODO i'm sure I can do this better
        auto range = tags.equal_range(name);
        auto lastit = range.first;
        for (auto it = range.first; it != range.second; ++it) {
            lastit = it;
        }
        return lastit->second;
    }

    void Monomer::cleanParams() {
        bonds.removeParamsNotInAtoms(atoms);
        angles.removeParamsNotInAtoms(atoms);
        dihedrals.removeParamsNotInAtoms(atoms);
        impropers.removeParamsNotInAtoms(atoms);
        pairs.removeParamsNotInAtoms(atoms);
        exclusions.removeParamsNotInAtoms(atoms);
    }

    void Monomer::removeParamsWithinAtomSetFromVector(std::set<Atom*> atomSet) {
        bonds.removeParamsWithinAtomSet(atomSet);
        pairs.removeParamsWithinAtomSet(atomSet);
        exclusions.removeParamsWithinAtomSet(atomSet);
        angles.removeParamsWithinAtomSet(atomSet);
        dihedrals.removeParamsWithinAtomSet(atomSet);
        impropers.removeParamsWithinAtomSet(atomSet);
    }

    std::size_t Monomer::addTag(std::string name, std::vector<Atom*> atomVector) {
        Tag *tag = new Tag(name, atomVector);
        return addTag(tag);
    }

    std::size_t Monomer::addTag(std::string name, std::vector<unsigned int> indices) {
        auto atomVector = getAtomsByIndices(indices);
        return addTag(name, atomVector);
    }

    std::size_t Monomer::addTag(Tag* tag) {
        tags.emplace(tag->name, tag);
        return tags.size();
    }

    std::size_t Monomer::copyTagFrom(Tag* tag) {
        std::vector<unsigned int> indices = tag->getAtomIndices();
        return addTag(tag->name, indices);
    }

    std::size_t Monomer::safeAddTag(Tag* tag) {
        if (containsAllAtoms(tag->atoms)) {
            tags.emplace(tag->name, tag);
        }
        return tags.size();
    }

    bool Monomer::containsAtom(Atom* atom) {
        if (std::find(atoms.begin(), atoms.end(), atom) == atoms.end()) {
            return false;
        }
        return true;
    }

    bool Monomer::containsAllAtoms(std::vector<Atom*> atomVector) {
        for (auto &el : atomVector) {
            if (std::find(atoms.begin(), atoms.end(), el) == atoms.end()) {
                return false;
            }
        }
        return true;
    }

    std::size_t Monomer::unsafeAddTagsFrom(Monomer* other) {
        for (auto &el : other->tags) {
            addTag(el.second);
        }
        return tags.size();
    }

    std::size_t Monomer::copyTagsFrom(Monomer* other) {
        for (auto &el : other->tags) {
            copyTagFrom(el.second);
        }
        return tags.size();
    }

    std::size_t Monomer::addTagsFrom(Monomer* other, bool checkContainsAtoms) {
        if (checkContainsAtoms) return safeAddTagsFrom(other);
        return unsafeAddTagsFrom(other);
    }

    std::size_t Monomer::safeAddTagsFrom(Monomer* other) {
        for (auto &el : other->tags) {
            safeAddTag(el.second);
        }
        return tags.size();
    }

    




    MonomerUnit::MonomerUnit(Monomer *mol) {
        // atoms
        rdMol = RDKit::RWMol(mol->rdMol);
        atoms.reserve(mol->atoms.size());
        for (auto &atom : mol->atoms) {
            Atom *newAtom = new Atom(*atom);
            atoms.push_back(newAtom);
        }

        copyParamsFrom(mol);
        copyTagsFrom(mol);
    };



};

#endif