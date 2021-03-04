#include <vector>
#include <utility>
#include <iostream>


#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/GraphMol.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
// #include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include <GraphMol/MolAlign/AlignMolecules.h>


#include "polymer.hpp"

namespace polytop {

    void Polymer::addMonomer(Monomer *mol,
                            std::string monomerTag,
                            std::string polymerTag,
                            bool replacePolymerAtoms) {
            MonomerUnit *unit = new MonomerUnit(mol);
            addMonomerUnit(unit, monomerTag, polymerTag,
                           replacePolymerAtoms);
        };

    void Polymer::addMonomer(Monomer *mol,
                            Tag* monomerTag,
                            Tag* polymerTag,
                            bool replacePolymerAtoms) {
            MonomerUnit *unit = new MonomerUnit(mol);
            addMonomerUnit(unit, monomerTag, polymerTag,
                           replacePolymerAtoms);
        };

    void Polymer::addMonomerUnit(MonomerUnit *unit,
                                std::string monomerTagName,
                                std::string polymerTagName,
                                bool replacePolymerAtoms) {
            Tag* monomerTag = unit->getLastTagByName(monomerTagName);
            Tag* polymerTag = getLastTagByName(polymerTagName);
            addMonomerUnit(unit, monomerTag, polymerTag, replacePolymerAtoms);
        }

    void Polymer::addMonomerUnit(MonomerUnit *unit,
                            Tag* monomerTag,
                            Tag* polymerTag,
                            bool replacePolymerAtoms) {
            
            if (monomerTag == nullptr || polymerTag == nullptr ) {
                addFromMoleculeUnit(unit);
                return;
            }
            MonoPolyAtomIndexVect *atomIndices = new MonoPolyAtomIndexVect();

            std::size_t nTagAtoms = monomerTag->minSize(polymerTag);
            for (std::size_t i = 0; i < nTagAtoms; i++) {
                auto pair = IndexPair(monomerTag->atoms[i]->index,
                                      polymerTag->atoms[i]->index);
                atomIndices->emplace_back(pair);
            }
            
            addMonomerUnit(unit, *atomIndices, replacePolymerAtoms);
        }


    void Polymer::addMonomerUnit(MonomerUnit *unit,
                                 MonoPolyAtomIndexVect atomIndices,
                                 bool replacePolymerAtoms) {
        // align positions
        RDKit::MolAlign::alignMol(unit->rdMol, rdMol, 0, 0, &atomIndices);
        
        // remove doubled-up atoms and topology parameters
        std::set<Atom*> toDelete;
        std::vector<std::pair<Atom*, Atom*>> atomVector;
        atomVector.reserve(atomIndices.size());

        if (replacePolymerAtoms) {
            for (auto &el : atomIndices) {
                auto polyAtom = atoms[el.second];
                auto monoAtom = unit->atoms[el.first];
                toDelete.insert(polyAtom);
                atomVector.emplace_back(AtomPair(monoAtom, polyAtom));
            }
            removeParamsWithinAtomSet(toDelete);
        }
        else {
            for (auto &el : atomIndices) {
                auto polyAtom = atoms[el.second];
                auto monoAtom = unit->atoms[el.first];
                toDelete.insert(monoAtom);
                atomVector.emplace_back(AtomPair(polyAtom, monoAtom));
            }
            unit->removeParamsWithinAtomSet(toDelete);
        }

        // move monomer atoms into polymer
        addMonomerUnitAtoms(unit);

        for (auto &el : atomVector) {
            replaceAtom(el.second, el.first, false);
        }

        delAtomsFromVector(toDelete);
        reindexAtoms();
        cleanUnitAtoms();
        addParamsFromMolecule(unit);
        cleanUnitParams();
        unsafeAddTagsFrom(unit);
        cleanTags();
    };

    void Polymer::addMonomerUnitAtoms(MonomerUnit *unit) {
        // add atoms
        std::size_t previousSize = rdMol.getNumAtoms();
        addUnit(unit);
        rdMol.insertMol(unit->rdMol);
        for (std::size_t i = 0; i < unit->atoms.size(); i++) {
            auto at = unit->atoms[i];
            at->rdAtom = rdMol.getAtomWithIdx(previousSize + i);
            // at->rdAtom->setMonomerInfo(at->resInfo);
            at->updateMonomerInfo();
        }
        atoms.insert(atoms.end(), unit->atoms.begin(), unit->atoms.end());
        
        reindexAtoms();
    };

    void Polymer::addUnit(MonomerUnit *unit) {
        units.emplace_back(unit);
        unit->setResNum(units.size());
    }

    void Polymer::addFromMoleculeUnit(MonomerUnit *unit) {
        // add atoms
        addMonomerUnitAtoms(unit);
        addParamsFromMolecule(unit);
        addTagsFrom(unit);
    };

    void Polymer::reorderAtomsByMonomerUnit() {
        atoms.clear();
        for (auto &mol : units) {
            atoms.insert(atoms.end(), mol->atoms.begin(), mol->atoms.end());
        }
        reindexAtoms();
    }

    // void Polymer::clearDeletedAtomsFromUnits() {
    //     for (auto &mol : units) {
    //         std::vector<Atom*> toDelete = {};
    //         for (auto &at : mol->atoms) {
    //             if (!containsAtom(at)) toDelete.emplace_back(at);
    //         }
    //         mol->delAtomsFromVector(toDelete);
    //     }

    // }

    void Polymer::cleanUnitAtoms() {
        for (auto &un : units) {
            for (std::size_t i = un->atoms.size(); i > 0; --i) {
                auto at = un->atoms[i];
                if (!containsAtom(at)) {
                    un->delAtomFromVector(at);
                }
            }
        }
    }

    RDKit::RWMol* Polymer::getRDUnit(MonomerUnit* unit, std::size_t numNeighbors) {
        std::vector<unsigned int> indices = unit->getAtomIndices();
        std::vector<unsigned int> neighbors = getNeighborInts(rdMol, indices, numNeighbors);
        
        indices.insert(indices.end(), neighbors.begin(), neighbors.end());
        return subsetRDMol(rdMol, indices);
    }

    void Polymer::cleanUnitParams() {
        for (auto un : units) {
            un->cleanParams();
        }
    }

    void Polymer::getCappedUnits(std::size_t numNeighbors) {
        reorderAtomsByMonomerUnit();
        cleanUnitAtoms();
        std::vector<RDKit::RWMol*> uniqueCappedMols;

        for (auto un : units) {
            auto capped = getRDUnit(un, numNeighbors);
        }


    }




    

    // int Polymer::addMonomer(Monomer monomer, bool useMonomerParams,
    //                         bool replacePolymerAtoms) {
    //     std::vector<RDKit::ROMOL_SPTR> mols;
    //     mols.emplace_back(monomer.rdMol);
    //     mols.emplace_back(rdMol);
    //     RDKit::MCSResult mcsre = findMCS(mols);

    //     if (mcsre.isCompleted()) {
    //         std::shared_ptr<RDKit::ROMol> match(RDKit::SmartsToMol(mcsre.SmartsString));
    //         // typedef std::vector<std::pair<int, int>> MatchVectType;
    //         RDKit::MatchVectType polyMatch, monoMatch;
    //         RDKit::SubstructMatch(rdMol, *match, polyMatch);
    //         RDKit::SubstructMatch(monomer.rdMol, *match, monoMatch);

    //         auto monoMatch = RDKit::SubstructMatch(monomer.rdMol, *match);

            
    //     } else {
    //         // just add monomer
    //     }

        
    // }
};