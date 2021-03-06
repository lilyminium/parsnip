
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include <GraphMol/MonomerInfo.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "utils.hpp"

namespace polytop {

std::vector<unsigned int> getNeighborInts(RDKit::RWMol rdMol,
                                              std::vector<unsigned int> indices,
                                              std::size_t numNeighbors) {
        std::vector<unsigned int> neighbors;
        std::set<unsigned int> seenIndices;
        for (auto ix : indices) {
            seenIndices.insert(ix);
        }

        for (auto ix : indices) {
            std::size_t currentNumNeighbors = numNeighbors;
            std::set<unsigned int> latestNeighborIndices = {ix};

            while (currentNumNeighbors > 0) {
                std::set<unsigned int> currentNeighborIndices = {};
                for (auto lix : latestNeighborIndices) {
                    RDKit::Atom *rdAtom = rdMol.getAtomWithIdx(lix);
                    auto rdNeighbors = rdMol.getAtomNeighbors(rdAtom);
                    auto range = boost::make_iterator_range(rdNeighbors);

                    for (const auto &nbri : range) {
                        const auto &nbr = rdMol[nbri];
                        auto nix = nbr->getIdx();
                        if ((currentNeighborIndices.find(nix) == currentNeighborIndices.end()) &&
                            (latestNeighborIndices.find(nix) == latestNeighborIndices.end()) &&
                            (seenIndices.find(nix) == seenIndices.end()) &&
                            (std::find(neighbors.begin(), neighbors.end(), nix) == neighbors.end())) {
                                currentNeighborIndices.insert(nbr->getIdx());
                            }
                        
                    }
                }
                for (auto el : currentNeighborIndices) {
                    neighbors.emplace_back(el);
                }
                latestNeighborIndices = currentNeighborIndices;
                currentNumNeighbors -= 1;
            }
        }
    return neighbors;
    }

    RDKit::RWMol* subsetRDMol(RDKit::RWMol rdMol, std::vector<unsigned int> indices,
                              std::vector<unsigned int> neighbors) {
        RDKit::RWMol *newMol = new RDKit::RWMol();
        if (!indices.size()) return newMol;

        RDKit::Conformer *newConf = new RDKit::Conformer(newMol->getNumAtoms());
        newMol->addConformer(newConf);
        auto rdConf = rdMol.getConformer();

        // add atoms and geometry
        for (std::size_t i = 0; i < indices.size(); i++) {
            auto rdAtom = rdMol.getAtomWithIdx(indices[i]);
            newMol->addAtom(rdAtom);
            newConf->setAtomPos(i, rdConf.getAtomPos(indices[i]));
        }

        auto length = indices.size();

        for (std::size_t i = 0; i < neighbors.size(); i++) {
            auto rdAtom = rdMol.getAtomWithIdx(neighbors[i]);
            newMol->addAtom(rdAtom);
            auto newAtom = newMol->getAtomWithIdx(length);
            auto info = getMonomerInfo(newAtom);
            info->setChainId("+");
            newAtom->setMonomerInfo(info);
            newConf->setAtomPos(length, rdConf.getAtomPos(neighbors[i]));
            length++;
        }

        std::vector<unsigned int> combined;
        combined.insert(combined.end(), indices.begin(), indices.end());
        combined.insert(combined.end(), neighbors.begin(), neighbors.end());
        
        // add bonds
        for (std::size_t i = 0; i < combined.size(); i++) {
            for (std::size_t j = i; j < combined.size(); j++) {
                auto newBond = rdMol.getBondBetweenAtoms(combined[i], combined[j]);
                if (newBond) {
                    newMol->addBond(i, j);
                }
            }
        }

        return newMol;
    };

    // std::vector<RDKit::RWMol*> getUniqueRDMols(std::vector<RDKit::RWMol*> rdMols) {
    //     std::unordered_map<std::string, RDKit::RWMol*> map;
    //     std::vector<std::string> smiles;
    //     // canonical smiles comparison is probably fine
    //     for (auto &mol : rdMols) {
    //         auto smi = RDKit::MolToSmiles(*mol);
    //         map[smi] = mol;
    //         smiles.emplace_back(smi);

    //         // for (auto it = mol->beginAtoms(); it != mol->endAtoms()) {
    //         //     // if (*it)->get
    //         // }
    //     }
    //     std::vector<RDKit::RWMol*> unique;
    //     for (auto el : map) {
    //         unique.emplace_back(el.second);
    //     }
    //     return unique;
    // };

    // RDKit::AtomPDBResidueInfo* getMonomerInfo(RDKit::Atom atom) {
    //     // This is annoying.
    //     auto mInfo = atom->getMonomerInfo();
    //     RDKit::AtomPDBResidueInfo* info = reinterpret_cast<RDKit::AtomPDBResidueInfo*>(mInfo);
    //     return info;
    // }




    std::vector<RDKit::RWMol*> getUniqueRDMols(std::vector<RDKit::RWMol*> rdMols) {
        std::unordered_map<std::string, RDKit::RWMol*> map;
        std::vector<std::string> smiles;
        // canonical smiles comparison is probably fine
        for (auto &mol : rdMols) {
            auto smi = RDKit::MolToSmiles(*mol);
            map[smi] = mol;
            smiles.emplace_back(smi);
        }

        std::vector<RDKit::RWMol*> unique;
        std::unordered_map<std::string, unsigned int> indices;
        for (auto el : map) {
            unique.emplace_back(el.second);
            indices[el.first] = indices.size();
        }

        for (auto &mol : unique) {
            for (auto at = mol->beginAtoms(); at != mol->endAtoms()) {
                auto info = getMonomerInfo(*at);
                if (info->getChainId() == "+") { // cap atoms
                    int index = info->getResidueNumber() - 1;
                    std::string smi = smiles[index];
                    int newIndex = indices[smi];
                    info->setResidueNumber(newIndex);
                }
            }
        }

        return unique;
    };

    std::vector<unsigned int> getUnitIndices(std::vector<RDKit::RWMol*> uniqueUnits,
                                             std::vector<RDKit::RWMol*> allUnits) {
        std::unordered_map<std::string, unsigned int> unitMap;
        for (std::size_t i = 0; i < uniqueUnits.size(); ++i) {
            auto smi = RDKit::MolToSmiles(uniqueUnits[i]);
            unitMap[smi] = i;
        }

        std::vector<unsigned int> unitIndices;
        unitIndices.reserve(allUnits.size());
        for (auto mol : allUnits) {
            unitIndices.push_back(unitMap[RDKit::MolToSmiles(mol)]);
        }
        return unitIndices;
    }
    


    std::vector<unsigned int> splitIntegerIntoRatio(unsigned int total,
                                                    std::vector<double> ratio) {
        double ratioTotal = std::accumulate(ratio.begin(), ratio.end(),
                                            decltype(ratio)::value_type(0));

        std::vector<double> currentFraction;
        std::vector<double> fraction;
        fraction.reserve(ratio.size());
        currentFraction.reserve(ratio.size());

        for (auto r : ratio) {
            fraction.push_back(r/ratioTotal);
            currentFraction.push_back(0);
        }

        
        std::vector<unsigned int> currentRatio;
        double smallest = *std::min_element(ratio.begin(), ratio.end());
        for (auto &el : ratio) {
            currentRatio.push_back(el/smallest);
        }
        unsigned int currentTotal = std::accumulate(currentRatio.begin(), currentRatio.end(),
                                                    decltype(currentRatio)::value_type(0));
        
        int quotient = (total / currentTotal) - 1;
        if (quotient) {
            for (std::size_t i = 0; i < currentRatio.size(); i++) {
                currentRatio[i] *= quotient;
            };
        };

        currentTotal = std::accumulate(currentRatio.begin(), currentRatio.end(),
                                       decltype(currentRatio)::value_type(0));
        
        while (currentTotal < total) {
            for (std::size_t i = 0; i < currentRatio.size(); i++ ) {
                currentFraction[i] = double(currentRatio[i]) / total;
                currentFraction[i] /= fraction[i];
            };

            auto it = std::min_element(currentFraction.begin(), currentFraction.end());
            unsigned int minIndex = std::distance(currentFraction.begin(), it);
            currentRatio[minIndex] += 1;
            currentTotal = std::accumulate(currentRatio.begin(), currentRatio.end(),
                                           decltype(currentRatio)::value_type(0));
        }

        return currentRatio;
    };

    RDKit::RWMol* copyRDMol(RDKit::RWMol rdMol) {
        RDKit::RWMol* newMol = new RDKit::RWMol(rdMol);
        for (std::size_t i = 0; i < newMol->getNumAtoms(); i++) {
            auto info = rdMol.getAtomWithIdx(i)->getMonomerInfo();
            // auto newInfo = new RDKit::AtomPDBResidueInfo(*info);
            newMol->getAtomWithIdx(i)->setMonomerInfo(info->copy());
        }
        return newMol;
    };

    // std::unordered_multimap<std::string, RDKit::RWMol*> mapRDMols(std::vector<RDKit::RWMol*> rdMols) {
    //     std::unordered_multimap<std::string, RDKit::RWMol*> map;
    //     for (auto &mol : rdMols) {
    //         auto smi = RDKit::MolToSmiles(*mol);
    //         map[smi] = mol;
    //     }
    //     return map;
    // }

    // std::vector<RDKit::RWMol*> getUniqueMols(std::unordered_multimap<std::string, RDKit::RWMol*> map,
    //                                          std::vector<RDKit::RWMol*> rdMols) {
        
    // }




};