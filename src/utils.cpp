
#include <iostream>
#include <vector>
#include <algorithm>

#include <GraphMol/MonomerInfo.h>

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

    RDKit::RWMol* subsetRDMol(RDKit::RWMol rdMol, std::vector<unsigned int> indices) {
        RDKit::RWMol *newMol = new RDKit::RWMol();
        if (!indices.size()) return newMol;
        auto length = indices.size();

        RDKit::Conformer *newConf = new RDKit::Conformer(newMol->getNumAtoms());
        newMol->addConformer(newConf);
        auto rdConf = rdMol.getConformer();

        // add atoms and geometry
        for (std::size_t i = 0; i < length; i++) {
            auto rdAtom = rdMol.getAtomWithIdx(indices[i]);
            newMol->addAtom(rdAtom);
            newConf->setAtomPos(i, rdConf.getAtomPos(indices[i]));
        }
        
        // add bonds
        for (std::size_t i = 0; i < length; i++) {
            for (std::size_t j = i; j < length; j++) {
                auto newBond = rdMol.getBondBetweenAtoms(indices[i], indices[j]);
                if (newBond) {
                    newMol->addBond(i, j);
                }
            }
        }

        return newMol;
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




};