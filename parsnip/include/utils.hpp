
#ifndef _PT_UTILS_H
#define _PT_UTILS_H

#include <GraphMol/GraphMol.h>


#include <vector>
#include <algorithm>
#include <numeric>

namespace polytop {

    std::vector<unsigned int> getNeighborInts(RDKit::RWMol rdMol,
                                              std::vector<unsigned int> indices,
                                              std::size_t numNeighbors=3);
            
    RDKit::RWMol* subsetRDMol(RDKit::RWMol rdMol, std::vector<unsigned int> indices,
                              std::vector<unsigned int> neighbors=std::vector<unsigned int>({}));

    RDKit::RWMol* copyRDMol(RDKit::RWMol rdMol);

    std::vector<RDKit::RWMol*> getUniqueRDMols(std::vector<RDKit::RWMol*>);


    std::vector<unsigned int> splitIntegerIntoRatio(unsigned int total,
                                                    std::vector<double> ratio);

};

#endif
