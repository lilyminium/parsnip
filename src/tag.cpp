
#include <vector>
#include <utility>

#include "tag.h"

namespace polytop {
    Tag::Tag(Monomer *mol, std::string tagName, std::vector<int> indices) {
        name = tagName;
        owningMol = mol;
        setAtomIndices(indices);
    }

    void Tag::setAtomIndices(std::vector<int> indices) {
        atomIndices = indices;
        atoms.clear();
        for (auto &el : indices) {
            atoms.emplace_back(owningMol->atoms[*el]);
        }
    }
}