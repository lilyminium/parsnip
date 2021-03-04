
#include <vector>
#include <utility>

#include "tag.hpp"

namespace polytop {
    Tag::Tag(std::string name, std::vector<Atom*> atoms) : name(name), atoms(atoms) {};

    std::vector<unsigned int> Tag::getAtomIndices() {
        std::vector<unsigned int> indices;
        for (auto &atom : atoms) {
            indices.push_back(atom->index);
        }
        return indices;
    };

    std::size_t Tag::size() {
        return atoms.size();
    }

    std::size_t Tag::minSize(Tag* other) {
        std::size_t value = atoms.size();
        if (other->size() < value) return other->size();
        return value;
    }
}