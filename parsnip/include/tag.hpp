#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "atom.hpp"
namespace polytop {
    class Tag {
        public:
            std::string name;
            std::vector<Atom*> atoms;

            Tag(std::string name, std::vector<Atom*> atoms);

            std::vector<unsigned int> getAtomIndices();
            std::size_t size();
            std::size_t minSize(Tag* other);

            // void setAtomIndices(std::vector<int> indices);


        // private:
        //     std::vector<int> atomIndices;

    };
};