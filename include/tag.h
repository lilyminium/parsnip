#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace polytop {

    class Monomer;
    class Atom;

    class Tag {
        public:
            std::string name;
            std::vector<Atom*> atoms;
            Monomer *owningMol;

            Tag(Monomer *mol, std::string tagName, std::vector<int> indices);

            void setAtomIndices(std::vector<int> indices);

        private:
            std::vector<int> atomIndices;

    }
}