
#include <unordered_map>
#include <set>
#include <vector>
#include <string>

// #include <GraphMol/GraphMol.h>

// #include "atom.h"
// #include "bond.h"

#ifndef _PT_MONOMER_H
#define _PT_MONOMER_H

namespace polytop {
    class Atom;

    class Monomer {
        // friend class Polymer;

        public:
            std::string resName;
            Monomer(std::string name="UNK");

            int addAtom(Atom &atom);
            // int addAtom(RDKit::Atom atom);

            int delAtom(Atom &atom);

            // RDKit::ROMol rdMol;
            std::vector<Atom*> atoms;
            
            
    };


    class MonomerUnit : Monomer {

        public:
            MonomerUnit(Monomer mol);

            void removeParamsWithinAtomSet(std::set<Atom*> atomSet);


            // RDKit::RWMol rdMol;
            std::vector<Atom*> atoms;
    };
};

#endif