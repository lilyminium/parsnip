#include <unordered_map>
#include <array>
#include <set>
#include <vector>
#include <string>
#include <functional>

#ifndef _PT_ATOM_H
#define _PT_ATOM_H

#include <GraphMol/GraphMol.h>

namespace polytop {

    class Atom {

        public:
            Atom(RDKit::Atom *atom, std::string name="X", double charge=0);

            int index = 0;
            double charge = 0;
            std::string name = "X";
            RDKit::Atom *rdAtom;

            void setRDAtom(RDKit::Atom *atom);


    };
};

#endif