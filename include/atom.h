#include <unordered_map>
#include <array>
#include <set>
#include <vector>
#include <string>
#include <functional>

#ifndef _PT_ATOM_H
#define _PT_ATOM_H

// #include <GraphMol/GraphMol.h>

// #include "bond.h"

#include "paramvector.h"

namespace polytop {

    typedef std::hash<std::string> hash_string;

    class Monomer;
    class MonomerUnit;
    class Atom {

        public:

            Atom();
            Atom(Monomer &mol);
            Atom(MonomerUnit &mol);
            // RDKit::Atom rdAtom;

            int index = 0;
            void replaceWithAtom(Atom newAtom);
            std::unique_ptr<Atom> copy();
            void removeParamsWithinAtomSet(std::set<Atom*> atomSet);
            void copyToMol(MonomerUnit mol);

            Monomer *owningMol;
            std::string name = "X";
            BondVector bonds = BondVector(this);
            PairVector pairs = PairVector(this);
            ExclusionVector exclusions = ExclusionVector(this);
            AngleVector angles = AngleVector(this);
            DihedralVector dihedrals = DihedralVector(this);
            ImproperVector impropers = ImproperVector(this);
            

    };
};

#endif