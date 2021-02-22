#include <unordered_map>
#include <array>
#include <set>
#include <vector>
#include <string>
#include <functional>

#ifndef _PT_ATOM_H
#define _PT_ATOM_H

#include <GraphMol/GraphMol.h>

#include "paramvector.h"

namespace polytop {
    class Monomer;
    class MonomerUnit;
    
    
    class Atom {

        public:

            // Atom();
            // Atom(const Atom &rhs);
            // Atom(Monomer &mol);
            // Atom(MonomerUnit &mol);
            Atom(Monomer &mol, int atomIndex);
            Atom(MonomerUnit &mol, int atomIndex);

            Monomer *owningMol;
            int index = 0;
            double charge = 0;
            std::string name = "X";

            BondVector bonds = BondVector(this);
            PairVector pairs = PairVector(this);
            ExclusionVector exclusions = ExclusionVector(this);
            AngleVector angles = AngleVector(this);
            DihedralVector dihedrals = DihedralVector(this);
            ImproperVector impropers = ImproperVector(this);

            void setRDAtom(RDKit::Atom *atom, bool updateCharge=false);

            void replaceWithAtom(Atom newAtom);
            void removeParamsWithinAtomSet(std::set<Atom*> atomSet);

    };
};

#endif