
#include <unordered_map>

#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#ifndef _PT_POLYMER_H
#define _PT_POLYMER_H

#include "monomer.h"

namespace polytop {

    using MonoPolyAtomIndexVect = RDKit::MatchVectType;

    class Polymer : MonomerUnit {
        // friend class Monomer;
        // friend class MonomerUnit;

        public:
            // int addMonomerUnit(MonomerUnit unit, MonoPolyAtomIndexVect atomIndices,
            //                    bool replacePolymerAtoms=true);

            // int addMonomerUnit(MonomerUnit unit);
            
            // int addMonomer(Monomer monomer, bool useMonomerParams=true,
            //                bool replacePolymerAtoms=true);

            // void reindexAtoms(bool setOwningMol=true);


        // protected:
            RDKit::RWMol rdMol;
            // std::vector<MonomerUnit> monomers;
            // std::vector<Atom> atoms;
            // std::unordered_multimap<std::array<int, 2>, Bond> bonds;

    };
};

#endif