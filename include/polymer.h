
#include <unordered_map>

#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#ifndef _PT_POLYMER_H
#define _PT_POLYMER_H

#include "monomer.h"

namespace polytop {

    using MonoPolyAtomIndexVect = RDKit::MatchVectType;

    class Polymer : public MonomerUnit {
        // friend class Monomer;
        // friend class MonomerUnit;

        public:
            Polymer(std::string name="UNK");
            // std::string name;
            void addMonomerUnit(MonomerUnit unit, MonoPolyAtomIndexVect atomIndices,
                               bool replacePolymerAtoms=true);

            void addMonomerUnit(MonomerUnit unit);
            // void removeParamsWithinAtomSet(std::set<Atom*> atomSet);
            
            // int addMonomer(Monomer monomer, bool useMonomerParams=true,
            //                bool replacePolymerAtoms=true);

            
            double charge = 0;

            // void reindexAtoms(bool setOwningMol=true);




            // RDKit::RWMol rdMol = RDKit::RWMol();
            // std::vector<MonomerUnit> monomers;
            // std::vector<Atom*> atoms;
            // std::unordered_multimap<std::array<int, 2>, Bond> bonds;

    };
};

#endif