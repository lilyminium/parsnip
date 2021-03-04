
#include <unordered_map>

#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#ifndef _PT_POLYMER_H
#define _PT_POLYMER_H

#include "monomer.hpp"
#include "utils.hpp"

namespace polytop {

    using MonoPolyAtomIndexVect = RDKit::MatchVectType;
    typedef std::pair<unsigned int, unsigned int> IndexPair;
    typedef std::pair<Atom*, Atom*> AtomPair;

    class Polymer : public Monomer {

        public:
            Polymer(std::string name_="UNK") : Monomer(name_) {}

            std::vector<MonomerUnit*> units;

            void addMonomerUnit(MonomerUnit *unit, MonoPolyAtomIndexVect atomIndices,
                               bool replacePolymerAtoms=true);

            void addMonomerUnitAtoms(MonomerUnit *unit);
            void addFromMoleculeUnit(MonomerUnit *unit);
            void addMonomer(Monomer* mol, Tag* monomerTag=nullptr,
                            Tag* polymerTag=nullptr,
                            bool replacePolymerAtoms=true);
            void addMonomer(Monomer* mol, std::string monomerTag,
                            std::string polymerTag,
                            bool replacePolymerAtoms=true);
            void addMonomerUnit(MonomerUnit* unit, Tag* monomerTag=nullptr,
                                Tag* polymerTag=nullptr,
                                bool replacePolymerAtoms=true);

            void addMonomerUnit(MonomerUnit* unit, std::string monomerTag,
                                std::string polymerTag,
                                bool replacePolymerAtoms=true);
            void addUnit(MonomerUnit *unit);

            void reorderAtomsByMonomerUnit();
            void clearDeletedAtomsFromUnits();
            void getCappedUnits(std::size_t numNeighbors=3);
            void cleanUnitAtoms();
            void cleanUnitParams();

            RDKit::RWMol* getRDUnit(MonomerUnit* unit, std::size_t numNeighbors);

            // int addMonomer(Monomer monomer, bool useMonomerParams=true,
            //                bool replacePolymerAtoms=true);
    };
};

#endif