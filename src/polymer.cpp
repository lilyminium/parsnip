#include <vector>

#include <GraphMol/FMCS/FMCS.h>
#include <GraphMol/GraphMol.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
// #include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolAlign/AlignMolecules.h>


#include "polymer.h"
#include "monomer.h"
#include "atom.h"

namespace polytop {

    Polymer::Polymer(std::string name) {
        name = name;
    };

    // int Polymer::addMonomerUnit(MonomerUnit unit, MonoPolyAtomIndexVect atomIndices,
    //                             bool replacePolymerAtoms) {
    //     // align positions
    //     RDKit::MolAlign::alignMol(unit.rdMol, rdMol, 0, 0, &atomIndices);
    //     // move monomer atoms into polymer
    //     addMonomerUnit(unit);

    //     // remove doubled-up atoms and topology parameters
    //     if (replacePolymerAtoms) {
    //         std::set<int> polyAtomIndices;
    //         for (auto& el : atomIndices) {
    //             polyAtomIndices.insert(el.second);
    //         }
    //         removeParamsWithinAtomSet(polyAtomIndices);

    //         for (auto& el : atomIndices) {
    //             int polyAtomIndex = el.second;
    //             int monoAtomIndex = el.first;
                
    //             atoms[polyAtomIndex].replaceWithAtom(unit.atoms[monoAtomIndex]);

    //         };
    //     };
    // };

    // int Polymer::addMonomerUnit(MonomerUnit unit) {
    //     // add atoms
    //     atoms.insert(atoms.end(), unit.atoms.begin(), unit.atoms.end());
    //     reindexAtoms(true);

    //     // add parameters and stuff
    //     bonds.insert(unit.bonds.begin(), unit.bonds.end());
    // };

    // void Polymer::reindexAtoms(bool setOwningMol) {
    //     for (size_t i = 0; i < atoms.size(); i++) {
    //         atoms[i].index = i;
    //         if (setOwningMol) {atoms[i].setOwningMol(&this); };
    //     }
    // }

    // int Polymer::addMonomer(Monomer monomer, bool useMonomerParams,
    //                         bool replacePolymerAtoms) {
    //     std::vector<RDKit::ROMOL_SPTR> mols;
    //     mols.emplace_back(monomer.rdMol);
    //     mols.emplace_back(rdMol);
    //     RDKit::MCSResult mcsre = findMCS(mols);

    //     if (mcsre.isCompleted()) {
    //         std::shared_ptr<RDKit::ROMol> match(RDKit::SmartsToMol(mcsre.SmartsString));
    //         // typedef std::vector<std::pair<int, int>> MatchVectType;
    //         RDKit::MatchVectType polyMatch, monoMatch;
    //         RDKit::SubstructMatch(rdMol, *match, polyMatch);
    //         RDKit::SubstructMatch(monomer.rdMol, *match, monoMatch);

    //         auto monoMatch = RDKit::SubstructMatch(monomer.rdMol, *match);

            
    //     } else {
    //         // just add monomer
    //     }

        
    // }
};