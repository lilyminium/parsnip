#include <vector>
#include <utility>

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

    Polymer::Polymer(std::string name_) {
        name = name_;
    };

    void Polymer::addMonomerUnit(MonomerUnit unit, MonoPolyAtomIndexVect atomIndices,
                                bool replacePolymerAtoms) {
        // align positions
        RDKit::MolAlign::alignMol(unit.rdMol, rdMol, 0, 0, &atomIndices);
        

        // remove doubled-up atoms and topology parameters
        std::set<Atom*> atomSet;
        std::vector<std::pair<Atom*, Atom*>> atomVector;
        atomVector.reserve(atomIndices.size());


        if (replacePolymerAtoms) {
            for (auto &el : atomIndices) {
                atomSet.insert(atoms[el.second]);
                atomVector.emplace_back(std::pair<Atom*, Atom*>(unit.atoms[el.first],
                                                                atoms[el.second]));
            }
            for (auto &atom : atoms) { atom->removeParamsWithinAtomSet(atomSet); }
        }
        else {
            for (auto &el : atomIndices) {
                atomSet.insert(unit.atoms[el.first]);
                atomVector.emplace_back(std::pair<Atom*, Atom*>(atoms[el.second],
                                                                unit.atoms[el.first]));
            }
            for (auto &atom : unit.atoms) { atom->removeParamsWithinAtomSet(atomSet); }
        }

        // move monomer atoms into polymer
        addMonomerUnit(unit);

        for (auto &el : atomVector) {
            el.second->replaceWithAtom(*(el.first));
        }








        // if (replacePolymerAtoms) {
        //     std::set<int> polyAtomIndices;
        //     for (auto& el : atomIndices) {
        //         polyAtomIndices.insert(el.second);
        //     }
        //     removeParamsWithinAtomSet(polyAtomIndices);

        //     for (auto& el : atomIndices) {
        //         int polyAtomIndex = el.second;
        //         int monoAtomIndex = el.first;
                
        //         atoms[polyAtomIndex].replaceWithAtom(unit.atoms[monoAtomIndex]);

        //     };
        // };
    };

    void Polymer::addMonomerUnit(MonomerUnit unit) {
        // add atoms
        atoms.insert(atoms.end(), unit.atoms.begin(), unit.atoms.end());
        rdMol.insertMol(unit.rdMol);
        this->reindexAtoms(true);

        std::cout << rdMol.getNumAtoms() << " heavy atoms in rdmol reindex 0" << std::endl;
        std::cout << atoms[atoms.size()-1]->owningMol->rdMol.getNumAtoms() << " heavy atoms in rdmol reindex" << std::endl;
    };

    

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