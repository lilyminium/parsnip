#include <iostream>
#include <string>


#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>


// #include "atom.h"

#include "atom.hpp"
#include "monomer.hpp"
#include "polymer.hpp"
// #include "param.h"

int main() {

    using namespace polytop;

    // read file into Monomer
    RDKit::ROMol *mol1 = RDKit::PDBFileToMol("/Users/lily/workspace/synthezymes/01_structures/pdb/m1_u4.pdb", false, false);
    std::cout << mol1->getNumConformers() << std::endl;
    Monomer* mol = new Monomer(*mol1);
    std::cout << mol->atoms.size() << " monomer atoms " << std::endl;
    std::cout << mol->rdMol.getNumAtoms() << " rd monomer atoms " << std::endl;
    std::cout << mol->bonds->data.size() << " monomer bonds " << std::endl;

    // add Tags
    mol->addTag("left", {43, 46, 44, 45});
    mol->addTag("right", {48, 51, 53, 52});
    std::cout << mol->tags.size() << " monomer tags " << std::endl;

    MonomerUnit *unit = new MonomerUnit(mol);
    std::cout << unit->atoms.size() << " unit atoms " << std::endl;
    std::cout << unit->tags.size() << " unit tags \n" << std::endl;

    Polymer pol = Polymer("PolyAla");
    std::cout << pol.name << std::endl;
    pol.addMonomerUnit(unit);
    std::cout << pol.atoms.size() << " polymer atoms " << std::endl;
    std::cout << pol.rdMol.getNumAtoms() << " rd polymer atoms " << std::endl;
    std::cout << pol.bonds->data.size() << " polymer bonds\n" << std::endl;
    std::cout << pol.atoms[0]->name << std::endl;

    // MonoPolyAtomIndexVect vec = { {43, 48}, {46, 51}, {44, 53}, {45, 52} };

    // MonomerUnit *unit2 = new MonomerUnit(mol);
    // pol.addMonomerUnit(unit2, vec, true);
    // std::cout << unit2->bonds->data.size() << " unit bonds " << std::endl;
    // std::cout << pol.atoms.size() << " polymer atoms " << std::endl;

    pol.addMonomer(mol, "left", "right", true);
    std::cout << pol.bonds->data.size() << " polymer bonds " << std::endl;
    std::cout << pol.tags.size() << " unit tags \n" << std::endl;

    pol.addMonomer(mol, "left", "right", true);
    std::cout << pol.tags.size() << " unit tags" << std::endl;
    std::cout << pol.atoms.size() << " polymer atoms " << std::endl;
    std::cout << pol.bonds->data.size() << " polymer bonds " << std::endl;

    RDKit::MolToPDBFile(pol.rdMol, "/Users/lily/Desktop/test.pdb");

    RDKit::RWMol *newMol = pol.getRDUnit(pol.units[1], 3);
    std::cout << newMol->getNumAtoms() << std::endl;
    std::cout << "new mol" << std::endl;
    RDKit::MolToPDBFile(*newMol, "/Users/lily/Desktop/capped.pdb");

    


    // RDKit::ROMol *mol1 = RDKit::SmilesToMol("Cc1ccccc1");
    
    // Monomer mol = Monomer(*mol1);
    // Polymer pol = Polymer("PolyAla");
    // // Atom atom = Atom(mol);
    // // mol->addAtom(atom);
    // std::cout << mol->resName << std::endl;
    // std::cout << mol->atoms.size() << std::endl;
    // std::cout << mol->atoms[0]->name << std::endl;
    // std::cout << mol->atoms[0]->bonds.data.size() << std::endl;
    // // mol->delAtom(atom);
    // // std::cout << mol->atoms.size() << std::endl;

    std::cout << "hmm" << std::endl;
}