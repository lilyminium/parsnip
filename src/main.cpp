#include <iostream>
#include <string>


#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>


// #include "atom.h"

#include "atom.h"
#include "monomer.h"
#include "polymer.h"
// #include "param.h"

int main() {

    using namespace polytop;

    RDKit::ROMol *mol1 = RDKit::PDBFileToMol("/Users/lily/workspace/synthezymes/01_structures/pdb/m1_u4.pdb");
    std::cout << mol1->getNumConformers() << std::endl;
    Monomer mol = Monomer(*mol1);
    std::cout << mol.atoms.size() << " monomer atoms " << std::endl;
    std::cout << mol.rdMol.getNumAtoms() << " rd monomer atoms " << std::endl;

    MonomerUnit unit = MonomerUnit(mol);
    std::cout << unit.atoms.size() << " unit atoms " << std::endl;
    std::cout << unit.rdMol.getNumAtoms() << " rd unit atoms " << std::endl;

    Polymer pol = Polymer("PolyAla");
    pol.addMonomerUnit(unit);
    std::cout << pol.atoms.size() << " polymer atoms " << std::endl;
    std::cout << pol.rdMol.getNumAtoms() << " rd polymer atoms " << std::endl;
    
    std::cout << pol.atoms[0]->name << std::endl;

    MonoPolyAtomIndexVect vec = { {43, 48}, {46, 51}, {44, 53}, {45, 52} };

    MonomerUnit unit2 = MonomerUnit(mol);
    std::cout << unit2.atoms.size() << " unit2 atoms " << std::endl;
    std::cout << unit2.rdMol.getNumAtoms() << " 2rd unit atoms " << std::endl;
    pol.addMonomerUnit(unit2, vec, true);

    std::cout << pol.atoms.size() << " polymer atoms " << std::endl;
    std::cout << pol.atoms[0]->name << std::endl;



    // RDKit::ROMol *mol1 = RDKit::SmilesToMol("Cc1ccccc1");
    
    // Monomer mol = Monomer(*mol1);
    // Polymer pol = Polymer("PolyAla");
    // // Atom atom = Atom(mol);
    // // mol.addAtom(atom);
    // std::cout << mol.resName << std::endl;
    // std::cout << mol.atoms.size() << std::endl;
    // std::cout << mol.atoms[0]->name << std::endl;
    // std::cout << mol.atoms[0]->bonds.data.size() << std::endl;
    // // mol.delAtom(atom);
    // // std::cout << mol.atoms.size() << std::endl;

    std::cout << "hmm" << std::endl;
}