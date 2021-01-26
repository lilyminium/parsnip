#include <iostream>
#include <string>

// #include "atom.h"

#include "atom.h"
#include "monomer.h"
#include "polymer.h"
// #include "param.h"

int main() {

    using namespace polytop;
    
    Monomer mol = Monomer("ALA");
    Atom atom = Atom(mol);
    mol.addAtom(atom);
    std::cout << mol.resName << std::endl;
    std::cout << mol.atoms.size() << std::endl;
    std::cout << mol.atoms[0]->name << std::endl;
    mol.delAtom(atom);
    std::cout << mol.atoms.size() << std::endl;
}