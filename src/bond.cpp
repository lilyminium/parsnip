
#include <array>

#include "atom.h"
#include "bond.h"
#include "monomer.h"

namespace polytop {

    Bond::Bond(Monomer mol, std::array<int, 2> indices) {
        atoms[0] = unit.atoms[indices[0]];
        atoms[1] = unit.atoms[indices[1]];
    }

    std::array<int, 2> Bond::getAtomIndices() {
        std::array<int, 2> indices = {atoms[0].index, atoms[1].index};
        return indices;
    }

    std::array<int, 2> Bond::updateAtom(Atom oldAtom, Atom newAtom) {
        for (size_t i = 0; i < atoms.size(); i++) {
            if (atoms[i] == oldAtom) {atoms[i] = newAtom; };
            return atoms;
        }
    }

    Bond Bond::copyTo(Monomer mol); {
        Bond bond = Bond(mol, getAtomIndices());
        for (auto& atom : bond.atoms) {
            atom.bonds.push_back(bond);
        }
        return bond;
    }



    
    
};