#include <unordered_map>
#include <array>
#include <set>
#include <vector>
#include <string>
#include <functional>

#ifndef _PT_ATOM_H
#define _PT_ATOM_H

#include <GraphMol/GraphMol.h>
#include <GraphMol/MonomerInfo.h>

namespace polytop {

    class Atom {

        public:
            Atom(RDKit::Atom *atom, std::string name="X", double charge=0);

            int index = 0;
            double charge = 0;
            std::string name = "X";
            RDKit::Atom *rdAtom;
            RDKit::AtomPDBResidueInfo *resInfo = new RDKit::AtomPDBResidueInfo();
            void setResName(std::string resName);
            std::string getResName();
            void setResNum(unsigned int resNum);
            int getResNum();
            void setIndex(unsigned int atomIndex);
            unsigned int getIndex();
            unsigned int getSerial();
            void updateMonomerInfo();

            void setRDAtom(RDKit::Atom *atom);


    };
};

#endif