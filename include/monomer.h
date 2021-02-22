
#include <unordered_map>
#include <set>
#include <vector>
#include <string>
#include <functional>
#include <iostream>


#include <GraphMol/GraphMol.h>
#ifndef _PT_MONOMER_H
#define _PT_MONOMER_H


namespace polytop {
    class Atom;
    class Monomer {
        // friend class Polymer;

        public:
            Monomer(std::string name="UNK") {resName = name;};
            Monomer(const RDKit::ROMol &mol, std::string name="UNK");

            RDKit::RWMol rdMol;
            std::string resName = "UNK";
            std::string name = "UNK";
            std::vector<Atom*> atoms;
            int numAtoms;

            void setRDMol(RDKit::RWMol mol);

            int addAtom(Atom &atom);
            int delAtom(Atom &atom);
            void reindexAtoms(bool setOwningMol=true);
            void removeParamsWithinAtomSet(std::set<Atom*> atomSet);

            template <typename param_, typename paramvector_>
            std::set<param_*> getParamSet() {
                std::set<param_*> els;
                if (atoms.size() == 0) return els;
                paramvector_ pv = paramvector_(atoms[0]);
                for (auto &atom : atoms) {
                    paramvector_ pv_ = pv.getCorrespondingParamVector(atom);
                    for (auto &el : pv_) {
                        els.emplace(el);
                    }
                }
                return els;
            }

            template <typename param_, typename paramvector_>
            void copyParams(Monomer mol) {
                auto params = mol.getParamSet<param_, paramvector_>();
                for (auto &el : params) {
                    auto newEl = param_(*this, *el);
                    for (auto &at : newEl.atoms) {
                        auto paramList = paramvector_(at).getCorrespondingParamVector(at);
                        paramList.addParam(&newEl);
                    }
                }
            }
    };


    class MonomerUnit : public Monomer {

        public:
            MonomerUnit(Monomer mol);
            MonomerUnit(){};
            
    };
};

#endif