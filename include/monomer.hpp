
#include <unordered_map>
#include <set>
#include <vector>
#include <string>
#include <functional>
#include <iostream>


#include <GraphMol/GraphMol.h>
#ifndef _PT_MONOMER_H
#define _PT_MONOMER_H

#include "atom.hpp"
#include "tag.hpp"
#include "paramvector.hpp"

namespace polytop {
    class Monomer {
        // friend class Polymer;

        public:
            Monomer(std::string name="UNK") : name(name) {};
            Monomer(const RDKit::ROMol &mol, std::string name="UNK");

            RDKit::RWMol rdMol;
            std::string name;
            std::vector<Atom*> atoms;
            std::unordered_multimap<std::string, Tag*> tags;

            BondVector bonds = BondVector();
            PairVector pairs = PairVector();
            ExclusionVector exclusions = ExclusionVector();
            AngleVector angles = AngleVector();
            DihedralVector dihedrals = DihedralVector();
            ImproperVector impropers = ImproperVector();

            BondVector getParamVector(Bond *param) {return bonds;};
            AngleVector getParamVector(Angle *param) {return angles;};
            DihedralVector getParamVector(Dihedral *param) {return dihedrals;};
            ImproperVector getParamVector(Improper *param) {return impropers;};
            PairVector getParamVector(Pair *param) {return pairs;};
            ExclusionVector getParamVector(Exclusion *param) {return exclusions;};

            BondVector getParamVector(BondVector param) {return bonds;};
            AngleVector getParamVector(AngleVector param) {return angles;};
            DihedralVector getParamVector(DihedralVector param) {return dihedrals;};
            ImproperVector getParamVector(ImproperVector param) {return impropers;};
            PairVector getParamVector(PairVector param) {return pairs;};
            ExclusionVector getParamVector(ExclusionVector param) {return exclusions;};

            void setRDMol(RDKit::RWMol mol);
            RDKit::RWMol* rdMolFromAtoms(std::size_t numNeighbors=3);
            void initParams();

            std::vector<unsigned int> getAtomIndices();

            int addAtom(Atom &atom);
            std::size_t delAtomByIndex(unsigned int index);
            std::size_t removeAtomByIndex(unsigned int index);
            std::size_t delAtomFromVector(Atom *atom);
            std::size_t delAtomsFromVector(std::vector<Atom*> toDelete);
            std::size_t delAtomsFromVector(std::set<Atom*> toDelete);
            std::size_t removeAtom(Atom *atom);
            void reindexAtoms();
            void replaceAtom(Atom *oldAtom, Atom *newAtom);
            void removeParamsWithinAtomSet(std::set<Atom*> atomSet);
            void removeParamsWithinAtomSetFromVector(std::set<Atom*> atomSet);
            void cleanParams();

            void copyParamsFrom(Monomer* other);
            void addParamsFromMolecule(Monomer* other);
            bool containsAtom(Atom* atom);
            bool containsAllAtoms(std::vector<Atom*> atomVector);
            
            std::size_t addTag(std::string name, std::vector<Atom*> atomVector);
            std::size_t addTag(std::string name, std::vector<unsigned int> indices);
            std::size_t addTag(Tag* tag);
            std::size_t safeAddTag(Tag* tag);
            std::size_t unsafeAddTagsFrom(Monomer* other);
            std::size_t copyTagsFrom(Monomer* other);
            std::size_t addTagsFrom(Monomer* other, bool checkContainsAtoms=true);
            std::size_t safeAddTagsFrom(Monomer* other);
            std::size_t cleanTags();

            Tag* getFirstTagByName(std::string name);
            Tag* getLastTagByName(std::string name);
            std::size_t copyTagFrom(Tag* tag);
            std::vector<Atom*> getAtomsByIndices(std::vector<unsigned int> indices);

            template <typename T>
            std::size_t addParam(T* param) {
                auto vec = getParamVector(param);
                vec.addParam(param);
                return vec.data.size();
            }

            template <typename T>
            std::size_t addParamFrom(T* param) {
                auto vec = getParamVector(param);
                vec.addParam(param);
                return vec.data.size();
            }

            template <typename T>
            std::size_t copyParamFrom(T* param) {
                auto vec = getParamVector(param);
                auto indices = param->getAtomIndices();
                auto atoms = getAtomsByIndices(indices);
                auto newParam = new T(atoms);
                vec.addParam(newParam);
                return vec.data.size();
            }

            template <typename T>
            std::size_t addParamsFromVector(T paramvector) {
                auto vec = getParamVector(paramvector);
                for (auto &el : paramvector) {
                    addParamFrom(el);
                }
                return vec.data.size();
            }

            template <typename T>
            std::size_t copyParamsFrom(T paramvector) {
                auto vec = getParamVector(paramvector);
                for (auto &el : paramvector) {
                    copyParamFrom(el);
                }
                return vec.data.size();
            }

            template <unsigned long nAtoms>
            std::array<Atom*, nAtoms> getAtomsByIndices(std::array<unsigned int, nAtoms> indices) {
                std::array<Atom*, nAtoms> atomArray;
                for (std::size_t i = 0; i < indices.size(); i++) {
                    atomArray[i] = atoms[indices[i]];
                }
                return atomArray;
            }

            // template <unsigned long nAtoms>
            // std::array<Atom*, nAtoms> getAtomsByIndices(std::array<int, nAtoms> indices) {
            //     std::array<Atom*, nAtoms> atomArray;
            //     for (std::size_t i = 0; i < indices.size(); i++) {
            //         atomArray[i] = atoms[indices[i]];
            //     }
            //     return atomArray;
            // }




            // template <typename param_, typename paramvector_>
            // std::set<param_*> getParamSet() {
            //     std::set<param_*> els;
            //     if (atoms.size() == 0) return els;
            //     paramvector_ pv = paramvector_(atoms[0]);
            //     for (auto &atom : atoms) {
            //         paramvector_ pv_ = pv.getCorrespondingParamVector(atom);
            //         for (auto &el : pv_) {
            //             els.emplace(el);
            //         }
            //     }
            //     return els;
            // }

            // template <typename param_, typename paramvector_>
            // void copyParams(Monomer mol) {
            //     auto params = mol.getParamSet<param_, paramvector_>();
            //     for (auto &el : params) {
            //         auto newEl = param_(*this, *el);
            //         for (auto &at : newEl.atoms) {
            //             auto paramList = paramvector_(at).getCorrespondingParamVector(at);
            //             paramList.addParam(&newEl);
            //         }
            //     }
            // }
    };


    class MonomerUnit : public Monomer {

        public:
            MonomerUnit(Monomer *mol);
            MonomerUnit(){};
            
    };
};

#endif