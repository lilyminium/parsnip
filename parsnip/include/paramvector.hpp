#include <array>
#include <vector>
#include <set>
#include <iostream>



#ifndef _PT_PARAMVECTOR_H
#define _PT_PARAMVECTOR_H

#include "param.hpp"

namespace polytop {

    template <typename T> class ParamVector {

        public:
            std::vector<T*> data;

            typename std::vector<T*>::iterator begin() {return data.begin();}
            typename std::vector<T*>::iterator end() {return data.end();}
            ParamVector<T> getCorrespondingParamVector(Atom *atom);
            std::size_t updateAtom(Atom *oldAtom, Atom *newAtom);
            std::size_t addParam(T *param);
            std::size_t removeParamsWithinAtomSet(std::set<Atom*>);
            std::size_t removeParamsNotInAtoms(std::vector<Atom*>);
            std::size_t size();
            void clear();

            T* &operator[](int i) {
                if ( i < 0 ) {
                    i = data.size() + i - 1;
                }
                return data[i];
            }
            
    };

    template <typename T>
    std::size_t ParamVector<T>::addParam(T *param) {
        data.emplace_back(param); //std::move(param));
        return data.size();
    };

    template <typename T>
    void ParamVector<T>::clear() { data.clear(); };

    template <typename T>
    std::size_t ParamVector<T>::size() { return data.size(); };
    
    template <typename T>
    std::size_t ParamVector<T>::updateAtom(Atom *oldAtom, Atom *newAtom) {
        for (auto &el : data) {
            el->updateAtom(oldAtom, newAtom);
        }
        return data.size();
    };

    template <typename T>
    std::size_t ParamVector<T>::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
        std::vector<std::size_t> toDelete;
        for (std::size_t i = 0; i < data.size(); i++) {
            if (data[i]->withinAtoms(atomSet)) {
                toDelete.insert(toDelete.begin(), i);
            }
        };

        for (auto &idx : toDelete) {
            data.erase(data.begin() + idx);
        };
        return data.size();
    };

    template <typename T>
    std::size_t ParamVector<T>::removeParamsNotInAtoms(std::vector<Atom*> atomVector) {
        std::vector<std::size_t> toDelete;
        for (std::size_t i = 0; i < data.size(); i++) {
            if (!data[i]->withinAtoms(atomVector)) {
                toDelete.insert(toDelete.begin(), i);
            }
        };

        for (auto &idx : toDelete) {
            data.erase(data.begin() + idx);
        };
        return data.size();
    }

    typedef ParamVector<Bond> BondVector;
    typedef ParamVector<Pair> PairVector;
    typedef ParamVector<Exclusion> ExclusionVector;
    typedef ParamVector<Angle> AngleVector;
    typedef ParamVector<Dihedral> DihedralVector;
    typedef ParamVector<Improper> ImproperVector;

};
#endif