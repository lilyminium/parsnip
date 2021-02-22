#include <array>
#include <vector>
#include <set>


#ifndef _PT_PARAMVECTOR_H
#define _PT_PARAMVECTOR_H

#include "param.h"

namespace polytop {
    class Atom;
    class Monomer;
    class MonomerUnit;

    template <typename T> class ParamVector {
        friend class Atom;

        public:
            ParamVector(Atom *atom);
            typename std::vector<T*>::iterator begin() {return data.begin();}
            typename std::vector<T*>::iterator end() {return data.end();}
            ParamVector<T> getCorrespondingParamVector(Atom *atom);
            std::size_t updateAtom(Atom &newAtom);
            std::size_t addParam(T *param);
            void copyParamsTo(Monomer &mol);
            void copyParamsTo(MonomerUnit &mol);
            std::size_t removeParamsWithinAtomSet(std::set<Atom*>);

            std::vector<T*> data;
            Atom *atom;
            

            
    };

    template <typename T>
    ParamVector<T>::ParamVector(Atom *newAtom) {
        atom = newAtom;
    }

    template <typename T>
    std::size_t ParamVector<T>::addParam(T *param) {
        data.push_back(std::move(param));
        return data.size();
    };


    template <typename T>
    std::size_t ParamVector<T>::updateAtom(Atom &newAtom) {
        std::vector<std::size_t> toDelete;
        ParamVector<T> vec = ParamVector<T>::getCorrespondingParamVector(&newAtom);
        for (std::size_t i = 0; i < data.size(); i++) {
            data[i]->updateAtom(atom, &newAtom);
            if (data[i]->containsAtom(&newAtom)) {
                toDelete.insert(toDelete.begin(), i);
                vec.addParam(data[i]);
            };
        };
        for (auto &idx : toDelete) {
            data.erase(data.begin() + idx);
        };
        return data.size();
    };

    template <typename T>
    std::size_t ParamVector<T>::removeParamsWithinAtomSet(std::set<Atom*> atomSet) {
        std::vector<std::size_t> toDelete;
        for (std::size_t i = 0; i < data.size(); i++) {
            if (data[i]->withinAtomSet(atomSet)) {
                toDelete.insert(toDelete.begin(), i);
            }
        };
        for (auto &idx : toDelete) {
            data.erase(data.begin() + idx);
        };
        return data.size();
    };

    typedef ParamVector<Bond> BondVector;
    typedef ParamVector<Pair> PairVector;
    typedef ParamVector<Exclusion> ExclusionVector;
    typedef ParamVector<Angle> AngleVector;
    typedef ParamVector<Dihedral> DihedralVector;
    typedef ParamVector<Improper> ImproperVector;

};
#endif