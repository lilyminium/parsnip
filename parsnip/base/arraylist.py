from collections import UserList

import numpy as np

from ..utils import cached, uncache
from .mixins import TypeMimicMixin, ContainsIndicesMixin, CachedIndicesMixin


class ArrayList(TypeMimicMixin, UserList):
    """Pretend to be a mutable array"""

    def __init__(self, items=None):
        self._cache = {}
        self._data = []
        if items is not None:
            self._data.extend(items)

    def __getitem__(self, key):
        val = self._array[key]
        if not isinstance(val, np.ndarray):
            return val
        return type(self)(val)

    def copy(self):
        return type(self)(self._data)


    @cached
    def _array(self):
        return np.array(self._data)

    # ========= mimic set =========

    def __sub__(self, other):
        return type(self)(self._array[~np.isin(self._array, other)])

    # ========= mimic list =========

    def __add__(self, other):
        return type(self)(list(self._data) + list(other._data))
    
    def __radd__(self, other):
        if not other:
            return type(self)(self._data)
        return other.__add__(self)

    @uncache("_array")
    def append(self, item):
        self._data.append(item)

    @uncache("_array")
    def extend(self, items):
        self._data.extend(list(items))

    @uncache("_array")
    def remove(self, item):
        self._data.remove(item)

    @uncache("_array")
    def clear(self):
        self._data = []


class ArrayListWithIndices(ContainsIndicesMixin, CachedIndicesMixin, ArrayList):
    """Pretend to be a mutable array of items with indices"""

    append = uncache("_indices",)(ArrayList.append)
    extend = uncache("_indices")(ArrayList.extend)
    remove = uncache("_indices")(ArrayList.remove)
    clear = uncache("_indices")(ArrayList.clear)
    

    @property
    def indices(self):
        return self._indices


class ParamList(ArrayListWithIndices):
    append = uncache("_atoms", "_indices")(ArrayList.append)
    remove = uncache("_atoms", "_indices")(ArrayList.remove)
    clear = uncache("_atoms", "_indices")(ArrayList.clear)

    @cached
    def _atoms(self):
        return np.array([x.atoms for x in self._data])

    @uncache("_array", "_atoms", "_indices")
    def extend(self, items):
        try:
            items = items._data
        except AttributeError:
            pass
        self._data.extend(items)

    @property
    def atoms(self):
        return self._atoms

    @uncache("_atoms", "_indices", "_array")
    def remove_within_atoms(self, atoms):
        if not len(self._data):
            return
        mask = np.all(np.isin(self.atoms, atoms), axis=-1)
        self._data = self._array[~mask].tolist()

    @uncache("_atoms", "_indices", "_array")
    def keep_within_atoms(self, atoms):
        if not len(self._data):
            return
        mask = np.all(np.isin(self.atoms, atoms), axis=-1)
        self._data = self._array[mask].tolist()
