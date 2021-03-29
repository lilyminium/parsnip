
from collections import UserList

import numpy as np

class ArrayList(UserList):

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
    
    @property
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
        self._data.extend(items)

    @uncache("_array")
    def remove(self, item):
        self._data.remove(item)