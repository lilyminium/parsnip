from collections import UserDict, defaultdict

import numpy as np

from ..utils import uncache, cached, asiterable
from .mixins import TypeMimicMixin, ContainsIndicesMixin, CachedIndicesMixin


class TagDict(TypeMimicMixin, ArrayDict):
    def __init__(self, *args, **kwargs):
        self._cache = {}
        self._data = defaultdict(list)
        if args:
            if len(args) == 1:
                self.extend(dict(args[0]))
            else:
                raise TypeError("ArrayDict expected at most 1 argument, got "
                                f"{len (args)}")
        self.extend(kwargs)

    def __contains__(self, item):
        return item in self._array or item in self._data.keys()

    @cached
    def _array(self):
        return np.array(list(self.itertags()))

    def itertags(self):
        for values in self.values():
            for item in values:
                yield item


    # ========= mimic dict =========

    def get(self, key, *default):
        return self._data.get(key, *default)

    @uncache("_array", "_indices")
    def update(self, other):
        self._data.update(other)
    
    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    @uncache("_array", "_indices")
    def extend(self, other):
        for key, value in other.items():
            self._data[key].extend(asiterable(other))

    def all_in_indices(self, indices):
        if len(self._indices) == 0:
            return type(self)()
        arr = np.isin(self._indices, indices)
        keys = np.array(tuple(self._data.keys()))
        keep = keys[np.all(arr, axis=-1)]
        return type(self)({k: self._data[k] for k in keep})

    def any_in_indices(self, indices):
        if len(self._indices) == 0:
            return type(self)()
        arr = np.isin(self._indices, indices)
        keys = np.array(tuple(self._data.keys()))
        keep = keys[np.any(arr, axis=-1)]
        return type(self)({k: self._data[k] for k in keep})

    def xor_in_indices(self, indices):
        if len(self._indices) == 0:
            return type(self)()
        arr = np.isin(self._indices, indices)
        keys = np.array(tuple(self._data.keys()))
        mat = np.any(arr, axis=-1) & ~np.all(arr, axis=-1)
        keep = keys[mat]
        return type(self)({k: self._data[k] for k in keep})
