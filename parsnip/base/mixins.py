
import numpy as np

from ..utils import uncache, cached


class TypeMimicMixin:

    def __contains__(self, item):
        return item in self._array

    def __repr__(self):
        return repr(self._data)

    def __iter__(self):
        for item in self._data:
            yield item

    def __len__(self):
        return len(self._data)

    @uncache("_array")
    def __delitem__(self, key):
        del self._data[key]

    @uncache("_array")
    def __setitem__(self, key, item):
        self._data[key] = item

    @uncache("_array")
    def pop(self, key, *default):
        return self._data.pop(key, *default)

    def __getitem__(self, key):
        return self._data[key]

class CachedIndicesMixin(TypeMimicMixin):

    @uncache("_indices")
    def __setitem__(self, key, item):
        super().__setitem__(key, item)

    @uncache("_indices")
    def __delitem__(self, key):
        super().__delitem__(key)

    @uncache("_indices")
    def pop(self, key, *default):
        super().pop(key, *default)

    @cached
    def _indices(self):
        return np.array([x.indices for x in self._array])


class ContainsIndicesMixin:

    def all_in_indices(self, indices):
        if not self._data:
            return type(self)([])
        arr = np.isin(self._indices, indices)
        new = type(self)(self._array[np.all(arr, axis=-1)])
        return new

    def any_in_indices(self, indices):
        if not self._data:
            return type(self)([])
        arr = np.isin(self._indices, indices)
        new = type(self)(self._array[np.any(arr, axis=-1)])
        return new

    def xor_in_indices(self, indices):
        if not self._data:
            return type(self)([])
        arr = np.isin(self._indices, indices)
        matrix = np.any(arr, axis=-1) & ~np.all(arr, axis=-1)
        new = type(self)(self._array[matrix])
        return new
