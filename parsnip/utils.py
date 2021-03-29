import functools


def isiterable(obj):
    """
    If an object is iterable and not a string
    """
    if isinstance(obj, str):
        return False
    if hasattr(obj, 'next'):
        return True
    try:
        len(obj)
    except (TypeError, AttributeError):
        return False
    return True

def asiterable(obj):
    if not isiterable(obj):
        obj = [obj]
    return obj



def cached(func):
    """Cache a property within a class.

    Requires the Class to have a cache dict called ``_cache``.

    Example
    -------
    How to add a cache for a variable to a class by using the `@cached`
    decorator::

       class A(object):
           def__init__(self):
               self._cache = dict()

           @property
           @cached('keyname')
           def size(self):
               # This code gets run only if the lookup of keyname fails
               # After this code has been ran once, the result is stored in
               # _cache with the key: 'keyname'
               return size

    .. note::

        Copied from MDAnalysis.
    """

    key = func.__name__
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        try:
            return self._cache[key]
        except KeyError:
            self._cache[key] = ret = func(self, *args, **kwargs)
            return ret
    return wrapper

def uncache(*keys):
    def refresh_cache(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            for k in keys:
                self._cache.pop(k, None)
            return func(self, *args, **kwargs)
        return wrapper
    return refresh_cache