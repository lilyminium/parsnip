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

           @cached('keyname')
           def size(self):
               # This code gets run only if the lookup of keyname fails
               # After this code has been ran once, the result is stored in
               # _cache with the key: 'keyname'
               return size

    .. note::

        Adapted from MDAnalysis.
    """

    key = func.__name__
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        try:
            return self._cache[key]
        except KeyError:
            self._cache[key] = ret = func(self, *args, **kwargs)
            return ret
    return property(wrapper)

def uncache(*keys):
    def refresh_cache(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            for k in keys:
                self._cache.pop(k, None)
            return func(self, *args, **kwargs)
        return wrapper
    return refresh_cache


# std::vector<unsigned int> splitIntegerIntoRatio(unsigned int total,
#                                                     std::vector<double> ratio) {
#         double ratioTotal = std::accumulate(ratio.begin(), ratio.end(),
#                                             decltype(ratio)::value_type(0));

#         std::vector<double> currentFraction;
#         std::vector<double> fraction;
#         fraction.reserve(ratio.size());
#         currentFraction.reserve(ratio.size());

#         for (auto r : ratio) {
#             fraction.push_back(r/ratioTotal);
#             currentFraction.push_back(0);
#         }

        
#         std::vector<unsigned int> currentRatio;
#         double smallest = *std::min_element(ratio.begin(), ratio.end());
#         for (auto &el : ratio) {
#             currentRatio.push_back(el/smallest);
#         }
#         unsigned int currentTotal = std::accumulate(currentRatio.begin(), currentRatio.end(),
#                                                     decltype(currentRatio)::value_type(0));
        
#         int quotient = (total / currentTotal) - 1;
#         if (quotient) {
#             for (std::size_t i = 0; i < currentRatio.size(); i++) {
#                 currentRatio[i] *= quotient;
#             };
#         };

#         currentTotal = std::accumulate(currentRatio.begin(), currentRatio.end(),
#                                        decltype(currentRatio)::value_type(0));
        
#         while (currentTotal < total) {
#             for (std::size_t i = 0; i < currentRatio.size(); i++ ) {
#                 currentFraction[i] = double(currentRatio[i]) / total;
#                 currentFraction[i] /= fraction[i];
#             };

#             auto it = std::min_element(currentFraction.begin(), currentFraction.end());
#             unsigned int minIndex = std::distance(currentFraction.begin(), it);
#             currentRatio[minIndex] += 1;
#             currentTotal = std::accumulate(currentRatio.begin(), currentRatio.end(),
#                                            decltype(currentRatio)::value_type(0));
#         }

#         return currentRatio;
#     };
