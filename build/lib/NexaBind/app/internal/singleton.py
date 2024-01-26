"""
Creates a subclass for all classes intended to be a singleton. This
maintains the correctness of instance is instance even following
pickling/unpickling
"""
class Singleton(object):
    _inst = None
    def __new__(cls):
        if cls._inst is None:
            cls._inst = super(Singleton, cls).__new__(cls)
        #Omid: If _inst already exists (is not None), it returns the existing instance, ensuring only one instance exists throughout the program's lifetime.
        return cls._inst
    #Note24: enforces the creation of a single instance and maintains the instance's correctness even when serialized (pickled) and deserialized (unpickled).
    def __reduce__(self):
        return repr(self)