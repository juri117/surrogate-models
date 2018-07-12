# ==============================================================================
# description     :this file is a helper to enable singletons in python... meaning, a class that uses Singleton as
#                  metaclass will only allow one instance, further tries to create an instance just returns the first
#                  instances pointer
# author          :Juri Bieler
# date            :2018-02-22
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

#helper class to support Singletons
class Singleton(type):
    #Define an Instance operation that lets clients access its unique instance.
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

    '''
    def __init__(cls, name, bases, attrs, **kwargs):
        super().__init__(name, bases, attrs)
        cls._instance = None

    def __call__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__call__(*args, **kwargs)
        return cls._instance
    '''