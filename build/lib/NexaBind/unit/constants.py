#!/bin/env python
"""
Module NexaBind.unit.constants

"""
from __future__ import print_function, division, absolute_import


from .unit_definitions import *

#################
### CONSTANTS ###
#################

# codata 2018
AVOGADRO_CONSTANT_NA = 6.02214076e23 / mole
BOLTZMANN_CONSTANT_kB = 1.380649e-23 * joule / kelvin
MOLAR_GAS_CONSTANT_R = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB
SPEED_OF_LIGHT_C = 2.99792458e8 * meter / second
GRAVITATIONAL_CONSTANT_G = 6.6743e-11 * newton * meter**2 / kilogram**2

# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])
