#!/bin/env python
"""
Module NexaBind.unit.standard_dimensions

Definition of principal dimensions: mass, length, time, etc.

"""
from __future__ import division, print_function, absolute_import


from .basedimension import BaseDimension

##################
### DIMENSIONS ###
##################

mass_dimension = BaseDimension('mass')
length_dimension = BaseDimension('length')
time_dimension = BaseDimension('time')
temperature_dimension = BaseDimension('temperature')
amount_dimension = BaseDimension('amount')
charge_dimension = BaseDimension('charge')
luminous_intensity_dimension = BaseDimension('luminous intensity')
angle_dimension = BaseDimension('angle')
information_dimension = BaseDimension('information')


# run module directly for testing
if __name__=='__main__':
    # Test the examples in the docstrings
    import doctest, sys
    doctest.testmod(sys.modules[__name__])

