"""
unitcell.py: Routines for converting between different representations of the periodic unit cell.

"""
from __future__ import absolute_import


from NexaBind import Vec3
from NexaBind.unit import nanometers, is_quantity, norm, dot, radians
import math


def computePeriodicBoxVectors(a_length, b_length, c_length, alpha, beta, gamma):
    """Convert lengths and angles to periodic box vectors.
    
    Lengths should be given in nanometers and angles in radians (or as Quantity
    instances)
    """

    if is_quantity(a_length): a_length = a_length.value_in_unit(nanometers)
    if is_quantity(b_length): b_length = b_length.value_in_unit(nanometers)
    if is_quantity(c_length): c_length = c_length.value_in_unit(nanometers)
    if is_quantity(alpha): alpha = alpha.value_in_unit(radians)
    if is_quantity(beta): beta = beta.value_in_unit(radians)
    if is_quantity(gamma): gamma = gamma.value_in_unit(radians)

    # Compute the vectors.

    a = [a_length, 0, 0]
    b = [b_length*math.cos(gamma), b_length*math.sin(gamma), 0]
    cx = c_length*math.cos(beta)
    cy = c_length*(math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma)
    cz = math.sqrt(c_length*c_length-cx*cx-cy*cy)
    c = [cx, cy, cz]

    # If any elements are very close to 0, set them to exactly 0.

    for i in range(3):
        if abs(a[i]) < 1e-6:
            a[i] = 0.0
        if abs(b[i]) < 1e-6:
            b[i] = 0.0
        if abs(c[i]) < 1e-6:
            c[i] = 0.0
    a = Vec3(*a)
    b = Vec3(*b)
    c = Vec3(*c)

    # Make sure they're in the reduced form required by NexaBind.

    c = c - b*round(c[1]/b[1])
    c = c - a*round(c[0]/a[0])
    b = b - a*round(b[0]/a[0])
    return (a, b, c)*nanometers

def reducePeriodicBoxVectors(periodicBoxVectors):
    """ Reduces the representation of the PBC. periodicBoxVectors is expected to
    be an unpackable iterable of length-3 iterables
    """
    if is_quantity(periodicBoxVectors):
        a, b, c = periodicBoxVectors.value_in_unit(nanometers)
    else:
        a, b, c = periodicBoxVectors
    a = Vec3(*a)
    b = Vec3(*b)
    c = Vec3(*c)

    c = c - b*round(c[1]/b[1])
    c = c - a*round(c[0]/a[0])
    b = b - a*round(b[0]/a[0])

    return (a, b, c) * nanometers

def computeLengthsAndAngles(periodicBoxVectors):
    """Convert periodic box vectors to lengths and angles.

    Lengths are returned in nanometers and angles in radians.
    """
    if is_quantity(periodicBoxVectors):
        (a, b, c) = periodicBoxVectors.value_in_unit(nanometers)
    else:
        a, b, c = periodicBoxVectors
    a_length = norm(a)
    b_length = norm(b)
    c_length = norm(c)
    alpha = math.acos(dot(b, c)/(b_length*c_length))
    beta = math.acos(dot(c, a)/(c_length*a_length))
    gamma = math.acos(dot(a, b)/(a_length*b_length))
    return (a_length, b_length, c_length, alpha, beta, gamma)
