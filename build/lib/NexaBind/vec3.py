"""
vec3.py: Defines the Vec3 class used by NexaBind

"""
from __future__ import absolute_import, division


from . import unit
from collections import namedtuple

class Vec3(namedtuple('Vec3', ['x', 'y', 'z'])):
    """Vec3 is a 3-element tuple that supports many math operations."""

    def __new__(cls, x, y, z):
        """Create a new Vec3."""
        return tuple.__new__(cls, (x, y, z))

    def __getnewargs__(self):
        "Support for pickle protocol 2: http://docs.python.org/2/library/pickle.html#pickling-and-unpickling-normal-class-instances"
        return self[0], self[1], self[2]

    def __add__(self, other):
        """Add two Vec3s."""
        return Vec3(self.x+other[0], self.y+other[1], self.z+other[2])

    def __radd__(self, other):
        """Add two Vec3s."""
        return Vec3(self.x+other[0], self.y+other[1], self.z+other[2])

    def __sub__(self, other):
        """Add two Vec3s."""
        return Vec3(self.x-other[0], self.y-other[1], self.z-other[2])

    def __rsub__(self, other):
        """Add two Vec3s."""
        return Vec3(other[0]-self.x, other[1]-self.y, other[2]-self.z)

    def __mul__(self, other):
        """Multiply a Vec3 by a constant."""
        if unit.is_unit(other):
            return unit.Quantity(self, other)
        return Vec3(other*self.x, other*self.y, other*self.z)

    def __rmul__(self, other):
        """Multiply a Vec3 by a constant."""
        if unit.is_unit(other):
            return unit.Quantity(self, other)
        return Vec3(other*self.x, other*self.y, other*self.z)

    def __div__(self, other):
        """Divide a Vec3 by a constant."""
        return Vec3(self.x/other, self.y/other, self.z/other)
    __truediv__ = __div__

    def __deepcopy__(self, memo):
        return Vec3(self.x, self.y, self.z)

    def __neg__(self):
        return Vec3(-self.x, -self.y, -self.z)
