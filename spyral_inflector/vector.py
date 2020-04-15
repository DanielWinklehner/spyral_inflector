import numpy as np

__author__ = "Daniel Winklehner"
__doc__ = "Simple vector class covering basic vector operations."


class Vector(object):
    def __init__(self, vector):

        assert(len(vector) in [2, 3]), "For now we can only handle vectors of length 2 or 3"

        self._vector = np.array(vector)
        self._dim = len(vector)
        self._length = self.calculate_length()

    def __str__(self):
        return str(self.vector)

    def __getitem__(self, item):
        return self._vector[item]

    def __mul__(self, other):
        return np.dot(self._vector, other[:])

    def __rmul__(self, other):
        return np.dot(self._vector, other[:])

    def __add__(self, other):
        try:

            return Vector(vector=(self._vector + other.vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    def __radd__(self, other):
        try:

            return Vector(vector=(self._vector + other.vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    def __sub__(self, other):
        try:

            return Vector(vector=(self._vector - other.vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    def __rsub__(self, other):
        try:

            return Vector(vector=(other.vector - self._vector))

        except AttributeError:

            print("Can only add and subtract vectors from vectors!")

    @property
    def length(self):
        return self._length

    @property
    def dim(self):
        return self._dim

    @property
    def vector(self):
        return self._vector

    def calculate_length(self):
        return np.sqrt(self * self)

    def normalized(self):
        return Vector(vector=(self._vector / self._length))

    def angle_with(self, other):
        return np.arccos(np.dot(self._vector, other[:])/(self._length * other.length))

    def cross(self, other):
        return Vector(vector=np.cross(self._vector, other.vector))


if __name__ == '__main__':
    # Test the vector class
    v1 = Vector([1, 0, 0])
    v2 = Vector([0, 1, 0])

    print(v1 + v2)
    print(v1 - v2)
    print(v1 * v2)
    print(v1.length)
    print(np.rad2deg(v1.angle_with(v2)))
    print((v1 + v2).normalized())
