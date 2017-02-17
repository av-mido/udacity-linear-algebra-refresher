import math
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def __add__(self, v):
        return Vector(tuple(map(sum, zip(self.coordinates, v.coordinates))))

    def __sub__(self, v):
        return Vector([a-b for a,b in zip(self.coordinates, v.coordinates)])

    def __mul__(self, v_or_s):
        if hasattr(v_or_s, "__len__"):
            raise ValueError("Vector multiplication not implemented yet.")
        else:
            # return Vector(tuple(map(lambda x: v_or_s * x, list(self.coordinates))))
            return Vector([Decimal(v_or_s) * x for x in self.coordinates])

    def __rmul__(self, v_or_s):
        if hasattr(v_or_s, "__len__"):
            raise ValueError("Vector multiplication not implemented yet.")
        else:
            # return Vector(tuple(map(lambda x: v_or_s * x, list(self.coordinates))))
            return Vector([Decimal(v_or_s) * x for x in self.coordinates])

    def magnitude(self):
        squares = [x**2 for x in self.coordinates]
        mag_squared = sum(squares) # reduce(lambda x,y:x+y, squares)
        return Decimal(math.sqrt(mag_squared))

    def normalization(self):
        mag = self.magnitude()
        # return Vector([x/mag for x in self.coordinates])
        return self.__mul__(Decimal('1.0')/mag)

    def dot(self, v):
        return sum([a*b for a,b in zip(self.coordinates, v.coordinates)])

    def angle_rad(self, v):
        # arg = self.dot(v) / (self.magnitude() * v.magnitude())
        arg = self.normalization().dot(v.normalization())
        return math.acos(arg)

    def angle_degrees(self, v):
        rads = self.angle_rad(v)
        return math.degrees(rads)