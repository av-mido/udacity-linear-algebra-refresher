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
        squares = [x**Decimal(2.0) for x in self.coordinates]
        mag_squared = sum(squares) # reduce(lambda x,y:x+y, squares)
        return Decimal(math.sqrt(mag_squared))

    def normalization(self):
        mag = self.magnitude()
        # return Vector([x/mag for x in self.coordinates])
        return self.__mul__(Decimal('1.0')/mag)

    def dot(self, v):
        return Decimal(sum([a*b for a,b in zip(self.coordinates, v.coordinates)]))

    def angle_rad(self, v):
        # arg = self.dot(v) / (self.magnitude() * v.magnitude())
        arg = self.normalization().dot(v.normalization())
        arg_bounded = min(1,max(arg,-1))
        return math.acos(arg_bounded)

    def angle_degrees(self, v):
        rads = self.angle_rad(v)
        return math.degrees(rads)

    def is_parallel(self, v, tolerance=0.001):
        return (self.is_zero_vector() or
                v.is_zero_vector() or
                (abs(self.angle_degrees(v)) < tolerance) or
                (180 + tolerance > self.angle_degrees(v) > 180-tolerance)
                )

    def is_zero_vector(self, tolerance=0.001):
        return self.magnitude() < tolerance

    def is_orthogonal(self, v, tolerance=0.001):
        return abs(self.dot(v)) < tolerance

    def projection(self, b):
        b_normalized = b.normalization()
        return self.dot(b_normalized) * b_normalized

    def projection_orthogonal(self, b):
        proj = self.projection(b)
        return self - proj

    def cross_product(self, w):
        # r1 = self.coordinates[1] * w.coordinates[2] - w.coordinates[1] * self.coordinates[2]
        # r2 = -1 * self.coordinates[0] * w.coordinates[2] + w.coordinates[0] * self.coordinates[2]
        # r3 = self.coordinates[0] * w.coordinates[1] - w.coordinates[0] * self.coordinates[1]
        x1, y1, z1 = self.coordinates
        x2, y2, z2 = w.coordinates
        r_coords = [y1*z2 - y2*z1, -1*x1*z2 + x2*z1, x1*y2 - x2*y1]
        return Vector(r_coords)

    def area_parallelogram(self, w):
        r = self.cross_product(w)
        return r.magnitude()

    def area_half_triangle(self, w):
        return self.area_parallelogram(w) / Decimal(2.0)

    def __iter__(self):
        return iter(self.coordinates)

    def __getitem__(self, key):
        return self.coordinates[key]