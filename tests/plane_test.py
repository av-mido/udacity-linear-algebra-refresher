from plane import Plane
from vector import Vector

p1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46)
p2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65)
print("planes are parallel? {}".format(p1.is_parallel_to(p2)))
print("planes are equal? {}".format(p1 == p2))


def test_parallel():
    assert p1.is_parallel_to(p2)


def test_equal():
    assert p1 == p2
