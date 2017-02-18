from vector import Vector
from line import Line

print("quiz coding functions for lines")
v = Vector([4.046, 2.836])
l1 = Line(normal_vector=v, constant_term=1.21)
l2 = Line(Vector([10.115, 7.09]), 3.025)
print("l1: {}".format(l1))
print("l1 basepoint: {}".format(l1.basepoint))
print("l2: {}".format(l2))
print("l2 basepoint: {}".format(l2.basepoint))
print("parallel: {}".format(l1.is_parallel_with(l2)))
print("is equal to: {}".format(l1 == l2))
print("intersection: {}".format(l1.find_intersection(l2)))

l1 = Line(Vector([7.204, 3.182]), constant_term=8.68)
l2 = Line(Vector([8.172, 4.114]), constant_term=9.883)
print("intersection 2: {}".format(l1.find_intersection(l2)))

l1 = Line(Vector([1.182, 5.562]), 6.744)
l2 = Line(Vector([1.773, 8.343]), 9.525)
print("intersection 3: {}".format(l1.find_intersection(l2)))