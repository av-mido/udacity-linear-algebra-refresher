import vector

v1 = vector.Vector([1,2,3])
v2 = vector.Vector([4,5,8])
v3 = vector.Vector([5,0,0])
print(v1)
print(v2)

print("\ntrying adding")
print(v1 + v2)
print("")

print("trying subtracting")
print(v1 - v2)
print("")

print("scalar mult")
print(v1*5)
print(5 * v1)

print("")
print("magnitude")
print(v1.magnitude())
print(v3.magnitude())

print("")
print("normalization")
print(v1.normalization())
print(v3.normalization())

print("")
print("dot product")
print(v1.dot(v2))

print("")
print("angle radians")
va = vector.Vector([3.183, -7.627])
vb = vector.Vector([-2.668, 5.319])
print(va.angle_rad(vb))

print("")
print("angle degrees")
va = vector.Vector([7.35, 0.221, 5.188])
vb = vector.Vector([2.751, 8.259, 3.985])
print(va.angle_degrees(vb))

print("Is parallel")
va = vector.Vector([-7.579, -7.88])
vb = vector.Vector([22.737, 23.64])
print(va.is_parallel(vb))
va = vector.Vector([-2.029, 9.97, 4.172])
vb = vector.Vector([-9.231, -6.639, -7.245])
print(va.is_parallel(vb))
va = vector.Vector([-2.328, -7.284, -1.214])
vb = vector.Vector([-1.821, 1.072, -2.94])
print(va.is_parallel(vb))
va = vector.Vector([2.118, 4.827])
vb = vector.Vector([0, 0])
print(va.is_parallel(vb))

print("Is orthogonal")
va = vector.Vector([-7.579, -7.88])
vb = vector.Vector([22.737, 23.64])
print(va.is_orthogonal(vb))
va = vector.Vector([-2.029, 9.97, 4.172])
vb = vector.Vector([-9.231, -6.639, -7.245])
print(va.is_orthogonal(vb))
va = vector.Vector([-2.328, -7.284, -1.214])
vb = vector.Vector([-1.821, 1.072, -2.94])
print(va.is_orthogonal(vb))
va = vector.Vector([2.118, 4.827])
vb = vector.Vector([0, 0])
print(va.is_orthogonal(vb))

print("")
print("projection")
v = vector.Vector([3.039, 1.879])
b = vector.Vector([0.825, 2.036])
print(v.projection(b))

print("")
print("orthogonal projection")
v = vector.Vector([-9.88, -3.264, -8.159])
b = vector.Vector([-2.155, -9.353, -9.473])
print(v.projection_orthogonal(b))

print("")
print("cross product")
v = vector.Vector([8.462, 7.893, -8.187])
w = vector.Vector([6.984, -5.975, 4.778])
print(v.cross_product(w))
v = vector.Vector([-8.987, -9.838, 5.031])
w = vector.Vector([-4.268, -1.861, -8.866])
print(v.area_parallelogram(w))
print(v.area_half_triangle(w))