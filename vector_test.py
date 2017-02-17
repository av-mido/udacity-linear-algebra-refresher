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