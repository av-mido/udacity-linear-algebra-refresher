from plane import Plane
from vector import Vector
from linsys import LinearSystem
from decimal import Decimal

p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

s = LinearSystem([p0, p1, p2, p3])

def test_case_1():
    s.swap_rows(0,1)
    assert not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3)

def test_case_2():
    s.swap_rows(1,3)
    assert not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0)

def test_case_3():
    s.swap_rows(3,1)
    assert not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3)

# print('linear system is now: ', s)

def test_case_4():
    s.multiply_coefficient_and_row(1,0)
    assert not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3)
        # print('test case 4 failed')

def test_case_5():
    s.multiply_coefficient_and_row(-2,2)
    # print('After test case 5:')
    # print('   s[2]: ]', s[2])
    # print('   times -1 should be: ', Plane(normal_vector=Vector(['-2','-2','2']), constant_term='-6'))
    # print('   s[2].normal_vector[0]: ', s[2].normal_vector[0])
    # print('   s[2].basepoint: ', s[2].basepoint)
    assert not (s[0] == p1 and
            s[1] == p0 and
            s[2] == Plane(normal_vector=Vector(['-2','-2','2']), constant_term='-6') and
            s[2].normal_vector[0] == -2 and s[2].normal_vector[1] == -2 and s[2].normal_vector[2] == 2 and
            s[2].constant_term == -6 and
            s[3] == p3)

s.multiply_coefficient_and_row(10,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 6 failed')

s.add_multiple_times_row_to_row(0,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 7 failed')

# print('linear system before case 8:', s)
s.add_multiple_times_row_to_row(1,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 8 failed')

s.add_multiple_times_row_to_row(-1,1,0)
if not (s[0] == Plane(normal_vector=Vector(['-10','-10','-10']), constant_term='-10') and
        s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 9 failed')


######################################## Triangular for tests ##################################


p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2):
    print('test case 1 failed')

p1 = Plane(normal_vector=Vector(['0','0','3']), constant_term='3')
p2 = Plane(normal_vector=Vector(['0','2','2']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
s = LinearSystem([p1,p2,p3])
t = s.compute_triangular_form()
# print(t)
if not (t[0] == p3 and
        t[1] == p2 and
        t[2] == p1):
    print('test case 1.1 failed')

# p1 = Plane(normal_vector=Vector(['0','0','3']), constant_term='3')
# p2 = Plane(normal_vector=Vector(['0','2','2']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
# s = LinearSystem([p1,p2,p3])
# t = s.compute_triangular_form()
# if not (t[0] == p3 and
#         t[1] == p2 and
#         t[2] == p1):
#     print('test case 1.2 failed')

# p1 = Plane(normal_vector=Vector(['0','0','3']), constant_term='3')
# p2 = Plane(normal_vector=Vector(['0','2','2']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['0','0','1']), constant_term='1')
# s = LinearSystem([p1,p2,p3])
# print('1.3')
# t = s.compute_triangular_form()
# print(t)
# if not (t[0] == p2 and
#         t[1] == p1 and
#         t[2] == p3):
#     print('test case 1.3 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == Plane(constant_term='1')):
    print('test case 2 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2 and
        t[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        t[3] == Plane()):
    print('test case 3 failed')

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
t = s.compute_triangular_form()
if not (t[0] == Plane(normal_vector=Vector(['1','-1','1']), constant_term='2') and
        t[1] == Plane(normal_vector=Vector(['0','1','1']), constant_term='1') and
        t[2] == Plane(normal_vector=Vector(['0','0','-9']), constant_term='-2')):
    print('test case 4 failed')



############################# Quiz: Coding RREF #################################
p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='-1') and
        r[1] == p2):
    print('rref test case 1 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term='1')):
    print('rref test case 2 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='0') and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        r[3] == Plane()):
    print('rref test case 3 failed')

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term=Decimal('23')/Decimal('9')) and
        r[1] == Plane(normal_vector=Vector(['0','1','0']), constant_term=Decimal('7')/Decimal('9')) and
        r[2] == Plane(normal_vector=Vector(['0','0','1']), constant_term=Decimal('2')/Decimal('9'))):
    print('rref test case 4 failed')


#################################### Quiz: Coding GE Solution ################
p1 = Plane(normal_vector=Vector(['5.862','1.178','-10.366']), constant_term='-8.15')
p2 = Plane(normal_vector=Vector(['-2.931','-0.589','5.183']), constant_term='-4.075')
s = LinearSystem([p1, p2])
try:
    ge = s.compute_gaussian_elimination()
    print('failed GE test 1')
except Exception as e:
    if str(e) == LinearSystem.NO_SOLUTIONS_MSG:
        pass # print('passed GE test 1')
    else:
        print('failed GE test 1')

p1 = Plane(normal_vector=Vector(['8.631','5.112','-1.816']), constant_term='-5.113')
p2 = Plane(normal_vector=Vector(['4.315','11.132','-5.27']), constant_term='-6.775')
p3 = Plane(normal_vector=Vector(['-2.158','3.01','-1.727']), constant_term='-0.831')
s = LinearSystem([p1, p2, p3])
try:
    ge = s.compute_gaussian_elimination()
    print('failed GE test 2')
except Exception as e:
    if str(e) == LinearSystem.INF_SOLUTIONS_MSG:
        pass # print('passed GE test 2')
    else:
        print('failed GE test 2')

p1 = Plane(Vector([5.262, 2.739, -9.878]), -3.441)
p2 = Plane(Vector([5.111, 6.358, 7.638]), -2.152)
p3 = Plane(Vector([2.016, -9.924, -1.367]), -9.278)
p4 = Plane(Vector([2.167, -13.543, -18.883]), -10.567)
s = LinearSystem([p1, p2, p3, p4])
try:
    ge = s.compute_gaussian_elimination()
    if not ge == Vector([Decimal('-1.17720187578995858313947665146'), Decimal('0.707150558138740933006474968216'), Decimal('-0.0826635849022828890650647196936')]):
        print('failed GE test 3')
except Exception as e:
    print('failed GE test 3')


##################################### Quiz: Coding Parametrization ###########################
# Hooray all of these pass.
# See https://classroom.udacity.com/courses/ud953/lessons/4624329808/concepts/48972686550923#
# TODO: turn these into unit tests.

# Fixed equation from github -
p1 = Plane(Vector([0.786, 0.786, 0.588]), -0.714)
p2 = Plane(Vector([-0.131, -0.131, 0.244]), 0.319)
s = LinearSystem([p1, p2])
param_version = s.compute_gaussian_elimination(True)
print("\nparam test 1:")
print(param_version)

p1 = Plane(Vector([8.631, 5.112, -1.816]), -5.113)
p2 = Plane(Vector([4.315, 11.132, -5.27]), -6.775)
p3 = Plane(Vector([-2.158, 3.01, -1.727]), -0.831)
s = LinearSystem([p1, p2, p3])
param_version = s.compute_gaussian_elimination(True)
print("\nparam test 2:")
print(param_version)

p1 = Plane(Vector([0.935, 1.76, -9.365]), -9.955)
p2 = Plane(Vector([0.187, 0.352, -1.873]), -1.991)
p3 = Plane(Vector([0.374, 0.704, -3.746]), -3.982)
p4 = Plane(Vector([-0.561, -1.056, 5.619]), 5.973)
s = LinearSystem([p1, p2, p3, p4])
param_version = s.compute_gaussian_elimination(True)
print("\nparam test 3:")
print(param_version)

######################