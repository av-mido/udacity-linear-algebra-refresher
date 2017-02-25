from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        self[row] = Plane(self[row].normal_vector * coefficient, self[row].constant_term * coefficient)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        plane_times_coef = Plane(self[row_to_add].normal_vector * coefficient, self[row_to_add].constant_term * coefficient)
        new_normal = self[row_to_be_added_to].normal_vector + plane_times_coef.normal_vector
        new_constant = self[row_to_be_added_to].constant_term + plane_times_coef.constant_term
        self[row_to_be_added_to] = Plane(new_normal, new_constant)

    def compute_triangular_form(self):
        system = deepcopy(self)
        for cur_row_idx in range(0, len(system.planes)):
            self.swap_if_necessary(cur_row_idx, system)

            if system.planes[cur_row_idx].normal_vector.is_zero_vector():
                return system

            # Do row subtractions
            first_nonzero_col = Plane.first_nonzero_index(system.planes[cur_row_idx].normal_vector)
            system.clear_below_rows(cur_row_idx, first_nonzero_col)
        return system

    def clear_below_rows(self, cur_row_idx, first_nonzero_col):
        for below_row in range(cur_row_idx + 1, len(self.planes)):
            a = self.planes[cur_row_idx].normal_vector[first_nonzero_col]
            b = self.planes[below_row].normal_vector[first_nonzero_col]
            coef = -b / a
            self.add_multiple_times_row_to_row(coef, cur_row_idx, below_row)

    def clear_above_rows(self, cur_row_idx, idx_of_first_nonzero_term):
        for above_row in range(cur_row_idx-1, -1, -1):
            # 'a' should be zero here
            a = self.planes[cur_row_idx].normal_vector[idx_of_first_nonzero_term]
            b = self.planes[above_row].normal_vector[idx_of_first_nonzero_term]
            coef = - b / a
            self.add_multiple_times_row_to_row(coef, cur_row_idx, above_row)

    @staticmethod
    def swap_if_necessary(cur_row_idx, system):
        # Do a swap if necessary...
        for cur_col_idx in range(0, system.dimension):
            cur_row_col_val = system.planes[cur_row_idx].normal_vector[cur_col_idx]
            if not MyDecimal(cur_row_col_val).is_near_zero():
                return
            else:
                for potential_swap_idx in range(cur_row_idx + 1, len(system.planes)):
                    potential_row_col_val = system.planes[potential_swap_idx].normal_vector[cur_col_idx]
                    if not MyDecimal(potential_row_col_val).is_near_zero():
                        system.swap_rows(cur_row_idx, potential_swap_idx)
                        return

    def compute_rref(self):
        tf = self.compute_triangular_form()
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()
        for row_idx in range(len(tf.planes)-1, -1, -1):
            idx_of_first_nonzero_term = pivot_indices[row_idx]
            if idx_of_first_nonzero_term < 0:
                continue
            first_nonzero_term = tf.planes[row_idx].normal_vector[idx_of_first_nonzero_term]
            if not (0.9999 < first_nonzero_term < 1.0001):
                tf.multiply_coefficient_and_row(Decimal(1)/first_nonzero_term, row_idx)
            tf.clear_above_rows(row_idx, idx_of_first_nonzero_term)
        return tf

    def compute_gaussian_elimination(self):
        rref = self.compute_rref()
        rref.raise_exception_if_contradictory_equation()
        rref.raise_exception_if_too_few_pivots()
        solution_coords = [rref.planes[i].constant_term for i in range(rref.dimension)]
        return Vector(solution_coords)

        #
        # pivot_indices = rref.indices_of_first_nonzero_terms_in_each_row()
        # xyz_values = [None] * rref.dimension
        # j = rref.dimension - 1
        # for row_idx in range(len(rref.planes)-1, -1, -1):
        #     idx_of_first_nonzero_term = pivot_indices[row_idx]
        #     if rref[row_idx].normal_vector[idx_of_first_nonzero_term] == -1 and not rref[row_idx].constant_term.is_near_zero():
        #         raise Exception("No Solution")
        #     if idx_of_first_nonzero_term == j:

    def raise_exception_if_contradictory_equation(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        for row_idx in range(0, len(self.planes)):
            idx_of_first_nonzero_term = pivot_indices[row_idx]
            if idx_of_first_nonzero_term == -1 and not self.planes[row_idx].constant_term.is_near_zero():
                raise Exception(self.NO_SOLUTIONS_MSG)

    def raise_exception_if_too_few_pivots(self):
        # col = 0
        # for row_idx in range(0, self.planes):
        #     if col >= self.dimension:
        #         return
        #     if self.planes[row_idx].normal_vector[col].is_near_zero():
        #         raise Exception("Infinite solutions")
        #     col += 1
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        if num_pivots < self.dimension:
            raise Exception(self.INF_SOLUTIONS_MSG)


    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret



class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


# p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
# p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
# p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
# p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
#
# s = LinearSystem([p0,p1,p2,p3])
#
# print('Indices of first nonzero terms: ', s.indices_of_first_nonzero_terms_in_each_row())
# print('{},  {},  {},  {}'.format(s[0], s[1], s[2], s[3]))
# print(len(s))
# print(s)
#
# s[0] = p1
# print(s)
#
# print(MyDecimal('1e-9').is_near_zero())
# print(MyDecimal('1e-11').is_near_zero())
