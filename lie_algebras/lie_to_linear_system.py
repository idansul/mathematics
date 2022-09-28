from sage.all import GF, matrix, vector, gap
import numpy as np
from itertools import product


class LieToLinearSystem:
    # A class for transforming the single bracket problem in a Lie algebra to a linear system solution: [a,b] = z
    
    def __init__(self, p, sc_table):
        # Initializes the characteristic of the field of algebra and the structure constants table
        self.p = p
        self.sc_table = sc_table
        self.dim = len(sc_table)
    
    def equations_matrix(self, y):
        # Creates a list of the equations: sum_i sum_j a_i*y_j*c_{ijk} for each k with respect to an element z = sum_k z_k*e_k
        equations = []  
        for k in range(self.dim):
            equations.append(
                list(map(lambda i: np.dot(y, list(map(lambda j: self.sc_table[i][j][k], range(self.dim)))), range(self.dim))))
        return matrix(GF(self.p), equations)
    
    def is_element_narrow(self, z):
        # Checks if an element is representable as a single bracket by checking if its underlying linear system is solvable
        z_vector = vector(z)
        for y in product(reversed(range(self.p)), repeat=self.dim):
            eq_mat = self.equations_matrix(y)
            if eq_mat.augment(z_vector).rank() == eq_mat.rank():
                return True
        return False
    
    def is_algebra_narrow(self):
        # Checks if the algebra is narrow by checking the solvability of the linear systems of all of its elements
        for element in product(range(self.p), repeat=self.dim):
            if not self.is_element_narrow(element):
                return False
        return True
    
    def engel_equations_matrix(self, y):
        # Creates the linear system that corresponds to the Engel's equation [[x,y],y] where y is given and x is unknown
        equations = []
        for m in range(self.dim):
            equations.append(list(map(lambda i:
                         sum([y[j]*y[k]*self.sc_table[i][j][l]*self.sc_table[l][k][m]
                              for j, k, l in product(range(self.dim), repeat=3)]), range(self.dim))))
        return matrix(GF(self.p), equations)
    
    def is_element_engel(self, z):
        # Checks if an element is representable as [[x,y],y] by checking if its underlying (Engel) linear system is solvable
        z_vector = vector(z)
        for y in product(reversed(range(self.p)), repeat=self.dim):
            eq_mat = self.engel_equations_matrix(y)
            if eq_mat.augment(z_vector).rank() == eq_mat.rank():
#                 print('Engel element:', z_vector, 'for y:', y)
                return True
        return False
    
    def is_algebra_engel(self):
        # Checks if the algebra is Engel by checking the solvability of the (Engel) linear systems of all of its elements
        for element in product(range(self.p), repeat=self.dim):
            if not self.is_element_engel(element):
#                 print('Non-Engel element:', element)
                return False
        return True


class LowDimensionalLieAlgebras:
    # Class for low dimensional Lie algebras calculations over GF(2) using the FinLie package in GAP
    
    def __init__(self):
        # Initializes the FinLie package from GAP
        gap.eval('LoadPackage("finlie");')
    
    def convert_to_numeric(self, matrix):
        # Converts the elements of an adjoint matrix from FinLie to numeric values
        return eval(matrix.replace('0*Z(2)', '0').replace('Z(2)^0', '1'))
    
    def get_admatrices(self, dim, index=1):
        # Returns the adjoint matrices of the algebra
        gap.eval(f'L := LieAlgebraByLibrary(2, {dim}, {index})')
        gap.eval('B := Basis(L)')
        matrix = gap.eval('List(B, x -> AdjointMatrix(B, x))')
        return self.convert_to_numeric(matrix)
    
    def get_all_admatrices_up_to_dimension(self, dimension):
        # Returns all adjoint matrices up to a given dimension
        admatrices = {}
        admatrices_of_dim = []
        for dim in range(1, dimension + 1):
            for index in range(1, 20):
                if gap.eval(f'L := LieAlgebraByLibrary(2, {dim}, {index})') == 'fail':
                    break
                admatrices_of_dim.append(self.get_admatrices(dim, index))
            if admatrices_of_dim:
                admatrices[dim] = admatrices_of_dim
                admatrices_of_dim = []
        return admatrices
