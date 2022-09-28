from sympy import Symbol, Eq, solve
from sage.all import PolynomialRing, GF
import numpy as np


class EqSolver:
    # A class for computing and solving symbolic equations
    
    def __init__(self, equations):
        # Initializes the equations and gamma variables
        self.eq = equations
        self.gamma = tuple([Symbol(f'gamma_{i}') for i in range(len(equations))])
    
    def assign_var(self, vars_vals):
        # Assigns variables to the equations
        equations = self.eq
        for var, val in vars_vals.items():
            equations = list(map(lambda eq: eq.subs(var, val), equations))
        return equations
    
    def solve_eq(self, eqs, eq_num):
        # Solves and equation
        solution = solve(Eq(eqs[eq_num], self.gamma[eq_num]))[0]
        if solution:
            return solution
        return 0
    
    def solve_all_eqs(self, eqs):
        # Solves a list of equations
        return [self.solve_eq(eqs, i) for i in range(len(self.eq))]


class SageGroebner:
    # A class for transforming the single bracket presentation problem to Groebner bases: [a,b] = z; a, b are unknown
    
    def __init__(self, algebra, sc_table):
        # Initializes the ring and the SC table
        self.algebra = algebra
        self.a = [Symbol(f'a_{i}') for i in range(algebra.dimension)]
        self.b = [Symbol(f'b_{i}') for i in range(algebra.dimension)]
        self.ring = PolynomialRing(GF(algebra.p), algebra.dimension*2, self.a + self.b)
        self.sc_table = sc_table
        self.dim = algebra.dimension
        
    def groebner(self, polynomials):
        # Computes the Gröbner basis of the ideal that is generated from the set of polynomials with respect to the ring
        ideal = self.ring.ideal(polynomials)
        return ideal.groebner_basis(algorithm='macaulay2:f4')
    
    def create_equations(self, z):
        # Creates a list of the equations: sum_i sum_j a_i*b_j*c_{ijk} - z_k for each k with respect to an element z
        equations = []
        z_vector = self.algebra.get_vector_presentation(z)
        for k in range(self.dim):
            temp_eq = 0
            for i in range(self.dim):
                temp_eq += self.a[i]*np.dot(self.b, list(map(lambda j: self.sc_table[i][j][k], range(self.dim))))
            equations.append(temp_eq - z_vector[k])
        return equations
    
    def is_narrow(self, z):
        # Checks if an element is presentable as a single bracket by checking if its underlying Gröbner basis is nontrivial
        return 1 not in self.groebner(self.create_equations(z))
