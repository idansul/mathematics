from sympy.functions.combinatorial.factorials import binomial
from sympy import Symbol, symbols
import numpy as np
from itertools import islice, product
import re
from truncated_polynomial_ring import TruncatedPolyRing

x_1, x_2, D_1, D_2 = symbols("x_1, x_2, D_1, D_2")


class JWitt(TruncatedPolyRing):
    # A class for calculations on the Jacobson-Witt Lie Algebras W(m;1) over a field of characteristic p
    
    def __init__(self, p, m):
        # Initializes the variables, derivations and basis
        super().__init__(p, m)
        self.D = [Symbol("D_" + str(i + 1)) for i in range(m)]
        self.basis = self.create_basis()
        self.dimension = m * p ** m
    
    def create_basis(self):
        # Creates the basis set of the algebra
        monoms = []
        monomials = []
        for x_i in self.x:
            monoms.append([x_i ** i for i in range(self.p)])
        for mon in product(*reversed(monoms)):
            monomials.append(np.prod(mon))
        return [monomial * D for D in self.D for monomial in monomials]
    
    def get_vector_presentation(self, element):
        # Gets the vector presentation of an element with respect to the basis
        element_vector = [0]*self.dimension
        for summand in self.get_summands(element):
            coefficient, basis_element = self.get_coef(summand)
            if basis_element:
                element_vector[self.basis.index(basis_element)] = coefficient
        return element_vector
    
    def vector_to_element(self, vector):
        # Transforms a vector presentation to its corresponding element
        return np.dot(vector, self.basis)
    
    def create_SC_table(self):
        # Creates structure constants table: if {e_i} is the basis, [e_i,e_j] = sum c_{ijk}*e_k; c_{ijk} are the structure constants
        sc_table = []
        for e_i in self.basis:
            c_i_row = []
            for e_j in self.basis:
                c_ij = self.lie(e_i, e_j)
                c_i_row.append(self.get_vector_presentation(c_ij))
            sc_table.append(c_i_row)
        self.sc_table = sc_table
    
    def extract_derivation(self, element):
        # Gets an element from the algebra and returns its polynomial and derivation components separately
        if element:
            d = eval(re.findall(f"D_[1-{self.m}]", str(element))[0])
            return element/d, d
        else:
            return 0, 0
    
    def decompose(self, element):
        # Decomposes an element into its coefficient, monomial exponent vector and derivation components
        coef, element = self.get_coef(element)
        monomial, d = self.extract_derivation(element)
        exp = self.get_exp(monomial)
        return coef, exp, d
    
    def alg_eliminate_exp_over_p(self, element):
        # Eliminates summands with exponents bigger than p - 1
        return sum(list(map(lambda x: x if max(self.decompose(x)[1]) < self.p else 0, self.get_summands(element))))
    
    def exp_to_monomial(self, exp):
        # Converts an exponent vector and a derivarion index into a corresponding element
        return np.prod(list(map(lambda i: self.x[i] ** exp[i], range(self.m))))

    def component_calc(self, a, b, epsilon, d):
        # Calculates a component in the Lie bracket
        d_index = int(str(epsilon)[-1])
        epsilon = [0 if i + 1 != d_index else 1 for i in range(self.m)]
        n = [a[i] + b[i] - epsilon[i] for i in range(self.m)]
        if any(i < 0 or i >= self.p for i in n):
            return 0
        binom = np.prod([binomial(n[i], a[i]) for i in range(self.m)]) % self.p
        return binom * self.exp_to_monomial(n) * d
    
    def lie_brackets(self, element1_exp, element1_d, element2_exp, element2_d):
        # Calculates the Lie brackets of 2 single elements
        return self.component_calc(element1_exp, element2_exp, element1_d, element2_d) \
                - self.component_calc(element2_exp, element1_exp, element2_d, element1_d)
        
    def lie(self, A, B):
        # Calculates the Lie brackets of 2 general element: [A,B]
        A_summands = self.get_summands(A)
        B_summands = self.get_summands(B)
        result = 0
        for element1 in A_summands:
            element1_coef, element1_exp, element1_d = self.decompose(element1)
            for element2 in B_summands:
                element2_coef, element2_exp, element2_d = self.decompose(element2)
                prod = self.lie_brackets(element1_exp, element1_d, element2_exp, element2_d) * element2_coef
                result += prod * element1_coef
        return self.coefs_mod_p(result)
    
    def search_bracket(self, elements, start=1, stop=None, verbose=False, verbose_result=False, removals=[None]):
        # Searches bracket presentations for a list of elements
        basis = [b for b in self.basis if b not in removals]
        M = [[i*bs for i in range(self.p)] for bs in basis]
        result = None
        results = {}
        for i in islice(product(*M), start, stop):
            left = sum(i)
            for j in islice(product(*M), 1, None):
                right = sum(j)
                result = self.lie(left, right)
                if verbose:
                    print(f"[{str(left)}, {str(right)}] =", result)
                if result in elements:
                    if verbose_result:
                        print(f"A bracket was found: [{left}, {right}] = {result}")
                    results[tuple([left, right])] = result
                    elements.remove(result)
                if not elements:
                    print(f"All of the brackets were found in the slice [start={start}, stop={stop}]")
                    return results
        print(f"The brackets for the following elements could not be found in the slice [start={start}, stop={stop}]: {elements}")
        if results:
            return results
    
    def get_all_elements_dict(self, start=1, stop=None, slices=None, left_removals=[None], right_removals=[None]):
        # Gets a dictionary of all unique elements that are presentable by a single bracket
        M_left = [[i*bs for i in range(self.p)] for bs in self.basis if bs not in left_removals]
        M_right = [[i*bs for i in range(self.p)] for bs in self.basis if bs not in right_removals]
        result = None
        results = {}
        for index, left in enumerate(islice(product(*M_left), start, stop)):
            if slices:
                if index + 1 not in slices:
                    continue
                else:
                    print(index + 1, len(results))
                if index + 1 > max(slices):
                    break
            left = sum(left)
            for right in islice(product(*M_right), 1, None):
                right = sum(right)
                result = self.lie(left, right)
                if result not in results.values():
                    results[tuple([left, right])] = result
                if len(results) == self.p ** self.dimension:
                    return results
        return results
    
    def epsilon_variations_calc(self, exponents):
        # Calculates [x^epsilon_i D_j, x^exp D_j] for all indices and the given exponents
        for exp in exponents:
            for d in product(self.D, repeat=2):
                for index in range(self.m):
                    epsilon = [0] * self.m
                    epsilon[index] = 1
                    left = self.exp_to_monomial(epsilon) * d[0]
                    right = self.exp_to_monomial(exp) * d[1]
                    print(f"[{str(left)}, {str(right)}] =", self.lie(left, right))
    
    def is_homogeneous(self, element):
        # Checks if an element is homogeneous or not (if all its summands have the same degree)
        summands = self.get_summands(element)
        degree = sum(self.decompose(summands[0])[1])
        for summand in summands[1:]:
            if sum(self.decompose(summand)[1]) != degree:
                return False
        return True
