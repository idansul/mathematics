import sympy
from sympy.functions.combinatorial.factorials import binomial
from sympy import Symbol, symbols, simplify, combsimp
import numpy
from itertools import islice, product
import re

x_1, x_2, D_1, D_2 = symbols("x_1, x_2, D_1, D_2")


class JWitt:
    # A class for calculations on the Jacobson-Witt Lie Algebras W(m;1) over a field of characteristic p
    
    def __init__(self, p, m):
        # Initializes the variables, derivations and basis
        self.p = p
        self.m = m
        self.x = [Symbol("x_" +  str(i + 1)) for i in range(m)]
        self.D = [Symbol("D_" + str(i + 1)) for i in range(m)]
        self.epsilon = [Symbol("epsilon_" + str(i + 1)) for i in range(m)]
        monoms = []
        for x_i in self.x:
            monoms.append([x_i ** i for i in range(self.p)])
        self.monomials = []
        for mon in product(*reversed(monoms)):
            self.monomials.append(numpy.prod(mon))
        self.basis = [monomial * D for D in self.D for monomial in self.monomials]
        self.dimension = m * self.p ** m
    
    def mod_p(self, n):
        # Takes the Modulo p of an integer n
        if type(n) is sympy.Integer:
            return n % self.p
        return n
    
    def coefs_mod_p(self, element):
        # Operates Modulo p on the coefficients of a general element in the algebra
        decomposed = re.split("(?=[-+])", str(element))
        if decomposed[0] == "":
            decomposed = decomposed[1:]
        decomposed = list(map(eval, decomposed))
        return sum(list(map(lambda x: numpy.prod(self.get_coef(x)), decomposed)))
    
    def element_conversion(self, exp, d):
        # Converts an exponent vector and a derivarion index into a corresponding element
        x = numpy.prod(list(map(lambda i: self.x[i] ** exp[i], range(self.m))))
        return x * self.D[d - 1]

    def component_calc(self, a, b, epsilon, d):
        # Calculates a component in the Lie bracket
        epsilon = [0 if i + 1 != epsilon else 1 for i in range(self.m)]
        n = [a[i] + b[i] - epsilon[i] for i in range(self.m)]
        if any(i < 0 or i >= self.p for i in n):
            return 0
        binom = self.mod_p(numpy.prod([binomial(n[i], a[i]) for i in range(self.m)]))
        return binom * self.element_conversion(n, d)
    
    def lie_brackets(self, element1_exp, element1_d, element2_exp, element2_d):
        # Calculates the Lie brackets of 2 single elements
        return self.component_calc(element1_exp, element2_exp, element1_d, element2_d) \
                - self.component_calc(element2_exp, element1_exp, element2_d, element1_d)
    
    def get_coef(self, element):
        # Gets coefficient of an element
        decomposed = str(element).split("*", 1)
        if isinstance(eval(decomposed[0]), int):
            return eval(decomposed[0]) % self.p, eval(decomposed[1])
        elif decomposed[0][0] == "-":
            return self.p - 1, -element
        return 1, element
    
    def get_exp(self, monomial):
        # Gets the exponent vector of a monomial
        monomials = re.split("\*(?!\*)", str(monomial))
        i = 0
        powers = [0] * self.m
        while i < len(monomials):
            index = int(monomials[i][2]) - 1
            if monomials[i][-1] == "*":
                powers[index] = int(monomials[i + 1])
                i += 2
            else:
                powers[index] = 1
                i += 1
        return powers
    
    def decompose(self, element):
        # Decomposes an element into its coefficient, monomial exponent vector and derivation components
        coef, element = self.get_coef(element)
        decomposed = str(element).split("*", 1)
        d = int(decomposed[0][-1])
        if len(decomposed) == 2:
            exp = self.get_exp(eval(decomposed[1]))
        else:
            exp = [0] * self.m
        return coef, exp, d
        
    def lie(self, A, B, params=False):
        # Calculates the Lie brackets of 2 general element
        result = 0
        for element1 in A:
            element1_coef, element1_exp, element1_d = self.decompose(element1)
            for element2 in B:
                element2_coef, element2_exp, element2_d = self.decompose(element2)
                prod = self.lie_brackets(element1_exp, element1_d, element2_exp, element2_d) * element2_coef
                result += prod * element1_coef
        if result != 0 and params is False:
            return self.coefs_mod_p(result)
        elif params is True:
            return result
        return 0
    
    def lie_simp(self, A, B):
        # Returns a simplified presentation of self.lie
        return simplify(self.lie(A, B))
    
    def lie_combsimp(self, A, B):
        # Returns a combinatorially simplified presentation of self.lie
        return combsimp(self.lie(A, B))
    
    def search_bracket(self, elements, start=1, stop=None, verbose=False, verbose_result=False, removals=[None]):
        # Searches bracket presentations for a list of elements
        basis = [b for b in self.basis if b not in removals]
        M = [[i*bs for i in range(self.p)] for bs in basis]
        result = None
        results = {}
        for i in islice(product(*M), start, stop):
            i = [t for t in i if t != 0]
            for j in islice(product(*M), 1, None):
                j = [t for t in j if t != 0]
                result = self.lie(i, j)
                if verbose:
                    print(f"[{str(i)}, {str(j)}] =", result)
                if result in elements:
                    if verbose_result:
                        left_element = " +".join(re.split(",", str(i)[1:-1]))
                        right_element = " +".join(re.split(",", str(j)[1:-1]))
                        print(f"A bracket was found: [{left_element}, {right_element}] = {result}")
                    results[(tuple(i), tuple(j))] = result
                    elements.remove(result)
                if not elements:
                    print(f"All of the brackets have been found in the slice [start={start}, stop={stop}]")
                    return results
        print(f"The brackets for the following elements could not be found in the slice [start={start}, stop={stop}]: {elements}")
        if results:
            return results
    
    def get_all_elements_dict(self, start=1, stop=None, slices=None):
        # Gets a dictionary of all unique elements that are presentable by a single bracket
        M = [[i*bs for i in range(self.p)] for bs in self.basis]
        result = None
        results = {}
        for index, left in enumerate(islice(product(*M), start, stop)):
            if slices:
                if index + 1 not in slices:
                    continue
                else:
                    print(index + 1, len(results))
                if index + 1 > max(slices):
                    break
            left = [element for element in left if element != 0]
            print(left)
            for right in islice(product(*M), 1, None):
                right = [element for element in right if element != 0]
                result = self.lie(left, right)
                if result not in results.values():
                    results[(tuple(left), tuple(right))] = result
        return results
