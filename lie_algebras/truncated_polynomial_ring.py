from sympy import Symbol, symbols, expand
from sympy.core.numbers import mod_inverse
import numpy as np
import re

x_1, x_2, D_1, D_2 = symbols("x_1, x_2, D_1, D_2")


class TruncatedPolyRing:
    
    def __init__(self, p, m):
        self.p = p
        self.m = m
        self.x = [Symbol("x_" +  str(i + 1)) for i in range(m)]
    
    def get_summands(self, element):
        # Takes an element and returns a list of its summands
        summands = re.split("(?=[-+])", str(element))
        if summands[0] == "":
            summands = summands[1:]
        return list(map(eval, summands))
    
    def get_coef(self, element):
        # Gets coefficient of an element
        decomposed = str(element).split("*", 1)
        if len(decomposed) == 2 and isinstance(eval(decomposed[0]), int):
            return eval(decomposed[0]) % self.p, eval(decomposed[1])
        elif isinstance(element, int):
            return 1, element % self.p
        elif str(element)[0] == "-":
            return self.p - 1, -element
        else:
            return 1, element
    
    def get_exp(self, monomial):
        # Gets the exponent vector of a monomial
        monomials = re.split("\*(?!\*)", str(monomial))
        powers = [0] * self.m
        try:
            int(monomial)
            return powers
        except:
            pass
        i = 0
        while i < len(monomials):
            index = int(monomials[i][2]) - 1
            if monomials[i][-1] == "*":
                powers[index] = int(monomials[i + 1])
                i += 2
            else:
                powers[index] = 1
                i += 1
        return powers
    
    def coefs_mod_p(self, element):
        # Operates Modulo p on the coefficients of a general element in the ring
        return sum(list(map(lambda x: np.prod(self.get_coef(x)), self.get_summands(element))))
    
    def eliminate_exp_over_p(self, element):
        # Eliminates summands with exponents bigger than p - 1
        return sum(list(map(lambda x: x if max(self.get_exp(self.get_coef(x)[1])) < self.p else 0, self.get_summands(expand(element)))))
    
    def is_invertible(self, element):
        # Checks if an element is invertible in the ring
        if isinstance(self.get_summands(element)[-1], int):
            return True
        return False
    
    def get_zero_degree_component(self, summands):
        # Gets the first integer in a list of summands
        for summand in summands:
            if isinstance(summand, int):
                return summand
    
    def get_inverse(self, element):
        # Gets the inverse of an element in the ring using Geometric progression
        summands = self.get_summands(element)
        zero_deg_comp = self.get_zero_degree_component(summands)
        zero_deg_comp_inv = mod_inverse(zero_deg_comp, self.p)
        positive_deg_comp = sum(summands) - zero_deg_comp
        if not positive_deg_comp:
            return zero_deg_comp_inv
        inverse = zero_deg_comp_inv*sum((-zero_deg_comp_inv*positive_deg_comp)**i for i in range(self.p*len(summands[:-1])))
        return self.coefs_mod_p(self.eliminate_exp_over_p(inverse))
    
    def fraction_to_inverse(self, element):
        # Transforms an element of the form a/b to a*b^{-1}
        numerator, denominator = list(map(eval, re.split("/", str(element))))
        return numerator * self.get_inverse(denominator)
    
    def fix_inv_list(self, inv_list):
        # Converts all fractions in a list to inverses
        return list(map(lambda x: self.fraction_to_inverse(x) if "/" in str(x) else x, inv_list))
