import sympy as sp
from sympy.functions.combinatorial.factorials import binomial
from sympy import Symbol, symbols, simplify, combsimp
from itertools import islice, product
import re

x, D, D_i, D_j, a, b, epsilon_i, epsilon_j = symbols("x, D, D_i, D_j, a, b, epsilon_i, epsilon_j")

class Lie:
    # A class for calculations on the Jacobson-Witt Lie Algebras W(m;n) over a field of characteristic p
    
    def __init__(self, p):
        # Initializes the variables, derivations and basis
        self.p = int(p)
        self.x = Symbol("x")
        self.D = Symbol("D")
        self.D_i = Symbol("D_i")
        self.D_j = Symbol("D_j")
        self.epsilon_i = Symbol("epsilon_i")
        self.epsilon_j = Symbol("epsilon_j")
        self.basis = [self.x ** i * self.D for i in range(self.p)]
    
    def mod_p(self, n):
        # Operates Modulo p on an integer
        if type(n) is int:
            return n % self.p
        return n
    
    def coefs_mod_p(self, element):
        # Operates Modulo p on the coefficients of a general element in the algebra
        s = re.split("([+-])", str(element))
        s = list(map(str.strip, s))
        if s[0] == "":
            s = s[1:]
        i = 0
        mod_s = 0
        if s[0] != "+" and s[0] != "-":
            s_coef, s_exp, s_d = self.decompose(s[0])
            mod_s += (s_coef % self.p) * x ** s_exp * s_d
            i += 1
        while i < len(s):
            s_coef, s_exp, s_d = self.decompose(s[i + 1])
            if s[i] == "+":
                mod_s += (s_coef % self.p) * x ** s_exp * s_d
            else:
                mod_s += (-s_coef % self.p) * x ** s_exp * s_d
            i += 2
        return mod_s
    
    def binomial_check(self, a, b, epsilon, k):
        # Checks if the binomial coefficient equals 0 or not
        if a == 0 and k == b or b == 0 and k == a:
            return 0
        return 1

    def component_calc(self, a, b, epsilon, k, d):
        # Calculates a component in the Lie bracket
        n = a + b - epsilon
        return self.binomial_check(a, b, epsilon, k) * binomial(n, k) * self.x ** self.mod_p(n) * d
    
    def eps(self, d):
        # Returns the epsilon that corresponds to a given derivation
        if d == self.D_i:
            return self.epsilon_i
        elif d == self.D_j:
            return self.epsilon_j
        return 1
    
    def lie_brackets(self, a, di, b, dj):
        # Calculates the Lie brackets of 2 single elements
        if isinstance(a, int) and isinstance(b, int):
            if a + b > self.p - 1:
                return 0
        return self.component_calc(a, b, self.eps(di), a, dj) - self.component_calc(a, b, self.eps(dj), b, di)
    
    def get_coef(self, element):
        # Gets coefficient of an element
        r = str(element).split("*", 1)
        if isinstance(eval(r[0]), int):
            return eval(r[0]), eval(r[1])
        return 1, element
    
    def get_exp(self, monomial):
        # Gets an exponent of a monomial
        if len(str(monomial)) != 1:
            return eval(str(monomial)[3:])
        return 1
    
    def decompose(self, element):
        # Decomposes an element into its coefficient, monomial exponent and derivation components
        coef, element = self.get_coef(element)
        r = str(element).split("*", 1)
        d = eval(r[0])
        if len(r) == 2:
            exp = self.get_exp(eval(r[1]))
        else:
            exp = 0
        return coef, exp, d
        
    def lie(self, A, B, params=False):
        # Calculates the Lie brackets of 2 general element
        result = 0
        for s in A:
            s_coef, s_exp, s_d = self.decompose(s)
            for t in B:
                t_coef, t_exp, t_d = self.decompose(t)
                prod = self.lie_brackets(s_exp, s_d, t_exp, t_d) * t_coef
                result += prod * s_coef
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
    
    def search_bracket(self, element, start=1, stop=None, verbose=False):
        # Searches a single bracket presentation of an element
        M = [[i*bs for i in range(self.p)] for bs in self.basis]
        for i in islice(product(M[0], M[1], M[2], M[3], M[4]), start, stop):
            i = [t for t in i if t != 0]
            if len(i) == self.p:
                continue
            for j in islice(product(M[0], M[1], M[2], M[3], M[4]), 1, None):
                j = [t for t in j if t != 0]
                if type(sum(i)/sum(j) + 6) not in {sp.Integer, sp.Float, sp.Rational}:
                    result = self.lie(i, j)
                    if verbose:
                        print("[" + str(i) + ",", str(j) + "] =", result)
                    if result == element:
                        print("The bracket was found!")
                        print("[" + str(i) + ",", str(j) + "] =", result)
                        break
            if result == element:
                break
        if result != element:
            print("The bracket for the element {} was not found in the slice [start={}, stop={}].".format(element, start, stop))


p = input("Enter p: ")
l = Lie(p)
element = x * D + x ** 4 * D
l.search_bracket(element, start=1500, stop=1700)

# A = [D_i, D_j]
# B = [x**(a + epsilon_i + epsilon_j)*D_i, x**(b + epsilon_i + epsilon_j)*D_j]
# print(l.lie(A, B))

# A = [D_i]
# B = [x**a*D_j]
# print(l.lie(A, B))

# cc = 'y'
# while cc == 'y':
#     print(l.lie(eval(input("Enter first: ")), eval(input("Enter second: ")), True))
#     cc = input("Calculate more? (y/n): ")