from sympy import symbols, Matrix, diff, Poly, expand
import numpy as np
from itertools import islice, product
from jacobson_witt_algebra import JWitt

x_1, x_2, D_1, D_2 = symbols("x_1, x_2, D_1, D_2")


class Aut:
    # Class for calculations of automorphisms of the Jacobson-Witt algebra
   
    def __init__(self, ring):
        # Initiailzes the truncated polynomial ring, the algebra and the automorphisms
        self.ring = ring
        self.algebra = JWitt(ring.p, ring.m)
        self.truncated_auto_list = self.load_autos()
    
    def load_autos(self):
        # Loads the truncated_autos methods
        autos = []
#         autos.extend(self.truncated_autos())
#         autos.extend(self.truncated_autos2())
#         autos.extend(self.truncated_autos3())
#         autos.extend(self.truncated_autos4())
        autos.extend(self.truncated_autos5())
        return autos
    
    def truncated_autos(self):
        # Returns all automorphisms of the truncated polynomial ring in the form of x_1^i_1 + x_2^i_2 + ... whose J matrix is invertible
        images = []
        autos = []
        for exp in product(range(self.ring.p), repeat=self.ring.m):
            images.append(sum(list(map(lambda x, i: x**i if i != 0 else 0, self.ring.x, exp))))
        for mapping in product(images, repeat=self.ring.m):
            if self.is_auto(mapping):
                autos.append(mapping)
        return autos
    
    def truncated_autos2(self):
        # More automorphisms of the forms x_i + x^{(i_1,...,i_m)} whose J matrix is invertible
        products = []
        autos = []
        for exp in product(range(1, self.ring.p), repeat=self.ring.m):
            products.append(np.prod(list(map(lambda x, i: x**i if i != 0 else 0, self.ring.x, exp))))
        for i in islice(product(range(self.ring.m), repeat=self.ring.m), 1, 2**self.ring.m):
            for prod in products:
                mapping = tuple(np.multiply(i, prod) + self.ring.x)
                if self.is_auto(mapping):
                    autos.append(mapping)
        return autos
    
    def truncated_autos3(self):
        # Elementary automorphisms of the form x_1 goes to x_1 + f(x_2) when f is without free terms and the rest go to themselves; and exchanging x_1 with x_2
        autos = []
        for prod in product(range(self.ring.p), repeat=self.ring.p - 1):
            x_1_add_fx_2 = [x_1 + np.dot(prod, [x_2**i for i in range(1, self.ring.p)])] + self.ring.x[1:]
            x_2_add_fx_1 = [x_2 + np.dot(prod, [x_1**i for i in range(1, self.ring.p)])] + self.ring.x[:1]
            if self.is_auto(x_1_add_fx_2):
                autos.append(tuple(x_1_add_fx_2))
            if self.is_auto(x_2_add_fx_1):
                autos.append(tuple(x_2_add_fx_1))
        return autos
    
    def truncated_autos4(self):
        # Triangular automorphisms of the form x_1 goes to x_1 * f(x_2); rest go to themselves
        autos = []
        for prod in product(range(self.ring.p), repeat=self.ring.p):
            x_1_multiply_fx_2 = [x_1 * np.dot(prod, [x_2**i for i in range(self.ring.p)])] + self.ring.x[1:]
            if self.is_auto(x_1_multiply_fx_2):
                autos.append(tuple(x_1_multiply_fx_2))
        return autos
    
    def truncated_autos5(self):
        # Automorphisms of the form x_1 goes to x_1 + f(x_2) and x_2 goes to x_2 + g(x_1) when f, g are without free terms
        autos = []
        for prod1 in product(range(self.ring.p), repeat=self.ring.p - 1):
            x_1_add_fx_2 = x_1 + np.dot(prod1, [x_2**i for i in range(1, self.ring.p)])
            for prod2 in product(range(self.ring.p), repeat=self.ring.p - 1):
                if random.random() > 0.9:
                    x_2_add_fx_1 = x_2 + np.dot(prod2, [x_1**i for i in range(1, self.ring.p)])
                    mapping = [x_1_add_fx_2, x_2_add_fx_1]
                    if self.is_auto(mapping):
                        autos.append(tuple(mapping))
        return autos
    
    def get_J_matrix(self, mapping, show_matrix=False):
        # Find the J matrix of an automorphism
        J = Matrix()
        for i in self.ring.x:
            J = Matrix([J, list(map(diff, mapping, [i] * self.ring.m))])
        if show_matrix:
            print(J)
        return J
    
    def is_auto(self, mapping, show_matrix=False):
        # Checks if the J matrix of a truncated polynomial mapping is invertible
        J = self.get_J_matrix(mapping, show_matrix)
        if Poly(J.det(), self.ring.x).TC() != 0:
            return True
        return False
    
    def auto_image(self, element, auto):
        # Returns the image of an element under an automorphism
        J_inv = Matrix([self.ring.fix_inv_list(self.get_J_matrix(auto).inv().row(i)) for i in range(self.ring.m)])
        result = 0
        summands = self.ring.get_summands(element)
        for summand in summands:
            coef, exp, d = self.algebra.decompose(summand)
            auto_element = np.prod(list(map(lambda x, i: x**i, auto, exp)))
            d_index = int(str(d)[-1])
            result += expand(coef * auto_element * np.dot(list(J_inv.row(d_index - 1)), self.algebra.D))
        return self.ring.coefs_mod_p(self.algebra.alg_eliminate_exp_over_p(result))
    
    def all_auto_image(self, element, is_homogeneous=False):
        # Returns a list of images of an element under all automorphisms
        is_homogen = ""
        for index, auto in enumerate(self.truncated_auto_list):
            image = self.auto_image(element, auto)
            if is_homogeneous:
                if self.algebra.is_homogeneous(image):
                    is_homogen = "(Homogeneous)"
                else:
                    is_homogen = "(Non-homogeneous)"
            print(f'{index + 1}. {image} {is_homogen}')
