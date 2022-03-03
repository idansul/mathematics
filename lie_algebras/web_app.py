import streamlit as st
from sympy import symbols, latex
import timeit
from jacobson_witt_algebra import JWitt
from truncated_polynomial_ring import TruncatedPolyRing
from automorphisms import Aut

x_1, x_2, D_1, D_2 = symbols('x_1, x_2, D_1, D_2')

start = timeit.default_timer()

def lie_brackets():
    st.latex(r'\text{Calculate Lie brackets $[A, B]$:}')
    c1, c2 = st.columns(2)
    with c1:
        A = st.text_input('A:')
    with c2:
        B = st.text_input('B:')
    if A and B:
        st.latex(A)
        st.latex(jwitt.lie(A, B))

def sidebar():
    if st.sidebar.button('Show dimension'):
        st.sidebar.latex(jwitt.dimension)
    if st.sidebar.button('Show basis'):
        st.sidebar.write(jwitt.basis)

def automorphism(p, m):
    st.header('Automorphisms')
    aut = Aut(TruncatedPolyRing(p, m))
    st.latex(r'\text{Check if a truncated polynomial mapping $\phi$ induces an algebra automorphism $\Phi$:}')
    indexes = st.columns(m)
    mapping = []
    for i, idx in enumerate(indexes):
        with idx:
            mapping.append(st.text_input(f'\phi(x_{i + 1}):'))
    if all(mapping):
        mapping = list(map(eval, mapping))
        if st.checkbox('Show J matrix'):
            st.latex(latex(aut.get_J_matrix(mapping)))
        if aut.is_auto(mapping):
            st.latex(r'\text{$\Phi$ is an Automorphism}')
            st.latex(r'\text{Compute an image $\Phi(-)$:}')
            element = st.text_input(r'Insert element:')
            if element:
                st.latex(aut.auto_image(element, mapping))
        else:
            st.latex(r'\text{$\Phi$ is not an Automorphism}')


st.set_page_config(page_title='Jacobson-Witt Algebras', page_icon=':heavy_plus_sign:')
st.title('Jacobson-Witt Algebras Computation')
st.latex(r'\text{Let $W(m)$ be the Jacobson-Witt algebra over a field $\mathbb{F}$ of characteristic $p>0$.}')

m = st.text_input('m:')
p = st.text_input('p:')

if p and m:
    m = int(m)
    p = int(p)
    jwitt = JWitt(p, m)
    lie_brackets()  
    sidebar()
    automorphism(p, m)
    
stop = timeit.default_timer()
st.text(f'Running time: {round(stop - start, 4)} sec')

st.markdown('[Github](https://github.com/idansul/mathematics/tree/main/lie_algebras)')
