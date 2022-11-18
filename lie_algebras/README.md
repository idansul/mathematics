## The Jacobson-Witt Lie Algebras and their Automorphisms
- The definition of the Jacobson-Witt Lie Algebras can be found in this paper: https://arxiv.org/pdf/math/0508373.pdf, page 47.
- The motivation of the scripts is to make Lie brackets calculations easier; to be able to search single brackets presentation of any element in the algebra; and to make automorphism computations on the algebra.
- `witt_algebra.py` is for the Witt Algebra W(1;1). `jacobson_witt_algebra.py` is for the general Jacobson-Witt Algebra W(m;1), m > 1.
- `truncated_polynomial_ring.py` is for the corresponding truncated polynomial ring.
- `automorphisms.py` is for automorphisms computations.
- `functions.py` contains some useful search functions.
- Libraries to install: SymPy, NumPy, Sage.
## Lie brackets to Linear Systems
- `lie_to_linear_system.py` contains tools for transforming Lie brackets to linear systems to be checked if they are narrow or have the Engel property. In addition, it contains a tool for getting the adjoint matrices of some low-dimensional Lie algebras over F_2 from the FinLie package in GAP, see: http://www.iaa.tu-bs.de/beick/so.html.
- `misc_solvers.py` contains some more tools to solve linear systems either with a more manual fashion or using Gröbner bases (the latter requires Sage).
- Most of the tools in this section require the package SageMath 9.3 which is used as a Python notebook interface, see: https://www.sagemath.org/download.html. When the computer algebra system GAP is used it is loaded from the Sage package.
### Web application
https://lie-algebras.streamlit.app
