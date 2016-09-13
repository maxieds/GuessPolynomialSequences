## A mock example of a polynomial sequence defined by a 
## finite summation formula over the binomial and Stirling number triangles 

from GuessPolynomialSequenceFunction import *
n, x, t = var('n x t')
S1_func = lambda n, k: S1(n, k)
S2_func = lambda n, k: S2(n, k)
def poly_seq_func(n): 
     poly = 0
     for i in range(0, n + 1):
          poly += (-5 ** i) * S2_func(n, i) * S1_func(n, i) * \
                  binomial(n, i) * (x ** i)
     ##
     return poly
poly_seq_func_hold = lambda n: poly_seq_func(n)
pseq_data = map(poly_seq_func_hold, range(1, 6))
guess_polynomial_sequence(pseq_data, x, seq_factors = ["S2", "S1", "Binom"], \
                          index_offset_params = [(0, 1)], 
                          index_offsets = [0, 1], index_offset = 1)

