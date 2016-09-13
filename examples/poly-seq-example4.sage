## First-order Eulerian numbers in the OGFs of the polylogarithm 
## functions, Li_{-m}(z), for natural numbers m >= 1

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

n, z = var('n z')
poly_seq_func = lambda m: expand(\
     factor(sum((n ** m) * (z ** n), n, 0, infinity)) *\
     ((1 - z) ** (m + 1))\
).subs(z = n)

pseq_data = map(poly_seq_func, range(1, 6))
guess_polynomial_sequence(pseq_data, n, seq_factors = ["E1"], index_offset = 1)

