## Rising factorial polynomials 

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

n = var('n')
poly_seq_func = lambda m: expand(simplify(binomial(n + m, m) * factorial(m)))
pseq_data = map(poly_seq_func, range(1, 4))
guess_polynomial_sequence(pseq_data, n, index_offset = 1)

