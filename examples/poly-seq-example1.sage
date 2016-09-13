## Falling factorial polynomials

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

n = var('n')
poly_seq_func = lambda k: expand(simplify(binomial(n, k) * factorial(k)))
pseq_data = map(poly_seq_func, range(1, 6))
guess_polynomial_sequence(pseq_data, n, index_offset = 1)

