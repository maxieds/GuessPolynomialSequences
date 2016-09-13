## Binomial squared difference operator identity

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

n, i, w = var('n i w')
poly_seq_func = lambda d: sum((binomial(d, i) ** 2) * factorial(d - i) * \
                              (n ** i), i, 0, d)
pseq_data = map(poly_seq_func, range(1, 8))
guess_polynomial_sequence(pseq_data, n, seq_factors = ["Binom2"], \
                          index_offset = 1)

