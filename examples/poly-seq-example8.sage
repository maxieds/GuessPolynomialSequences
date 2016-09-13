## Exponential generating functions for the symmetric--indexed 
## binomial coefficients

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

n, z = var('n z')
poly_seq_func = lambda k: \
     sum(binomial(n + k, k) * (z ** n) / factorial(n), n, 0, infinity) * exp(-z)
user_guess_func = lambda n, k: 1 / factorial(k)
pseq_data = map(poly_seq_func, range(1, 6))
guess_polynomial_sequence(pseq_data, z, seq_factors = ["Binom2"], \
                          user_guess_func = user_guess_func, index_offset = 1)

