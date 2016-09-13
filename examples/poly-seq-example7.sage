## Another exponential falling factorial polynomial example 
## with a user guess function 

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

n, z = var('n z')
poly_seq_func = lambda k: binomial(n, k)
user_guess_func = lambda n, k: 1 / factorial(k)
pseq_data = map(poly_seq_func, range(1, 6))
guess_polynomial_sequence(pseq_data, n, \
                          user_guess_func = user_guess_func, index_offset = 1)

