## A Legendre polynomial identity involving squared powers of the 
## binomial coefficients

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

n, z = var('n z')
poly_seq_func = lambda m: \
     expand( factor( legendre_P(m, (1+z)/(1-z)) * ((1-z) ** m) ) )

pseq_data = map(poly_seq_func, range(1, 6))
guess_polynomial_sequence(pseq_data, z, seq_factors = ["Binom2"])

