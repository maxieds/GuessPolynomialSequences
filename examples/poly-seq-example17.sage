## Finite summation formula for the non-exponential Bell polynomials

from GuessPolynomialSequenceFunction import guess_polynomial_sequence

def series_coefficient_zpow(f, fvar, ncoeff): \
     return f.taylor(fvar, 0, ncoeff) - f.taylor(fvar, 0, ncoeff - 1)
def series_coefficient(f, fvar, ncoeff): 
     sc = series_coefficient_zpow(f, fvar, ncoeff)
     return sc.subs_expr(fvar == 1) 

n, x, t = var('n x t')
spoly_ogf = exp( (exp(t) - 1) * x)

poly_seq_func = lambda n: series_coefficient(spoly_ogf, t, n) * factorial(n)
pseq_data = map(poly_seq_func, range(1, 6))
guess_polynomial_sequence(pseq_data, x, seq_factors = ["S2"])

