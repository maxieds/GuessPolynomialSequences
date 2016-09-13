## GuessSequenceFunctions.py 
 #
 # @author Maxie D. Schmidt
 # @since 2016.05.19
##

from sage.all import *
from sage.symbolic.function_factory import function_factory
from sage.symbolic.function_factory import eval_on_operands
from sage.symbolic.expression import Expression as SageExpr

import sympy
import numpy as np
from SpecialSequences import pochhammer

def take_ratio(p, q): 
     return p / q if q != 0 else 0
## def

def constant_array(alen, fill_value):
     a = np.empty(alen)
     a.fill(fill_value)
     return map(lambda e: int(e), list(a))
##

def is_integer(poly):
     return isinstance(poly, int) or \
            type(poly) == sage.rings.integer.Integer or\
            type(poly) == sage.rings.rational.Rational and poly.is_integer() or\
            type(poly) == sage.symbolic.expression.Expression and\
            poly.is_integer()
## def

def is_zero(poly):
     return isinstance(poly, int) and poly == 0 or\
            not isinstance(poly, int) and poly.is_zero()
## def

def to_polynomial(coeffs, pvar):
     poly = 0
     for i in range(0, len(coeffs)):
          poly += coeffs[i] * (pvar ** i)
     ##
     return poly
## def 

def reverse(lst):
     return list(reversed(lst))
## def

## Note: pvar not used yet ... 
def get_polynomial_exponents(poly, pvar):
     
     if is_integer(poly):
          return [0, 0]
     ##
     tpoly = poly.polynomial(QQ)
     poly_coeffs = tpoly.coefficients(sparse = False)
     for eidx in range(0, len(poly_coeffs)):
          if poly_coeffs[eidx] != 0:
               return [eidx, len(poly_coeffs) - 1]
          ##
     ##
     return None

## def 

def guess_sequence_function(input_seq, user_guess_func = None, \
                            use_sigma = False): 

     if len(input_seq) < 2:
          return None
     ##
     
     # apply the user guess function: 
     gseq = input_seq
     if user_guess_func != None:
          guess_map_func = lambda se, idx: take_ratio(se, user_guess_func(idx)) 
          gseq = map(guess_map_func, input_seq)
     ## 
     gseq_func = lambda gidx: 0 if gidx < 0 or gidx >= len(gseq) \
                              else gseq[gidx]

     # define variables for the specific forms we are looking for:
     a, n, r, c, a0, a1, a2 = var('a n r c a0 a1 a2')

     # check if the sequence is a constant power, a factorial / Pochhammer 
     # symbol, or some product of both:
     ratio_eqns = []
     for n in range(1, len(gseq)): 
          ratio = take_ratio(gseq_func(n), gseq_func(n - 1)) 
          ratio_eqns += [ratio == c * (a * (n - 1) + 1)]
     ##
     if ratio_eqns == constant_array(len(ratio_eqns), 1):
          return one_func()
     ##
     
     solve_vars = solve(ratio_eqns, c, a)
     if solve_vars != []: 
          
          solve_vars = solve_vars[0]
          c = c.subs(solve_vars[0]) 
          a = a.subs(solve_vars[1])
          constant_term = gseq_func(0)
          const_pow = c
          if a == 0:
               const_pow = c
               seq_func = lambda n: gseq_func(0) * (const_pow ** n)
               return seq_func
          #elif a == 1:
          #     seq_func = lambda n: gseq_func(0) * (const_pow ** n) * factorial(n)
          #     return seq_func
          #else:
          #     const_pow = c * a # from the Pochhammer symbol product rep.
          #     pochhammer_seq_func = lambda n: pochhammer(1 / a, n)
          #     seq_func = lambda n: gseq_func(0) * \
          #                          (const_pow ** n if const_pow != 0 else 0) * \
          #                          pochhammer_seq_func(n)
          #     return seq_func

          ##
     ## 

     # next, see if the sequence satisfies a recurrence of the form 
     # f[n] == f[n-1] + a[0] + a[1]*n + a[2]*n^2, i.e., the sequence 
     # corresponds to some polynomial function of small degree:
     poly_eqns = []
     for n in range(1, len(gseq)):
          poly_eqns += [gseq_func(n) == gseq_func(n - 1) + a0 + a1 * n]
     ## 
     solve_vars = solve(poly_eqns, a0, a1)
     if solve_vars != []:
          solve_vars = solve_vars[0]
          a0, a1, i = a0.subs(solve_vars[0]), a1.subs(solve_vars[1]), var('i')
          const_term = gseq_func(0)
          seq_func = lambda n: const_term + \
                     sum(a0 + a1 * i, i, 1, n)
          return seq_func
     ##

     return None

## def 

def check_sequence_diff(seq, seq_func):
     for n in range(0, len(seq)):
          if seq[n] != seq_func(n):
               return False
          ##
     ##
     return True
## def

def guess_linear_index_func(input_seq, quick = False):
     
     nrange = [0, 1]
     if not quick: 
          nrange = range(0, len(input_seq))
     ## 

     a, b = var('a b')
     linear_func_eqns = []
     for n in nrange:
          linear_func_eqns += [input_seq[n] == a * n + b]
     ##
     solve_vars = solve(linear_func_eqns, a, b)
     if solve_vars == []:
          return False, None, None
     else:
          return True, solve_vars[0][0].rhs(), solve_vars[0][1].rhs()

## def

def transform_polynomial_variable(poly_seq, pvar, pvar_offset, next_pvar): 
     tpoly_seq = map(lambda p: sympy.sympify(p), poly_seq)
     tpoly_map_func = lambda tp: \
               expand(tp.subs(pvar, next_pvar - pvar_offset))
     tpoly_seq = map(tpoly_map_func, tpoly_seq)
     return tpoly_seq
## def 

def transform_ffact_polynomial_basis(poly, pvar, ffact_pow_var):
     m = var('m')
     tpoly = expand(sympy.sympify(poly))
     ffact_poly_func = lambda p, w: sum(S2(p, m) * (w ** m), m, 0, p)
     [min_exponent, max_exponent] = get_polynomial_exponents(poly)
     for ppow in reverse(range(min_exponent, max_exponent + 1)):
          tpoly = tpoly.subs(pvar ** ppow, ffact_poly_func(ppow, ffact_pow_var))
     ## 
     return expand(tpoly)
## def 

def transform_rfact_polynomial_basis(poly, pvar, rfact_pow_var):
     m = var('m')
     tpoly = expand(sympy.sympify(poly))
     rfact_poly_func = lambda p, w: \
                       sum(S2(p, m) * (-1 ** p-m) * (w ** m), m, 0, p)
     [min_exponent, max_exponent] = get_polynomial_exponents(poly)
     for ppow in list(reversed(min_exponent, max_exponent + 1)):
          tpoly = tpoly.subs(pvar ** ppow, rfact_poly_func(ppow, rfact_pow_var))
     ## 
     return expand(tpoly)
## def 

