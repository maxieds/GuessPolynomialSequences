## GuessPolynomialSequenceFunction.py 
 # 
 # @author Maxie D. Schmidt
 # @since 2016.05.19
##

from SpecialSequences import *
from GuessSequenceFunctions import *

SPFN_NAME, SPFN_UINDEX, SPFN_LINDEX = 0, 1, 2
UIDXA, UIDXB, LIDXA, LIDXB = 1, 2, 3, 4
FMD_FACTORFUNC, FMD_EXPFUNC, FMD_SPFNDATA, FMD_REMSEQ = 0, 1, 2, 3

def extract_formula_match_data(fmatch, data_param):
     return fmatch[data_param]
## def 

def get_index_func(a, b, index_param):
     if index_param == -1:
          return lambda j, i: a * (j - i) + b
     elif index_param == 1:
          return lambda j, i: a * i + b
     elif index_param == 0:
          return lambda j, i: a * j + b
     else:
          return None
## def 

def preprocess_polynomial_sequence(seq, pvar, user_guess_func, index_offset):
     
     if user_guess_func == None:
          user_guess_func = one_func()
     ##

     # first, determine any sequences of polynomial multiples:
     poly_factor_exps, poly_exps = [], []
     for (pidx, poly) in enumerate(seq):
          if pidx == 0 and is_integer(poly):
               continue
          [min_exp, max_exp] = get_polynomial_exponents(poly, pvar)
          poly_factor_exps += [min_exp]
          poly_exps += [max_exp - min_exp]
     ## 

     poly_factor_func = guess_sequence_function(poly_factor_exps)
     poly_factor_func = zero_func() if poly_factor_func == None \
                        else poly_factor_func
     poly_exp_func = guess_sequence_function(poly_exps)
     poly_exp_func = identity_func() if poly_exp_func == None\
                     else poly_exp_func

     remove_factors_func = lambda pidx: \
          seq[pidx] if is_integer(seq[pidx]) else \
          expand(seq[pidx] * (pvar ** (-1 * poly_factor_func(pidx)))) 
     tseq = map(remove_factors_func, range(0, len(seq)))

     # next, preprocess the coefficients by removing factors of the 
     # user guess function:
     for (pidx, tpoly) in enumerate(tseq):
          if is_integer(tpoly): 
               continue
          tpoly2 = tpoly.polynomial(QQ)
          tpoly_coeffs = tpoly2.coefficients(sparse = False)
          coeff_transform_func = lambda cfidx: \
               take_ratio(tpoly_coeffs[cfidx], \
                          user_guess_func(cfidx, pidx + index_offset))
          tpoly_coeffs = map(coeff_transform_func, range(0, len(tpoly_coeffs)))
          tseq[pidx] = to_polynomial(tpoly_coeffs, pvar)
     ##

     return tseq, poly_factor_func, poly_exp_func

## def 

def get_remaining_sequence(poly, n, px, seq_factor_func, ff, integerQ = True):
     
     if is_integer(poly):
          remterm = take_ratio(poly, seq_factor_func(n, 0 + ff(0)))
          if not is_integer(remterm) or\
             (seq_factor_func(n, 0 + ff(0)) == 0 and poly != 0):
               return None
          ##
          return [remterm]
     ##
     
     poly2 = poly.polynomial(QQ)
     poly_coeffs = poly2.coefficients(sparse = False)
     get_remainder_func = lambda cidx: \
                   Rational(take_ratio(poly_coeffs[cidx], \
                            seq_factor_func(n, cidx + ff(cidx))))\
                   if not is_zero(seq_factor_func(n, cidx + ff(cidx))) or\
                      is_zero(poly_coeffs[cidx]) else Integer(0)
     rseq = map(get_remainder_func, range(0, len(poly_coeffs)))
     #print "   >>> (poly) ", poly
     #print "   >>> (polycfs) ", poly_coeffs
     #print "   >>> ", map(lambda i: seq_factor_func(n, i + ff(i)), range(0, len(poly_coeffs)))
     for (tidx, term) in enumerate(rseq):
          #print "   >>> ", term, ", ", type(term), ", ", term.is_integer()
          if integerQ and not is_integer(term):
               return None
          ##
     ##
     return rseq

## def

def verify_formula_matches_raw(seq, px, exp_func, spfactor_func, ff, index_offset):

     cur_rseq, revQ = [], None
     for is_reversed in [True, False]:
       for (nidx, poly) in enumerate(seq):
          
          nidx += index_offset
          rseq = get_remaining_sequence(poly, nidx, px, spfactor_func, ff)
          if rseq == None and is_reversed:
               break
          elif rseq == None:
               #print "   > Returning here II ... "
               return False, None, None, None
          elif cur_rseq == []:
               cur_rseq = rseq
               continue
          ##
          trunc_rseq = rseq[0:exp_func(nidx - 1) + 1]
          rev_trunc_rseq = reverse(rseq)[0:exp_func(nidx - 1) + 1]
          
          #print "   > (revQ) ", is_reversed
          #print "   > (currseq) ", cur_rseq
          #print "   > (rseq) ", rseq
          #print "   > (trseq) ", trunc_rseq
          #print "   > (revtrseq) ", rev_trunc_rseq

          #if is_reversed == None and \
          #     check_sequence_diff(trunc_rseq, sequence_func(cur_rseq)):
          #     is_reversed = False
          #elif is_reversed == None:
          #     is_reversed = True
          if trunc_rseq != [] and len(cur_rseq) > 1 and not is_reversed and \
               (not check_sequence_diff(trunc_rseq, lambda n: cur_rseq[n]) or\
               is_reversed and not check_sequence_diff(\
               rev_trunc_rseq, lambda n: reverse(trunc_rseq)[n])):
               if is_reversed: 
                    break
               else: 
                    #print "   > Returning here I ... "
                    return False, None, None, None
          ##
          cur_rseq = rseq

       ## 
       if is_reversed:
            cur_rseq = reverse(cur_rseq)
       ##
       revQ = is_reversed
       if len(cur_rseq) == seq[len(seq) - 1].degree(px) + 1: 
            break
     ##
     #print revQ, cur_rseq, cur_rseq[0:seq[0].degree(px) + 1], constant_array(seq[0].degree(px) + 1, 0)

     poly0_degree = seq[0].degree(px)
     poly0_min_exponent = seq[0].polynomial(QQ).exponents()[0]
     if cur_rseq == constant_array(len(cur_rseq), 0) or\
        cur_rseq[0:poly0_degree + 1] == constant_array(poly0_degree + 1, 0) or\
        revQ and reverse(cur_rseq)[0:poly0_degree + 1] == \
        constant_array(poly0_degree + 1, 0) or\
        poly0_min_exponent == 0 and cur_rseq[0] == 0 or\
        poly0_min_exponent > 0 and cur_rseq[0:poly0_min_exponent] == \
        constant_array(poly0_min_exponent, 0):
          return False, spfactor_func, cur_rseq, is_reversed
     else:
          return True, spfactor_func, cur_rseq, revQ

## def 

def get_sequence_formula_func(factor_func, exp_func, spfn_func, rem_seq, revQ):
     i = var('i')
     rseq_func = sequence_func(rem_seq if not revQ else reverse(rem_seq))
     seq_func = lambda n, x: \
                sum(spfn_func(n, i) * rseq_func(i) * \
                    (x ** (factor_func(i) + exp_func(i))), i, 0, exp_func(n))
     return seq_func
## def

def get_spfactor_func_latex(spfn_data): 
     sdata_latex = ""
     for (spfn_name, uidx_func, lidx_func) in spfn_data: 
          sdata_latex += "%s(%s, %s) " %\
                         (spfn_name, eval_on_operands(uidx_func)(n, i), \
                          eval_on_operands(lidx_func)(n, i))
     ##
     return sdata_latex
## def 

def print_sequence_match(seq, factor_func, exp_func, spfn_data, \
                         user_guess_func, rem_seq, is_reversed):
     
     x, n, i = var('x n i')
     factor_func_print = eval_on_operands(factor_func)
     exp_func_print = eval_on_operands(exp_func)
     
     print ""
     print "[ ============================================================= ]"
     print ""
     print seq
     print ""
     print "Polynomial Factor Function:          x^{%s}" % factor_func_print(n)
     print "Upper Sum Index:                    ", exp_func_print(n)
     
     print "User Guess Func:                    ", \
           None if user_guess_func == None else user_guess_func(n, i)

     
     spseq_factor_str = get_spfactor_func_latex(spfn_data)
     print "Special Sequence Factor:             %s" % spseq_factor_str
     
     print "Remaining Sequence Terms:           ", \
           rem_seq if not is_reversed else reverse(rem_seq)
     
     rseq_func_index = i if not is_reversed else exp_func_print(n) - identity_func()(i)
     rseq_func = guess_sequence_function(rem_seq)
     guess_rseq_func = eval_on_operands(rseq_func)
     print "Remaining Sequence Formula (Guess): ", \
           guess_rseq_func(n).subs(n = rseq_func_index) if rseq_func != None\
                                                     else None
     
     #rseq_str = latex(guess_rseq_func(n).subs(n = rseq_func_index)) \
     #           if rseq_func != None \
     #           else "RemSeq(%s)" % rseq_func_index
     rseq_str = "RemSeq(%s)" % rseq_func_index
     sdata_latex = get_spfactor_func_latex(spfn_data)
     latex_sum= "\\sum_{i=0}^{%s} %s%s %s \\times x^{i+%s}" %\
                (exp_func_print(n), sdata_latex, rseq_str, \
                 "" if user_guess_func == None else user_guess_func(n, i), \
                 factor_func_print(n))
     print "LaTeX Formula for the Sequence:      $%s$" % latex_sum 
     print ""
     print "[ ============================================================= ]"
     print ""

## def

def get_sequence_factor_func(fidx, spfn_name, ipdata): 
     uidx_func = get_index_func(ipdata[UIDXA + 5 * fidx], \
                                ipdata[UIDXB + 5 * fidx], \
                                ipdata[5 * fidx][0])
     lidx_func = get_index_func(ipdata[LIDXA + 5 * fidx], \
                                ipdata[LIDXB + 5 * fidx], \
                                ipdata[5 * fidx][1])
     spfactor_func = lambda n, i: \
       get_special_function_by_name(spfn_name)(uidx_func(n, i), lidx_func(n, i))
     return uidx_func, lidx_func, spfactor_func
## def 

def guess_polynomial_sequence(seq, px, seq_factors = ["S1"], \
                              index_offset_params = [(0, 1), (0, -1)], \
                              user_guess_func = None, \
                              index_multiples = [1], \
                              index_offsets = [-1, 0, 1],\
                              num_spfunc_rows = 12, index_offset = 1, \
                              print_summary = True):

     n, k, z, x = var('n k z x')
     seq, factor_func, exp_func = \
          preprocess_polynomial_sequence(seq, px, user_guess_func, index_offset)
     #print seq
     
     # update the exp_func formula according to index_offset:
     exp_poly = eval_on_operands(exp_func)(n)
     exp_poly = exp_poly.polynomial(QQ)
     [exp_poly] = transform_polynomial_variable([exp_poly], n, index_offset, n)
     exp_func = lambda x: exp_poly.subs({n: x})

     formula_matches = []

     # a straightforward lazy implementation (to refactor and make more like 
     # the original Mathematica package code later):
     [I1, I2, I3] = map(lambda I: FiniteEnumeratedSet(I), \
                        [index_offset_params, index_multiples, index_offsets])
     
     cpindex_sets = [I1, I2, I3, I2, I3]
     cpsets = []
     for seq_factor in seq_factors:
          cpsets += cpindex_sets
     ##
     index_params = cartesian_product(cpsets)

     for (idx, ipdata) in enumerate(index_params.list()):
          
          #if idx > 0: break
          uidx_funcs, lidx_funcs, spfactor_funcs = [], [], []
          for (fidx, spfn_name) in enumerate(seq_factors):
               uidx_func, lidx_func, spfactor_func = \
                          get_sequence_factor_func(fidx, spfn_name, ipdata)
               uidx_funcs += [uidx_func]
               lidx_funcs += [lidx_func]
               spfactor_funcs += [spfactor_func]
          ## 

          seq_factor_func = one_func()
          if len(seq_factors) == 1:
               seq_factor_func = lambda n, i: spfactor_funcs[0](n, i)
          elif len(seq_factors) == 2:
               seq_factor_func = lambda n, i: \
                          spfactor_funcs[0](n, i) * spfactor_funcs[1](n, i)
          ##
          elif len(seq_factors) == 3:
               seq_factor_func = lambda n, i: \
                          spfactor_funcs[0](n, i) * spfactor_funcs[1](n, i) *\
                          spfactor_funcs[2](n, i)
          ##
          elif len(seq_factors) == 4:
               seq_factor_func = lambda n, i: \
                          spfactor_funcs[0](n, i) * spfactor_funcs[1](n, i) *\
                          spfactor_funcs[2](n, i) * spfactor_funcs[3](n, i)
          else:
               print "We only support 1 <= m <= 4 sequence factors ... "
          ##

          spfn_data_map_func = lambda i: \
                              (seq_factors[i], uidx_funcs[i], lidx_funcs[i])
          spfn_data = map(spfn_data_map_func, range(0, len(seq_factors)))

          is_valid, spfactor_func, rem_seq, is_reversed = \
               verify_formula_matches_raw(seq, px, exp_func, seq_factor_func, \
                                          factor_func, index_offset = index_offset)
          
          if not is_valid:
               continue
          ##
          formula_matches += [ (factor_func, exp_func, spfn_data, \
                                seq_factor_func, [rem_seq, is_reversed]) ]
     ##
     if formula_matches == []:
          return []
     ##

     formula_funcs = []
     for fmd in formula_matches:
          (ffunc, efunc, sdata, sfunc, [rseq, revQ]) = fmd
          if print_summary:
               print_sequence_match(seq, ffunc, efunc, sdata, \
                                    user_guess_func, rseq, revQ)
          ##
          formula_funcs += [(get_sequence_formula_func(ffunc, efunc, sfunc, \
                             rseq, revQ), sdata)] ## TODO
     ##
     return formula_funcs

## def


