## SpecialSequences.py 
 # Implements several well-known triangular sequences in two indices 
 # and a couple of special functions still not stock in Sage
 # @author Maxie D. Schmidt
 # @since 2016.05.19
## 

import sys
import numpy as np

from sage.all import *
from sage.symbolic.function_factory import function_factory

def pochhammer(x, n):
     product_factor_func = lambda i: x + i
     prod_func = prod(map(product_factor_func, range(0, n)))
     return prod(factor_list)
## def

def harmonic_number(n, r = 1):
     k = var('k')
     return sum(k ** -r, k, 1, n)
## def

def identity_func(): return lambda n: n
def zero_func(): return lambda n: 0
def one_func(): return lambda n, k = 1: 1
def constant_func(constant): return lambda n: constant
def linear_func(a, b): return lambda n: a * n + b
def sequence_func(seq_data):
     return lambda n: seq_data[n] if n >= 0 and n < len(seq_data) else 0

class triangular_sequence(object):

     def __init__(self, a, b, g, ap, bp, gp, max_rows = 10):
          self.alpha = a 
          self.beta = b
          self.gamma = g
          self.alpha_prime = ap
          self.beta_prime = bp
          self.gamma_prime = gp
          self.stored_rows = max_rows
          self.rec_data = [[1]]
          self.generate_rows(1, max_rows)
     ## def 

     def generate_rows(self, row_min, row_max):
          [a, b, g, ap, bp, gp] = [self.alpha, self.beta, self.gamma, \
                                   self.alpha_prime, self.beta_prime, \
                                   self.gamma_prime]
          for n in range(row_min, row_max + 1): 
               row_data = []
               for k in range(0, n + 1): 
                    recnk = (a * n + b * k + g) * self.get_data(n - 1, k) + \
                            (ap * n + bp * k + gp) * self.get_data(n - 1, k - 1)
                    row_data += [recnk]
               ##
               self.rec_data += [row_data]
          ##
          self.stored_rows = row_max
     ## def 

     def get_data(self, n, k): 
          if n < 0 or k < 0 or k > n or n > self.stored_rows:
               return 0
          elif n > self.stored_rows:
               self.generate_rows(self.stored_rows + 1, n)
          ##
          return self.rec_data[n][k]
     ## def 

     def print_table_rows(self, row_min, row_max):
          #max_element = np.max(self.rec_data)
          for n in range(row_min, row_max + 1):
               for k in range(0, n + 1): 
                    sys.stdout.write("% 10d" % self.get_data(n, k))
               ## 
               print ""
          ##
     ## def 

     def print_table(self):
          self.print_table_rows(0, self.stored_rows)
     ## def

## class 
     
S1Triangle = triangular_sequence(1, 0, -1, 0, 0, 1)
S2Triangle = triangular_sequence(0, 1, 0, 0, 0, 1)
E1Triangle = triangular_sequence(0, 1, 1, 1, -1, 0)
E2Triangle = triangular_sequence(0, 1, 1, 2, -1, -1)
BinomTriangle = triangular_sequence(0, 0, 1, 0, 0, 1)

def S1(n, k): return S1Triangle.get_data(n, k)
def S2(n, k): return S2Triangle.get_data(n, k)
def E1(n, k): return E1Triangle.get_data(n, k)
def E2(n, k): return E2Triangle.get_data(n, k)
def Binom(n, k): return BinomTriangle.get_data(n, k)
def BinomSymmetric(n, k): return BinomTriangle.get_data(n + k, k)
def Binom2(n, k): return BinomTriangle.get_data(n, k) ** 2

def get_special_function_by_name(str_name):
     if str_name == "Identity":
          return lambda n, k: 1
     if str_name == "S1":
          return lambda n, k: S1(n, k)
     elif str_name == "S2":
          return lambda n, k: S2(n, k)
     elif str_name == "E1":
          return lambda n, m: E1(n, m)
     elif str_name == "E2":
          return lambda n, k: E2(n, k)
     elif str_name == "Binom":
          return lambda n, k: Binom(n, k)
     elif str_name == "BinomSym":
          return lambda n, k: Binom(n + k, k)
     elif str_name == "Binom2":
          return lambda n, k: Binom2(n, k)
     else: 
         return None
## def 

