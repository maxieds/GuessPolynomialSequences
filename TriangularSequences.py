## TriangularSequences.py 
 # Implements several well-known triangular sequences in two indices 
 # @author Maxie D. Schmidt
 # @since 2016.05.19
## 

import sys
import numpy as np

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
     ## def 

     def get_data(self, n, k): 
          if n < 0 or k < 0 or k > n:
               return 0
          elif n > self.stored_rows:
               self.generate_rows(self.stored_rows + 1, n)
          ##
          return self.rec_data[n][k]
     ## def 

     def print_table_rows(self, row_min, row_max):
          #max_element = np.max(self.rec_data)
          for n in range(row_min, row_max + 1):
               for k in range(1, n + 1): 
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

