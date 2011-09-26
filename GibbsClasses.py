#!/usr/bin/env python

from scipy import *
import scipy.optimize

class Model:
    R = 8.314
    name = 'Generic'
       
class DWPM(Model):
    def __init__(self, Params):
        self.name = 'DWPM'
        self.c12 = Params[0]
        self.c21 = Params[1]
        self.c13 = Params[2]
        self.c31 = Params[3]        
        self.c23 = Params[4]
        self.c32 = Params[5]
        
        self.s_1 = Params[6]
        self.s_2 = Params[7]
        self.s_3 = Params[8]
        
                         
    def deltaGmix(self, x, C):
        
        x1 = x[0]
        x2 = x[1]
        
        Lambda_12 = self.c12/C[0]
        Lambda_21 = self.c21/C[1]
        Lambda_13 = self.c13/C[0]
        Lambda_31 = self.c31/C[2]
        Lambda_23 = self.c23/C[1]
        Lambda_32 = self.c32/C[2]

        #if (x1==0.0 and x2==0.0) or (x1==0.0 and x2==1.0) or (x1==1.0 and x2==0.0) or (x1==1.0 and x2==1.0):
        #    delGmix = 0.0
        #else:
        delGmix = -(-x2-x1+1)*log(Lambda_32**self.s_3*x2-x2+Lambda_31**self.s_3*x1-x1+1)/self.s_3-x1*log(Lambda_12**self.s_1*x2+Lambda_13**self.s_1*(-x2-x1+1)+x1)/self.s_1-x2*log(x2+Lambda_23**self.s_2*(-x2-x1+1)+Lambda_21**self.s_2*x1)/self.s_2+x2*log(x2)+log(-x2-x1+1)*(-x2-x1+1)+x1*log(x1)
          
        return delGmix
            
  #  def FirstDerivative(self, x, T, c, M):
  #      
  #      dGmix_dx = GibbsFunctions.gibbsfunctions.d_dwpm_dx(x, self.s1, self.s2, self.c12, self.c21, c)
  #      
  #      return dGmix_dx
 
