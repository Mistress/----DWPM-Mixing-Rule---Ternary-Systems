#!/usr/bin/env python

from scipy import *
import scipy.optimize
from scipy.linalg import eig

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
        
        x1, x2 = x[0:2]

        if x1 + x2 > 1:
            return NaN

        else:
            
            Lambda_12 = self.c12/C[0]
            Lambda_21 = self.c21/C[1]
            Lambda_13 = self.c13/C[0]
            Lambda_31 = self.c31/C[2]
            Lambda_23 = self.c23/C[1]
            Lambda_32 = self.c32/C[2]

       
            delGmix = -(-x2-x1+1)*log(Lambda_32**self.s_3*x2-x2+Lambda_31**self.s_3*x1-x1+1)/self.s_3-x1*log(Lambda_12**self.s_1*x2+Lambda_13**self.s_1*(-x2-x1+1)+x1)/self.s_1-x2*log(x2+Lambda_23**self.s_2*(-x2-x1+1)+Lambda_21**self.s_2*x1)/self.s_2+x2*log(x2)+log(-x2-x1+1)*(-x2-x1+1)+x1*log(x1)
            
            return float(delGmix)
            
    def Gradient(self, x, C):
        
        x1, x2 = x[0:2]

        if x1 + x2 > 1:
            return NaN

        else:
            
            Lambda_12 = self.c12/C[0]
            Lambda_21 = self.c21/C[1]
            Lambda_13 = self.c13/C[0]
            Lambda_31 = self.c31/C[2]
            Lambda_23 = self.c23/C[1]
            Lambda_32 = self.c32/C[2]

            dGmix_dx1 = log(Lambda_32**self.s_3*x2-x2+Lambda_31**self.s_3*x1-x1+1)/self.s_3-log(Lambda_12**self.s_1*x2+Lambda_13**self.s_1*(-x2-x1+1)+x1)/self.s_1-(Lambda_31**self.s_3-1)*(-x2-x1+1)/(self.s_3*(Lambda_32**self.s_3*x2-x2+Lambda_31**self.s_3*x1-x1+1))-(1-Lambda_13**self.s_1)*x1/(self.s_1*(Lambda_12**self.s_1*x2+Lambda_13**self.s_1*(-x2-x1+1)+x1))-(Lambda_21**self.s_2-Lambda_23**self.s_2)*x2/(self.s_2*(x2+Lambda_23**self.s_2*(-x2-x1+1)+Lambda_21**self.s_2*x1))-log(-x2-x1+1)+log(x1)

            dGmix_dx2 = log(Lambda_32**self.s_3*x2-x2+Lambda_31**self.s_3*x1-x1+1)/self.s_3-log(x2+Lambda_23**self.s_2*(-x2-x1+1)+Lambda_21**self.s_2*x1)/self.s_2+log(x2)-(Lambda_32**self.s_3-1)*(-x2-x1+1)/(self.s_3*(Lambda_32**self.s_3*x2-x2+Lambda_31**self.s_3*x1-x1+1))-(Lambda_12**self.s_1-Lambda_13**self.s_1)*x1/(self.s_1*(Lambda_12**self.s_1*x2+Lambda_13**self.s_1*(-x2-x1+1)+x1))-(1-Lambda_23**self.s_2)*x2/(self.s_2*(x2+Lambda_23**self.s_2*(-x2-x1+1)+Lambda_21**self.s_2*x1))-log(-x2-x1+1)
            return array([dGmix_dx1, dGmix_dx2])

    def Hessian(self, x, C):
        
        x1, x2 = x[0:2]
        
        if x1 + x2 > 1:
            return NaN
        
        else:
            
            Lambda_12 = self.c12/C[0]
            Lambda_21 = self.c21/C[1]
            Lambda_13 = self.c13/C[0]
            Lambda_31 = self.c31/C[2]
            Lambda_23 = self.c23/C[1]
            Lambda_32 = self.c32/C[2]

            Hessian = empty((2, 2), float)

            Hessian[0, 0] = 2*(Lambda_31**s_3-1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1- x1+1))+(Lambda_31**s_3-1)**2*(-x2-x1+1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1)**2)-2*(1-Lambda_13**s_1)/(s_1*(Lambda_12**s_1*x2+Lambda_13**s_1*(-x2-x1+1)+x1))+(1-Lambda_13**s_1)**2*x1/(s_1*(Lambda_12**s_1*x2+Lambda_13**s_1*(-x2-x1+1)+x1)**2)+(Lambda_21**s_2-Lambda_23**s_2)**2*x2/(s_2*(x2+Lambda_23**s_2*(-x2-x1+1)+Lambda_21**s_2*x1)**2)+1/(-x2-x1+1)+1/x1

            Hessian[0, 1] = (Lambda_32**s_3-1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1))+(Lambda_31**s_3-1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1))+(Lambda_31**s_3-1)*(Lambda_32**s_3-1)*(-x2-x1+1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1)**2)-(Lambda_12**s_1-Lambda_13**s_1)/(s_1*(Lambda_12**s_1*x2+Lambda_13**s_1*(-x2-x1+1)+x1))+(1-Lambda_13**s_1)*(Lambda_12**s_1-Lambda_13**s_1)*x1/(s_1*(Lambda_12**s_1*x2+Lambda_13**s_1*(-x2-x1+1)+x1)**2)-(Lambda_21**s_2-Lambda_23**s_2)/(s_2*(x2+Lambda_23**s_2*(-x2-x1+1)+Lambda_21**s_2*x1))+(1-Lambda_23**s_2)*(Lambda_21**s_2-Lambda_23**s_2)*x2/(s_2*(x2+Lambda_23**s_2*(-x2-x1+1)+Lambda_21**s_2*x1)**2)+1/(-x2-x1+1)

            Hessian[1, 0] = (Lambda_32**s_3-1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1))+(Lambda_31**s_3-1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1))+(Lambda_31**s_3-1)*(Lambda_32**s_3-1)*(-x2-x1+1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1)**2)-(Lambda_12**s_1-Lambda_13**s_1)/(s_1*(Lambda_12**s_1*x2+Lambda_13**s_1*(-x2-x1+1)+x1))+(1-Lambda_13**s_1)*(Lambda_12**s_1-Lambda_13**s_1)*x1/(s_1*(Lambda_12**s_1*x2+Lambda_13**s_1*(-x2-x1+1)+x1)**2)-(Lambda_21**s_2-Lambda_23**s_2)/(s_2*(x2+Lambda_23**s_2*(-x2-x1+1)+Lambda_21**s_2*x1))+(1-Lambda_23**s_2)*(Lambda_21**s_2-Lambda_23**s_2)*x2/(s_2*(x2+Lambda_23**s_2*(-x2-x1+1)+Lambda_21**s_2*x1)**2)+1/(-x2-x1+1)

            Hessian[1, 1] = 2*(Lambda_32**s_3-1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1))+(Lambda_32**s_3-1)**2*(-x2-x1+1)/(s_3*(Lambda_32**s_3*x2-x2+Lambda_31**s_3*x1-x1+1)**2)+(Lambda_12**s_1-Lambda_13**s_1)**2*x1/(s_1*(Lambda_12**s_1*x2+Lambda_13**s_1*(-x2-x1+1)+x1)**2)-2*(1-Lambda_23**s_2)/(s_2*(x2+Lambda_23**s_2*(-x2-x1+1)+Lambda_21**s_2*x1))+(1-Lambda_23**s_2)**2*x2/(s_2*(x2+Lambda_23**s_2*(-x2-x1+1)+Lambda_21**s_2*x1)**2)+1/x2+1/(-x2-x1+1)

            HessEigens = eig(Hessian, left=False, right=False)

            if all(HessEigens>0):
                return 1
            else:
                return 0 
