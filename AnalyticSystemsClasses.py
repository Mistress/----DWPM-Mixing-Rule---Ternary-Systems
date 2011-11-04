#!/usr/bin/env python

from scipy import *
import numpy

class SystemDWPM:
    
    def __init__(self, Params, Cell_s, C):
        self.s_1 = Cell_s[0]
        self.s_2 = Cell_s[1]
        self.s_3 = Cell_s[2]
        
        self.Lambda_12 = Params[0]/C[0]
        self.Lambda_21 = Params[1]/C[1]
    
        self.Lambda_23 = Params[4]/C[1]
        self.Lambda_32 = Params[5]/C[2]

        self.Lambda_13 = Params[2]/C[0]
        self.Lambda_31 = Params[3]/C[2]

    def NewN(self, x_alpha, x_Beta, m, step):
        
        direction = +1.
        tieline3 = numpy.append(x_Beta - x_alpha, 0)
        zaxis = numpy.array([0., 0., 1.])

        UnitDirection = numpy.cross(tieline3, zaxis)[:2]
        UnitDirection /= numpy.linalg.norm(UnitDirection)

        return m + direction*step*UnitDirection

    def SystemEquations(self, Predicted, n):

        x1_a, x2_a, x1_b, x2_b, m1, n1, c1 = Predicted
        xm_1, xm_2 = n
        
        Eq1 =  -(-x2_a-x1_a+1)*log(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-\
        x1_a+1)/self.s_3-x1_a*log(self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1*(-x2_a-\
        x1_a+1)+x1_a)/self.s_1-x2_a*log(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+\
        self.Lambda_21**self.s_2*x1_a)/self.s_2+x2_a*log(x2_a)-n1*x2_a+log(-x2_a-x1_a+1)\
        *(-x2_a-x1_a+1)+x1_a*log(x1_a)-m1*x1_a-c1

        Eq2 = -(-x2_b-x1_b+1)*log(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-\
        x1_b+1)/self.s_3-x1_b*log(self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-\
        x1_b+1)+x1_b)/self.s_1-x2_b*log(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+\
        self.Lambda_21**self.s_2*x1_b)/self.s_2+x2_b*log(x2_b)-n1*x2_b+log(-x2_b-x1_b+1)\
        *(-x2_b-x1_b+1)+x1_b*log(x1_b)-m1*x1_b-c1

        Eq3 = log(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1)/self.s_3-log(\
        self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a)/self.s_1-(\
        self.Lambda_31**self.s_3-1)*(-x2_a-x1_a+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+\
        self.Lambda_31**self.s_3*x1_a-x1_a+1))-(1-self.Lambda_13**self.s_1)*x1_a/(self.s_1*(self.Lambda_12**self.s_1\
        *x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a))-(self.Lambda_21**self.s_2\
        -self.Lambda_23**self.s_2)*x2_a/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)\
        +self.Lambda_21**self.s_2*x1_a))-log(-x2_a-x1_a+1)+log(x1_a)-m1

        Eq4 = log(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1)/self.s_3-log(\
        self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b)/self.s_1-(self.Lambda_31**self.s_3\
        -1)*(-x2_b-x1_b+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3\
        *x1_b-x1_b+1))-(1-self.Lambda_13**self.s_1)*x1_b/(self.s_1*(self.Lambda_12**self.s_1\
        *x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b))-(self.Lambda_21**self.s_2\
        -self.Lambda_23**self.s_2)*x2_b/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)\
        +self.Lambda_21**self.s_2*x1_b))-log(-x2_b-x1_b+1)+log(x1_b)-m1
        
        Eq5 = log(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1)/self.s_3-log(\
        x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+self.Lambda_21**self.s_2*x1_a)/self.s_2+log(\
        x2_a)-(self.Lambda_32**self.s_3-1)*(-x2_a-x1_a+1)/(self.s_3*(self.Lambda_32**self.s_3*\
        x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1))-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1\
        )*x1_a/(self.s_1*(self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1\
        )+x1_a))-(1-self.Lambda_23**self.s_2)*x2_a/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(\
        -x2_a-x1_a+1)+self.Lambda_21**self.s_2*x1_a))-log(-x2_a-x1_a+1)-n1

        Eq6 = log(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1)/self.s_3-log(\
        x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+self.Lambda_21**self.s_2*x1_b)/self.s_2+log(\
        x2_b)-(self.Lambda_32**self.s_3-1)*(-x2_b-x1_b+1)/(self.s_3*(self.Lambda_32**self.s_3*\
        x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1))-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1\
        )*x1_b/(self.s_1*(self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1\
        )+x1_b))-(1-self.Lambda_23**self.s_2)*x2_b/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(\
        -x2_b-x1_b+1)+self.Lambda_21**self.s_2*x1_b))-log(-x2_b-x1_b+1)-n1

        #Eq7 = (xm_1-x1_a)/(x1_b-x1_a)-(xm_2-x2_a)/(x2_b-x2_a)
        Eq7 = -n1*xm_2+(n1*x2_b-n1*x2_a+m1*x1_b-m1*x1_a)*(xm_1-x1_a)/(x1_b-x1_a)-m1*xm_1+n1*x2_a+m1*x1_a

        return real(array([Eq1, Eq2, Eq3, Eq4, Eq5, Eq6, Eq7]))

    def SystemEquationsJac(self, Predicted, n):

        x1_a, x2_a, x1_b, x2_b, m1, n1, c1 = Predicted
        xm_1, xm_2 = n

        Jac1 = [log(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1)/self.s_3-log(\
        self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a)/self.s_1-(\
        self.Lambda_31**self.s_3-1)*(-x2_a-x1_a+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+\
        self.Lambda_31**self.s_3*x1_a-x1_a+1))-(1-self.Lambda_13**self.s_1)*x1_a/(self.s_1*(self.Lambda_12**self.s_1\
        *x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a))-(self.Lambda_21**self.s_2\
        -self.Lambda_23**self.s_2)*x2_a/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a\
        +1)+self.Lambda_21**self.s_2*x1_a))-log(-x2_a-x1_a+1)+log(x1_a)-m1,log(self.Lambda_32**self.s_3\
        *x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1)/self.s_3-log(x2_a\
        +self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+self.Lambda_21**self.s_2*x1_a)/self.s_2+log(x2_a\
        )-(self.Lambda_32**self.s_3-1)*(-x2_a-x1_a+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-\
        x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1))-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1\
        )*x1_a/(self.s_1*(self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)\
        +x1_a))-(1-self.Lambda_23**self.s_2)*x2_a/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a\
        -x1_a+1)+self.Lambda_21**self.s_2*x1_a))-log(-x2_a-x1_a+1)-n1,0,0,-x1_a,\
        -x2_a,-1]

        Jac2 = [0,0,log(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1)/self.s_3-\
        log(self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b)/self.s_1\
        -(self.Lambda_31**self.s_3-1)*(-x2_b-x1_b+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b\
        +self.Lambda_31**self.s_3*x1_b-x1_b+1))-(1-self.Lambda_13**self.s_1)*x1_b/(self.s_1*(\
        self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b))-(self.Lambda_21**self.s_2\
        -self.Lambda_23**self.s_2)*x2_b/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-\
        x1_b+1)+self.Lambda_21**self.s_2*x1_b))-log(-x2_b-x1_b+1)+log(x1_b)-m1,\
        log(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1)/self.s_3-\
        log(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+self.Lambda_21**self.s_2*x1_b)/self.s_2+\
        log(x2_b)-(self.Lambda_32**self.s_3-1)*(-x2_b-x1_b+1)/(self.s_3*(self.Lambda_32**self.s_3*\
        x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1))-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1\
        )*x1_b/(self.s_1*(self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b\
        +1)+x1_b))-(1-self.Lambda_23**self.s_2)*x2_b/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(\
        -x2_b-x1_b+1)+self.Lambda_21**self.s_2*x1_b))-log(-x2_b-x1_b+1)-n1,-x1_b,\
        -x2_b,-1]

        Jac3 = [2*(self.Lambda_31**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3\
        *x1_a-x1_a+1))+(self.Lambda_31**self.s_3-1)**2*(-x2_a-x1_a+1)/(self.s_3*(self.Lambda_32**self.s_3\
        *x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1)**2)-2*(1-self.Lambda_13**self.s_1\
        )/(self.s_1*(self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a\
        +1)+x1_a))+(1-self.Lambda_13**self.s_1)**2*x1_a/(self.s_1*(self.Lambda_12**self.s_1*x2_a\
        +self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a)**2)+(self.Lambda_21**self.s_2-\
        self.Lambda_23**self.s_2)**2*x2_a/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+\
        self.Lambda_21**self.s_2*x1_a)**2)+1/(-x2_a-x1_a+1)+1/x1_a,(self.Lambda_32**self.s_3-1)\
        /(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1))+(\
        self.Lambda_31**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*\
        x1_a-x1_a+1))+(self.Lambda_31**self.s_3-1)*(self.Lambda_32**self.s_3-1)*(-x2_a-x1_a\
        +1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1)**2\
        )-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)/(self.s_1*(self.Lambda_12**self.s_1*x2_a+\
        self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a))+(1-self.Lambda_13**self.s_1)*(\
        self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)*x1_a/(self.s_1*(self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1\
        *(-x2_a-x1_a+1)+x1_a)**2)-(self.Lambda_21**self.s_2-self.Lambda_23**self.s_2)\
        /(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+self.Lambda_21**self.s_2*x1_a))\
        +(1-self.Lambda_23**self.s_2)*(self.Lambda_21**self.s_2-self.Lambda_23**self.s_2)*x2_a/(self.s_2*(\
        x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+self.Lambda_21**self.s_2*x1_a)**2)+1/(\
        -x2_a-x1_a+1),0,0,-1,0,0]

        Jac4 = [0,0,2*(self.Lambda_31**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3\
        *x1_b-x1_b+1))+(self.Lambda_31**self.s_3-1)**2*(-x2_b-x1_b+1)/(self.s_3*(\
        self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1)**2)-2*(1-\
        self.Lambda_13**self.s_1)/(self.s_1*(self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-\
        x1_b+1)+x1_b))+(1-self.Lambda_13**self.s_1)**2*x1_b/(self.s_1*(self.Lambda_12**self.s_1*\
        x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b)**2)+(self.Lambda_21**self.s_2-\
        self.Lambda_23**self.s_2)**2*x2_b/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+\
        self.Lambda_21**self.s_2*x1_b)**2)+1/(-x2_b-x1_b+1)+1/x1_b,(self.Lambda_32**self.s_3\
        -1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1)\
        )+(self.Lambda_31**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3\
        *x1_b-x1_b+1))+(self.Lambda_31**self.s_3-1)*(self.Lambda_32**self.s_3-1)*(-x2_b-\
        x1_b+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b\
        +1)**2)-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)/(self.s_1*(self.Lambda_12**self.s_1*\
        x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b))+(1-self.Lambda_13**self.s_1)*(self.Lambda_12**self.s_1\
        -self.Lambda_13**self.s_1)*x1_b/(self.s_1*(self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1\
        *(-x2_b-x1_b+1)+x1_b)**2)-(self.Lambda_21**self.s_2-self.Lambda_23**self.s_2\
        )/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+self.Lambda_21**self.s_2*x1_b\
        ))+(1-self.Lambda_23**self.s_2)*(self.Lambda_21**self.s_2-self.Lambda_23**self.s_2)*x2_b/(\
        self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+self.Lambda_21**self.s_2*x1_b)**2)\
        +1/(-x2_b-x1_b+1),-1,0,0]
        
        Jac5 = [(self.Lambda_32**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3*\
        x1_a-x1_a+1))+(self.Lambda_31**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a\
        +self.Lambda_31**self.s_3*x1_a-x1_a+1))+(self.Lambda_31**self.s_3-1)*(self.Lambda_32**self.s_3\
        -1)*(-x2_a-x1_a+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3\
        *x1_a-x1_a+1)**2)-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)/(self.s_1*(self.Lambda_12**self.s_1\
        *x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a))+(1-self.Lambda_13**self.s_1\
        )*(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)*x1_a/(self.s_1*(self.Lambda_12**self.s_1\
        *x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a)**2)-(self.Lambda_21**self.s_2\
        -self.Lambda_23**self.s_2)/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+self.Lambda_21**self.s_2\
        *x1_a))+(1-self.Lambda_23**self.s_2)*(self.Lambda_21**self.s_2-self.Lambda_23**self.s_2\
        )*x2_a/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)+self.Lambda_21**self.s_2\
        *x1_a)**2)+1/(-x2_a-x1_a+1),2*(self.Lambda_32**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3\
        *x2_a-x2_a+self.Lambda_31**self.s_3*x1_a-x1_a+1))+(self.Lambda_32**self.s_3\
        -1)**2*(-x2_a-x1_a+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_a-x2_a+self.Lambda_31**self.s_3\
        *x1_a-x1_a+1)**2)+(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)**2*x1_a/(self.s_1\
        *(self.Lambda_12**self.s_1*x2_a+self.Lambda_13**self.s_1*(-x2_a-x1_a+1)+x1_a)**2\
        )-2*(1-self.Lambda_23**self.s_2)/(self.s_2*(x2_a+self.Lambda_23**self.s_2*(-x2_a-x1_a+1)\
        +self.Lambda_21**self.s_2*x1_a))+(1-self.Lambda_23**self.s_2)**2*x2_a/(self.s_2*(x2_a+self.Lambda_23**self.s_2\
        *(-x2_a-x1_a+1)+self.Lambda_21**self.s_2*x1_a)**2)+1/x2_a+1/(\
        -x2_a-x1_a+1),0,0,0,-1,0]

        Jac6 = [0,0,(self.Lambda_32**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3\
        *x1_b-x1_b+1))+(self.Lambda_31**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-\
        x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1))+(self.Lambda_31**self.s_3-1)*(self.Lambda_32**self.s_3\
        -1)*(-x2_b-x1_b+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3\
        *x1_b-x1_b+1)**2)-(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)/(self.s_1*(self.Lambda_12**self.s_1\
        *x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b))+(1-self.Lambda_13**self.s_1\
        )*(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)*x1_b/(self.s_1*(self.Lambda_12**self.s_1\
        *x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b)**2)-(self.Lambda_21**self.s_2\
        -self.Lambda_23**self.s_2)/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+\
        self.Lambda_21**self.s_2*x1_b))+(1-self.Lambda_23**self.s_2)*(self.Lambda_21**self.s_2-self.Lambda_23**self.s_2\
        )*x2_b/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+self.Lambda_21**self.s_2\
        *x1_b)**2)+1/(-x2_b-x1_b+1),2*(self.Lambda_32**self.s_3-1)/(self.s_3*(self.Lambda_32**self.s_3\
        *x2_b-x2_b+self.Lambda_31**self.s_3*x1_b-x1_b+1))+(self.Lambda_32**self.s_3\
        -1)**2*(-x2_b-x1_b+1)/(self.s_3*(self.Lambda_32**self.s_3*x2_b-x2_b+self.Lambda_31**self.s_3\
        *x1_b-x1_b+1)**2)+(self.Lambda_12**self.s_1-self.Lambda_13**self.s_1)**2*x1_b\
        /(self.s_1*(self.Lambda_12**self.s_1*x2_b+self.Lambda_13**self.s_1*(-x2_b-x1_b+1)+x1_b\
        )**2)-2*(1-self.Lambda_23**self.s_2)/(self.s_2*(x2_b+self.Lambda_23**self.s_2*(-x2_b-x1_b\
        +1)+self.Lambda_21**self.s_2*x1_b))+(1-self.Lambda_23**self.s_2)**2*x2_b/(self.s_2*(x2_b\
        +self.Lambda_23**self.s_2*(-x2_b-x1_b+1)+self.Lambda_21**self.s_2*x1_b)**2)+1/x2_b\
        +1/(-x2_b-x1_b+1),0,-1,0]

        Jac7 = [(n1*x2_b-n1*x2_a+m1*x1_b-m1*x1_a)*(xm_1-x1_a)/(x1_b-x1_a)**2-m1*(\
        xm_1-x1_a)/(x1_b-x1_a)-(n1*x2_b-n1*x2_a+m1*x1_b-m1*x1_a)/(x1_b-\
        x1_a)+m1,n1-n1*(xm_1-x1_a)/(x1_b-x1_a),m1*(xm_1-x1_a)/(x1_b-x1_a\
        )-(n1*x2_b-n1*x2_a+m1*x1_b-m1*x1_a)*(xm_1-x1_a)/(x1_b-x1_a)**2\
        ,n1*(xm_1-x1_a)/(x1_b-x1_a),0,-xm_2+(x2_b-x2_a)*(xm_1-x1_a)/(x1_b\
        -x1_a)+x2_a,0]

        return real(array([Jac1, Jac2, Jac3, Jac4, Jac5, Jac6, Jac7]))
        
        
    
    
    
