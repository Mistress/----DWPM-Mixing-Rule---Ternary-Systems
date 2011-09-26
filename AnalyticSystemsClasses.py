#!/usr/bin/env python

from scipy import *


class SystemDWPM:
    def __init__(self, Cell_s):
        self.s_1 = Cell_s[0]
        self.s_2 = Cell_s[1]
        self.s_3 = Cell_s[2]

    def SystemEquations(self, Params, C, Actual, R, T):
    
        Lambda_12 = Params[0]/C[0]
        Lambda_21 = Params[1]/C[1]

        Lambda_13 = Params[2]/C[0]
        Lambda_31 = Params[3]/C[2]
        
        Lambda_23 = Params[4]/C[1]
        Lambda_32 = Params[5]/C[2]
                
        m1 = Params[6]
        n1 = Params[7]
        c1 = Params[8]

        m2 = Params[9]
        n2 = Params[10]
        c2 = Params[11]
        
        x1_a1 = Actual[0, 0]
        x2_a1 = Actual[0, 1]
        x1_b1 = Actual[0, 2]
        x2_b1 = Actual[0, 3]
        
        x1_a2 = Actual[1, 0]
        x2_a2 = Actual[1, 1]
        x1_b2 = Actual[1, 2]
        x2_b2 = Actual[1, 3]

        Eq1 = float(-(-x2_a1-x1_a1+1)*log(Lambda_32**self.s_3*x2_a1-x2_a1+Lambda_31**self.s_3*x1_a1-x1_a1+1)/self.s_3-x1_a1*log(Lambda_12**self.s_1*x2_a1+Lambda_13**self.s_1*(-x2_a1-x1_a1+1)+x1_a1)/self.s_1-x2_a1*log(x2_a1+Lambda_23**self.s_2*(-x2_a1-x1_a1+1)+Lambda_21**self.s_2*x1_a1)/self.s_2+x2_a1*log(x2_a1)-n1*x2_a1+log(-x2_a1-x1_a1+1)*(-x2_a1-x1_a1+1)+x1_a1*log(x1_a1)-m1*x1_a1-c1 )

        Eq2 = float(-(-x2_b1-x1_b1+1)*log(Lambda_32**self.s_3*x2_b1-x2_b1+Lambda_31**self.s_3*x1_b1-x1_b1+1)/self.s_3-x1_b1*log(Lambda_12**self.s_1*x2_b1+Lambda_13**self.s_1*(-x2_b1-x1_b1+1)+x1_b1)/self.s_1-x2_b1*log(x2_b1+Lambda_23**self.s_2*(-x2_b1-x1_b1+1)+Lambda_21**self.s_2*x1_b1)/self.s_2+x2_b1*log(x2_b1)-n1*x2_b1+log(-x2_b1-x1_b1+1)*(-x2_b1-x1_b1+1)+x1_b1*log(x1_b1)-m1*x1_b1-c1)

        Eq3 = float(-(-x2_a2-x1_a2+1)*log(Lambda_32**self.s_3*x2_a2-x2_a2+Lambda_31**self.s_3*x1_a2-x1_a2+1)/self.s_3-x1_a2*log(Lambda_12**self.s_1*x2_a2+Lambda_13**self.s_1*(-x2_a2-x1_a2+1)+x1_a2)/self.s_1-x2_a2*log(x2_a2+Lambda_23**self.s_2*(-x2_a2-x1_a2+1)+Lambda_21**self.s_2*x1_a2)/self.s_2+x2_a2*log(x2_a2)-n2*x2_a2+log(-x2_a2-x1_a2+1)*(-x2_a2-x1_a2+1)+x1_a2*log(x1_a2)-m2*x1_a2-c2)
        
        Eq4 = float(-(-x2_b2-x1_b2+1)*log(Lambda_32**self.s_3*x2_b2-x2_b2+Lambda_31**self.s_3*x1_b2-x1_b2+1)/self.s_3-x1_b2*log(Lambda_12**self.s_1*x2_b2+Lambda_13**self.s_1*(-x2_b2-x1_b2+1)+x1_b2)/self.s_1-x2_b2*log(x2_b2+Lambda_23**self.s_2*(-x2_b2-x1_b2+1)+Lambda_21**self.s_2*x1_b2)/self.s_2+x2_b2*log(x2_b2)-n2*x2_b2+log(-x2_b2-x1_b2+1)*(-x2_b2-x1_b2+1)+x1_b2*log(x1_b2)-m2*x1_b2-c2)
        
        Eq5 = float(log(Lambda_32**self.s_3*x2_a1-x2_a1+Lambda_31**self.s_3*x1_a1-x1_a1+1)/self.s_3-log(Lambda_12**self.s_1*x2_a1+Lambda_13**self.s_1*(-x2_a1-x1_a1+1)+x1_a1)/self.s_1-(Lambda_31**self.s_3-1)*(-x2_a1-x1_a1+1)/(self.s_3*(Lambda_32**self.s_3*x2_a1-x2_a1+Lambda_31**self.s_3*x1_a1-x1_a1+1))-(1-Lambda_13**self.s_1)*x1_a1/(self.s_1*(Lambda_12**self.s_1*x2_a1+Lambda_13**self.s_1*(-x2_a1-x1_a1+1)+x1_a1))-(Lambda_21**self.s_2-Lambda_23**self.s_2)*x2_a1/(self.s_2*(x2_a1+Lambda_23**self.s_2*(-x2_a1-x1_a1+1)+Lambda_21**self.s_2*x1_a1))-log(-x2_a1-x1_a1+1)+log(x1_a1)-m1)

        Eq6 = float(log(Lambda_32**self.s_3*x2_b1-x2_b1+Lambda_31**self.s_3*x1_b1-x1_b1+1)/self.s_3-log(Lambda_12**self.s_1*x2_b1+Lambda_13**self.s_1*(-x2_b1-x1_b1+1)+x1_b1)/self.s_1-(Lambda_31**self.s_3-1)*(-x2_b1-x1_b1+1)/(self.s_3*(Lambda_32**self.s_3*x2_b1-x2_b1+Lambda_31**self.s_3*x1_b1-x1_b1+1))-(1-Lambda_13**self.s_1)*x1_b1/(self.s_1*(Lambda_12**self.s_1*x2_b1+Lambda_13**self.s_1*(-x2_b1-x1_b1+1)+x1_b1))-(Lambda_21**self.s_2-Lambda_23**self.s_2)*x2_b1/(self.s_2*(x2_b1+Lambda_23**self.s_2*(-x2_b1-x1_b1+1)+Lambda_21**self.s_2*x1_b1))-log(-x2_b1-x1_b1+1)+log(x1_b1)-m1)

        Eq7 = float(log(Lambda_32**self.s_3*x2_a2-x2_a2+Lambda_31**self.s_3*x1_a2-x1_a2+1)/self.s_3-log(Lambda_12**self.s_1*x2_a2+Lambda_13**self.s_1*(-x2_a2-x1_a2+1)+x1_a2)/self.s_1-(Lambda_31**self.s_3-1)*(-x2_a2-x1_a2+1)/(self.s_3*(Lambda_32**self.s_3*x2_a2-x2_a2+Lambda_31**self.s_3*x1_a2-x1_a2+1))-(1-Lambda_13**self.s_1)*x1_a2/(self.s_1*(Lambda_12**self.s_1*x2_a2+Lambda_13**self.s_1*(-x2_a2-x1_a2+1)+x1_a2))-(Lambda_21**self.s_2-Lambda_23**self.s_2)*x2_a2/(self.s_2*(x2_a2+Lambda_23**self.s_2*(-x2_a2-x1_a2+1)+Lambda_21**self.s_2*x1_a2))-log(-x2_a2-x1_a2+1)+log(x1_a2)-m2)

        Eq8 = float(log(Lambda_32**self.s_3*x2_b2-x2_b2+Lambda_31**self.s_3*x1_b2-x1_b2+1)/self.s_3-log(Lambda_12**self.s_1*x2_b2+Lambda_13**self.s_1*(-x2_b2-x1_b2+1)+x1_b2)/self.s_1-(Lambda_31**self.s_3-1)*(-x2_b2-x1_b2+1)/(self.s_3*(Lambda_32**self.s_3*x2_b2-x2_b2+Lambda_31**self.s_3*x1_b2-x1_b2+1))-(1-Lambda_13**self.s_1)*x1_b2/(self.s_1*(Lambda_12**self.s_1*x2_b2+Lambda_13**self.s_1*(-x2_b2-x1_b2+1)+x1_b2))-(Lambda_21**self.s_2-Lambda_23**self.s_2)*x2_b2/(self.s_2*(x2_b2+Lambda_23**self.s_2*(-x2_b2-x1_b2+1)+Lambda_21**self.s_2*x1_b2))-log(-x2_b2-x1_b2+1)+log(x1_b2)-m2)

        Eq9 = float(log(Lambda_32**self.s_3*x2_a1-x2_a1+Lambda_31**self.s_3*x1_a1-x1_a1+1)/self.s_3-log(x2_a1+Lambda_23**self.s_2*(-x2_a1-x1_a1+1)+Lambda_21**self.s_2*x1_a1)/self.s_2+log(x2_a1)-(Lambda_32**self.s_3-1)*(-x2_a1-x1_a1+1)/(self.s_3*(Lambda_32**self.s_3*x2_a1-x2_a1+Lambda_31**self.s_3*x1_a1-x1_a1+1))-(Lambda_12**self.s_1-Lambda_13**self.s_1)*x1_a1/(self.s_1*(Lambda_12**self.s_1*x2_a1+Lambda_13**self.s_1*(-x2_a1-x1_a1+1)+x1_a1))-(1-Lambda_23**self.s_2)*x2_a1/(self.s_2*(x2_a1+Lambda_23**self.s_2*(-x2_a1-x1_a1+1)+Lambda_21**self.s_2*x1_a1))-log(-x2_a1-x1_a1+1)-n1)

        Eq10 = float(log(Lambda_32**self.s_3*x2_b1-x2_b1+Lambda_31**self.s_3*x1_b1-x1_b1+1)/self.s_3-log(x2_b1+Lambda_23**self.s_2*(-x2_b1-x1_b1+1)+Lambda_21**self.s_2*x1_b1)/self.s_2+log(x2_b1)-(Lambda_32**self.s_3-1)*(-x2_b1-x1_b1+1)/(self.s_3*(Lambda_32**self.s_3*x2_b1-x2_b1+Lambda_31**self.s_3*x1_b1-x1_b1+1))-(Lambda_12**self.s_1-Lambda_13**self.s_1)*x1_b1/(self.s_1*(Lambda_12**self.s_1*x2_b1+Lambda_13**self.s_1*(-x2_b1-x1_b1+1)+x1_b1))-(1-Lambda_23**self.s_2)*x2_b1/(self.s_2*(x2_b1+Lambda_23**self.s_2*(-x2_b1-x1_b1+1)+Lambda_21**self.s_2*x1_b1))-log(-x2_b1-x1_b1+1)-n1)

        Eq11 = float(log(Lambda_32**self.s_3*x2_a2-x2_a2+Lambda_31**self.s_3*x1_a2-x1_a2+1)/self.s_3-log(x2_a2+Lambda_23**self.s_2*(-x2_a2-x1_a2+1)+Lambda_21**self.s_2*x1_a2)/self.s_2+log(x2_a2)-(Lambda_32**self.s_3-1)*(-x2_a2-x1_a2+1)/(self.s_3*(Lambda_32**self.s_3*x2_a2-x2_a2+Lambda_31**self.s_3*x1_a2-x1_a2+1))-(Lambda_12**self.s_1-Lambda_13**self.s_1)*x1_a2/(self.s_1*(Lambda_12**self.s_1*x2_a2+Lambda_13**self.s_1*(-x2_a2-x1_a2+1)+x1_a2))-(1-Lambda_23**self.s_2)*x2_a2/(self.s_2*(x2_a2+Lambda_23**self.s_2*(-x2_a2-x1_a2+1)+Lambda_21**self.s_2*x1_a2))-log(-x2_a2-x1_a2+1)-n2)

        Eq12 = float(log(Lambda_32**self.s_3*x2_b2-x2_b2+Lambda_31**self.s_3*x1_b2-x1_b2+1)/self.s_3-log(x2_b2+Lambda_23**self.s_2*(-x2_b2-x1_b2+1)+Lambda_21**self.s_2*x1_b2)/self.s_2+log(x2_b2)-(Lambda_32**self.s_3-1)*(-x2_b2-x1_b2+1)/(self.s_3*(Lambda_32**self.s_3*x2_b2-x2_b2+Lambda_31**self.s_3*x1_b2-x1_b2+1))-(Lambda_12**self.s_1-Lambda_13**self.s_1)*x1_b2/(self.s_1*(Lambda_12**self.s_1*x2_b2+Lambda_13**self.s_1*(-x2_b2-x1_b2+1)+x1_b2))-(1-Lambda_23**self.s_2)*x2_b2/(self.s_2*(x2_b2+Lambda_23**self.s_2*(-x2_b2-x1_b2+1)+Lambda_21**self.s_2*x1_b2))-log(-x2_b2-x1_b2+1)-n2)
       
        return array([Eq1, Eq2, Eq3, Eq4, Eq5, Eq6, Eq7, Eq8, Eq9, Eq10, Eq11, Eq12])

    def SystemEquationsJac(self, Params, C, Actual, R, T):

        c_12 = Params[0]
        c_21 = Params[1]

        c_13 = Params[2]
        c_31 = Params[3]
        
        c_23 = Params[4]
        c_32 = Params[5]

        c_1 = float(C[0])
        c_2 = float(C[1])
        c_3 = float(C[2])
                
        m1 = Params[6]
        n1 = Params[7]
        c1 = Params[8]

        m2 = Params[9]
        n2 = Params[10]
        c2 = Params[11]
        
        x1_a1 = Actual[0, 0]
        x2_a1 = Actual[0, 1]
        x1_b1 = Actual[0, 2]
        x2_b1 = Actual[0, 3]
        
        x1_a2 = Actual[1, 0]
        x2_a2 = Actual[1, 1]
        x1_b2 = Actual[1, 2]
        x2_b2 = Actual[1, 3]

        Row1 = [-(c_12/c_1)**self.s_1*x1_a1*x2_a1/(c_12*((c_12/c_1)**self.s_1*x2_a1+(c_13/\
        c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1)),-(c_21/c_2)**self.s_2*x1_a1*x2_a1/\
        (c_21*(x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1+1)+(c_21/c_2)**self.s_2*\
        x1_a1)),-(c_13/c_1)**self.s_1*x1_a1*(-x2_a1-x1_a1+1)/(c_13*((c_12/c_1\
        )**self.s_1*x2_a1+(c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1)),-(c_31/c_3\
        )**self.s_3*x1_a1*(-x2_a1-x1_a1+1)/(c_31*((c_32/c_3)**self.s_3*x2_a1-x2_a1\
        +(c_31/c_3)**self.s_3*x1_a1-x1_a1+1)),-(c_23/c_2)**self.s_2*(-x2_a1-x1_a1\
        +1)*x2_a1/(c_23*(x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1+1)+(c_21\
        /c_2)**self.s_2*x1_a1)),-(c_32/c_3)**self.s_3*(-x2_a1-x1_a1+1)*x2_a1/(c_32\
        *((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/c_3)**self.s_3*x1_a1-x1_a1+1)),\
        -x1_a1,-x2_a1,-1,0,0,0]

        Row2 = [-(c_12/c_1)**self.s_1*x1_b1*x2_b1/(c_12*((c_12/c_1)**self.s_1*x2_b1+(c_13/\
        c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1)),-(c_21/c_2)**self.s_2*x1_b1*x2_b1/\
        (c_21*(x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1+1)+(c_21/c_2)**self.s_2*\
        x1_b1)),-(c_13/c_1)**self.s_1*x1_b1*(-x2_b1-x1_b1+1)/(c_13*((c_12/c_1\
        )**self.s_1*x2_b1+(c_13/c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1)),-(c_31/c_3\
        )**self.s_3*x1_b1*(-x2_b1-x1_b1+1)/(c_31*((c_32/c_3)**self.s_3*x2_b1-x2_b1\
        +(c_31/c_3)**self.s_3*x1_b1-x1_b1+1)),-(c_23/c_2)**self.s_2*(-x2_b1-x1_b1\
        +1)*x2_b1/(c_23*(x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1+1)+(c_21\
        /c_2)**self.s_2*x1_b1)),-(c_32/c_3)**self.s_3*(-x2_b1-x1_b1+1)*x2_b1/(c_32\
        *((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/c_3)**self.s_3*x1_b1-x1_b1+1)),\
        -x1_b1,-x2_b1,-1,0,0,0]

        Row3 = [-(c_12/c_1)**self.s_1*x1_a2*x2_a2/(c_12*((c_12/c_1)**self.s_1*x2_a2+(c_13/\
        c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2)),-(c_21/c_2)**self.s_2*x1_a2*x2_a2/\
        (c_21*(x2_a2+(c_23/c_2)**self.s_2*(-x2_a2-x1_a2+1)+(c_21/c_2)**self.s_2*\
        x1_a2)),-(c_13/c_1)**self.s_1*x1_a2*(-x2_a2-x1_a2+1)/(c_13*((c_12/c_1\
        )**self.s_1*x2_a2+(c_13/c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2)),-(c_31/c_3\
        )**self.s_3*x1_a2*(-x2_a2-x1_a2+1)/(c_31*((c_32/c_3)**self.s_3*x2_a2-x2_a2\
        +(c_31/c_3)**self.s_3*x1_a2-x1_a2+1)),-(c_23/c_2)**self.s_2*(-x2_a2-x1_a2\
        +1)*x2_a2/(c_23*(x2_a2+(c_23/c_2)**self.s_2*(-x2_a2-x1_a2+1)+(c_21\
        /c_2)**self.s_2*x1_a2)),-(c_32/c_3)**self.s_3*(-x2_a2-x1_a2+1)*x2_a2/(c_32\
        *((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/c_3)**self.s_3*x1_a2-x1_a2+1)),\
        0,0,0,-x1_a2,-x2_a2,-1]

        Row4 = [-(c_12/c_1)**self.s_1*x1_b2*x2_b2/(c_12*((c_12/c_1)**self.s_1*x2_b2+(c_13/\
        c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2)),-(c_21/c_2)**self.s_2*x1_b2*x2_b2/\
        (c_21*(x2_b2+(c_23/c_2)**self.s_2*(-x2_b2-x1_b2+1)+(c_21/c_2)**self.s_2*\
        x1_b2)),-(c_13/c_1)**self.s_1*x1_b2*(-x2_b2-x1_b2+1)/(c_13*((c_12/c_1\
        )**self.s_1*x2_b2+(c_13/c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2)),-(c_31/c_3\
        )**self.s_3*x1_b2*(-x2_b2-x1_b2+1)/(c_31*((c_32/c_3)**self.s_3*x2_b2-x2_b2\
        +(c_31/c_3)**self.s_3*x1_b2-x1_b2+1)),-(c_23/c_2)**self.s_2*(-x2_b2-x1_b2\
        +1)*x2_b2/(c_23*(x2_b2+(c_23/c_2)**self.s_2*(-x2_b2-x1_b2+1)+(c_21\
        /c_2)**self.s_2*x1_b2)),-(c_32/c_3)**self.s_3*(-x2_b2-x1_b2+1)*x2_b2/(c_32\
        *((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/c_3)**self.s_3*x1_b2-x1_b2+1)),\
        0,0,0,-x1_b2,-x2_b2,-1]

        Row5 = [(c_12/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_a1*x2_a1/(c_12*((c_12/c_1)**self.s_1\
        *x2_a1+(c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1)**2)-(c_12/\
        c_1)**self.s_1*x2_a1/(c_12*((c_12/c_1)**self.s_1*x2_a1+(c_13/c_1)**self.s_1*(\
        -x2_a1-x1_a1+1)+x1_a1)),(c_21/c_2)**self.s_2*((c_21/c_2)**self.s_2-(c_23/c_2\
        )**self.s_2)*x1_a1*x2_a1/(c_21*(x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1\
        +1)+(c_21/c_2)**self.s_2*x1_a1)**2)-(c_21/c_2)**self.s_2*x2_a1/(c_21*(\
        x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1+1)+(c_21/c_2)**self.s_2*x1_a1)),-(\
        c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)/(c_13*((c_12/c_1)**self.s_1*x2_a1+(\
        c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1))+(c_13/c_1)**self.s_1*x1_a1/(\
        c_13*((c_12/c_1)**self.s_1*x2_a1+(c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1\
        ))+(c_13/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_a1*(-x2_a1-x1_a1+1)/\
        (c_13*((c_12/c_1)**self.s_1*x2_a1+(c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)+\
        x1_a1)**2),-(c_31/c_3)**self.s_3*(-x2_a1-x1_a1+1)/(c_31*((c_32/c_3)**self.s_3\
        *x2_a1-x2_a1+(c_31/c_3)**self.s_3*x1_a1-x1_a1+1))+(c_31/c_3)**self.s_3\
        *x1_a1/(c_31*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/c_3)**self.s_3*x1_a1\
        -x1_a1+1))+(c_31/c_3)**self.s_3*((c_31/c_3)**self.s_3-1)*x1_a1*(-x2_a1-\
        x1_a1+1)/(c_31*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/c_3)**self.s_3*x1_a1\
        -x1_a1+1)**2),(c_23/c_2)**self.s_2*x2_a1/(c_23*(x2_a1+(c_23/c_2)**self.s_2\
        *(-x2_a1-x1_a1+1)+(c_21/c_2)**self.s_2*x1_a1))+(c_23/c_2)**self.s_2*((c_21\
        /c_2)**self.s_2-(c_23/c_2)**self.s_2)*(-x2_a1-x1_a1+1)*x2_a1/(c_23*(x2_a1\
        +(c_23/c_2)**self.s_2*(-x2_a1-x1_a1+1)+(c_21/c_2)**self.s_2*x1_a1)**2),\
        (c_32/c_3)**self.s_3*x2_a1/(c_32*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/\
        c_3)**self.s_3*x1_a1-x1_a1+1))+((c_31/c_3)**self.s_3-1)*(c_32/c_3)**self.s_3*(\
        -x2_a1-x1_a1+1)*x2_a1/(c_32*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/\
        c_3)**self.s_3*x1_a1-x1_a1+1)**2),-1,0,0,0,0,0]

        Row6 = [(c_12/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_b1*x2_b1/(c_12*((c_12/c_1)**self.s_1\
        *x2_b1+(c_13/c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1)**2)-(c_12/\
        c1)**self.s_1*x2_b1/(c_12*((c_12/c_1)**self.s_1*x2_b1+(c_13/c_1)**self.s_1*(\
        -x2_b1-x1_b1+1)+x1_b1)),(c_21/c_2)**self.s_2*((c_21/c_2)**self.s_2-(c_23/\
        c_2)**self.s_2)*x1_b1*x2_b1/(c_21*(x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1\
        +1)+(c_21/c_2)**self.s_2*x1_b1)**2)-(c_21/c_2)**self.s_2*x2_b1/(c_21*(\
        x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1+1)+(c_21/c_2)**self.s_2*x1_b1)),-(\
        c_13/c_1)**self.s_1*(-x2_b1-x1_b1+1)/(c_13*((c_12/c_1)**self.s_1*x2_b1+(c_13\
        /c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1))+(c_13/c_1)**self.s_1*x1_b1/(c_13\
        *((c_12/c_1)**self.s_1*x2_b1+(c_13/c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1\
        ))+(c_13/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_b1*(-x2_b1-x1_b1+1)/\
        (c_13*((c_12/c_1)**self.s_1*x2_b1+(c_13/c_1)**self.s_1*(-x2_b1-x1_b1+1)+\
        x1_b1)**2),-(c_31/c_3)**self.s_3*(-x2_b1-x1_b1+1)/(c_31*((c_32/c_3)**self.s_3\
        *x2_b1-x2_b1+(c_31/c_3)**self.s_3*x1_b1-x1_b1+1))+(c_31/c_3)**self.s_3\
        *x1_b1/(c_31*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/c_3)**self.s_3*x1_b1\
        -x1_b1+1))+(c_31/c_3)**self.s_3*((c_31/c_3)**self.s_3-1)*x1_b1*(-x2_b1-\
        x1_b1+1)/(c_31*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/c_3)**self.s_3*x1_b1\
        -x1_b1+1)**2),(c_23/c_2)**self.s_2*x2_b1/(c_23*(x2_b1+(c_23/c_2)**self.s_2\
        *(-x2_b1-x1_b1+1)+(c_21/c_2)**self.s_2*x1_b1))+(c_23/c_2)**self.s_2*((c_21\
        /c_2)**self.s_2-(c_23/c_2)**self.s_2)*(-x2_b1-x1_b1+1)*x2_b1/(c_23*(x2_b1\
        +(c_23/c_2)**self.s_2*(-x2_b1-x1_b1+1)+(c_21/c_2)**self.s_2*x1_b1)**2),\
        (c_32/c_3)**self.s_3*x2_b1/(c_32*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/\
        c_3)**self.s_3*x1_b1-x1_b1+1))+((c_31/c_3)**self.s_3-1)*(c_32/c_3)**self.s_3*(\
        -x2_b1-x1_b1+1)*x2_b1/(c_32*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/\
        c_3)**self.s_3*x1_b1-x1_b1+1)**2),-1,0,0,0,0,0]

        Row7 = [(c_12/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_a2*x2_a2/(c_12*((c_12/c_1)**self.s_1\
        *x2_a2+(c_13/c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2)**2)-(c_12/c_1\
        )**self.s_1*x2_a2/(c_12*((c_12/c_1)**self.s_1*x2_a2+(c_13/c_1)**self.s_1*(\
        -x2_a2-x1_a2+1)+x1_a2)),(c_21/c_2)**self.s_2*((c_21/c_2)**self.s_2-(c_23/c_2\
        )**self.s_2)*x1_a2*x2_a2/(c_21*(x2_a2+(c_23/c_2)**self.s_2*(-x2_a2-x1_a2\
        +1)+(c_21/c_2)**self.s_2*x1_a2)**2)-(c_21/c_2)**self.s_2*x2_a2/(c_21*(x2_a2\
        +(c_23/c_2)**self.s_2*(-x2_a2-x1_a2+1)+(c_21/c_2)**self.s_2*x1_a2)),-(\
        c_13/c_1)**self.s_1*(-x2_a2-x1_a2+1)/(c_13*((c_12/c_1)**self.s_1*x2_a2+(c_13\
        /c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2))+(c_13/c_1)**self.s_1*x1_a2/(c_13\
        *((c_12/c_1)**self.s_1*x2_a2+(c_13/c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2\
        ))+(c_13/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_a2*(-x2_a2-x1_a2+1)/\
        (c_13*((c_12/c_1)**self.s_1*x2_a2+(c_13/c_1)**self.s_1*(-x2_a2-x1_a2+1)+\
        x1_a2)**2),-(c_31/c_3)**self.s_3*(-x2_a2-x1_a2+1)/(c_31*((c_32/c_3)**self.s_3\
        *x2_a2-x2_a2+(c_31/c_3)**self.s_3*x1_a2-x1_a2+1))+(c_31/c_3)**self.s_3\
        *x1_a2/(c_31*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/c_3)**self.s_3*x1_a2\
        -x1_a2+1))+(c_31/c_3)**self.s_3*((c_31/c_3)**self.s_3-1)*x1_a2*(-x2_a2-\
        x1_a2+1)/(c_31*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/c_3)**self.s_3*x1_a2\
        -x1_a2+1)**2),(c_23/c_2)**self.s_2*x2_a2/(c_23*(x2_a2+(c_23/c_2)**self.s_2\
        *(-x2_a2-x1_a2+1)+(c_21/c_2)**self.s_2*x1_a2))+(c_23/c_2)**self.s_2*((c_21\
        /c_2)**self.s_2-(c_23/c_2)**self.s_2)*(-x2_a2-x1_a2+1)*x2_a2/(c_23*(x2_a2\
        +(c_23/c_2)**self.s_2*(-x2_a2-x1_a2+1)+(c_21/c_2)**self.s_2*x1_a2)**2),\
        (c_32/c_3)**self.s_3*x2_a2/(c_32*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/\
        c_3)**self.s_3*x1_a2-x1_a2+1))+((c_31/c_3)**self.s_3-1)*(c_32/c_3)**self.s_3*(\
        -x2_a2-x1_a2+1)*x2_a2/(c_32*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/\
        c_3)**self.s_3*x1_a2-x1_a2+1)**2),0,0,0,-1,0,0]

        Row8 = [(c_12/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_b2*x2_b2/(c_12*((c_12/c_1)**self.s_1\
        *x2_b2+(c_13/c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2)**2)-(c_12/\
        c_1)**self.s_1*x2_b2/(c_12*((c_12/c_1)**self.s_1*x2_b2+(c_13/c_1)**self.s_1*(\
        -x2_b2-x1_b2+1)+x1_b2)),(c_21/c_2)**self.s_2*((c_21/c_2)**self.s_2-(c_23/c_2\
        )**self.s_2)*x1_b2*x2_b2/(c_21*(x2_b2+(c_23/c_2)**self.s_2*(-x2_b2-x1_b2\
        +1)+(c_21/c_2)**self.s_2*x1_b2)**2)-(c_21/c_2)**self.s_2*x2_b2/(c_21*(x2_b2\
        +(c_23/c_2)**self.s_2*(-x2_b2-x1_b2+1)+(c_21/c_2)**self.s_2*x1_b2)),-(c_13\
        /c_1)**self.s_1*(-x2_b2-x1_b2+1)/(c_13*((c_12/c_1)**self.s_1*x2_b2+(c_13\
        /c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2))+(c_13/c_1)**self.s_1*x1_b2/(c_13\
        *((c_12/c_1)**self.s_1*x2_b2+(c_13/c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2\
        ))+(c_13/c_1)**self.s_1*(1-(c_13/c_1)**self.s_1)*x1_b2*(-x2_b2-x1_b2+1)/\
        (c_13*((c_12/c_1)**self.s_1*x2_b2+(c_13/c_1)**self.s_1*(-x2_b2-x1_b2+1)+\
        x1_b2)**2),-(c_31/c_3)**self.s_3*(-x2_b2-x1_b2+1)/(c_31*((c_32/c_3)**self.s_3\
        *x2_b2-x2_b2+(c_31/c_3)**self.s_3*x1_b2-x1_b2+1))+(c_31/c_3)**self.s_3\
        *x1_b2/(c_31*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/c_3)**self.s_3*x1_b2\
        -x1_b2+1))+(c_31/c_3)**self.s_3*((c_31/c_3)**self.s_3-1)*x1_b2*(-x2_b2-x1_b2\
        +1)/(c_31*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/c_3)**self.s_3*x1_b2\
        -x1_b2+1)**2),(c_23/c_2)**self.s_2*x2_b2/(c_23*(x2_b2+(c_23/c_2)**self.s_2\
        *(-x2_b2-x1_b2+1)+(c_21/c_2)**self.s_2*x1_b2))+(c_23/c_2)**self.s_2*((c_21\
        /c_2)**self.s_2-(c_23/c_2)**self.s_2)*(-x2_b2-x1_b2+1)*x2_b2/(c_23*(x2_b2\
        +(c_23/c_2)**self.s_2*(-x2_b2-x1_b2+1)+(c_21/c_2)**self.s_2*x1_b2)**2),\
        (c_32/c_3)**self.s_3*x2_b2/(c_32*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/\
        c_3)**self.s_3*x1_b2-x1_b2+1))+((c_31/c_3)**self.s_3-1)*(c_32/c_3)**self.s_3*(\
        -x2_b2-x1_b2+1)*x2_b2/(c_32*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/\
        c_3)**self.s_3*x1_b2-x1_b2+1)**2),0,0,0,-1,0,0]

        Row9 =  [(c_12/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13/c_1)**self.s_1)*x1_a1*x2_a1/(c_12\
        *((c_12/c_1)**self.s_1*x2_a1+(c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1\
        )**2)-(c_12/c_1)**self.s_1*x1_a1/(c_12*((c_12/c_1)**self.s_1*x2_a1+(c_13\
        /c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1)),(c_21/c_2)**self.s_2*(1-(c_23/c_2\
        )**self.s_2)*x1_a1*x2_a1/(c_21*(x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1\
        +1)+(c_21/c_2)**self.s_2*x1_a1)**2)-(c_21/c_2)**self.s_2*x1_a1/(c_21*(\
        x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1+1)+(c_21/c_2)**self.s_2*x1_a1)),(\
        c_13/c_1)**self.s_1*x1_a1/(c_13*((c_12/c_1)**self.s_1*x2_a1+(c_13/c_1)**self.s_1\
        *(-x2_a1-x1_a1+1)+x1_a1))+(c_13/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13\
        /c_1)**self.s_1)*x1_a1*(-x2_a1-x1_a1+1)/(c_13*((c_12/c_1)**self.s_1*x2_a1\
        +(c_13/c_1)**self.s_1*(-x2_a1-x1_a1+1)+x1_a1)**2),(c_31/c_3)**self.s_3*\
        x1_a1/(c_31*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/c_3)**self.s_3*x1_a1-\
        x1_a1+1))+(c_31/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*x1_a1*(-x2_a1-x1_a1\
        +1)/(c_31*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/c_3)**self.s_3*x1_a1-\
        x1_a1+1)**2),(c_23/c_2)**self.s_2*x2_a1/(c_23*(x2_a1+(c_23/c_2)**self.s_2*\
        (-x2_a1-x1_a1+1)+(c_21/c_2)**self.s_2*x1_a1))-(c_23/c_2)**self.s_2*(-x2_a1\
        -x1_a1+1)/(c_23*(x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1+1)+(c_21/\
        c_2)**self.s_2*x1_a1))+(c_23/c_2)**self.s_2*(1-(c_23/c_2)**self.s_2)*(-x2_a1-\
        x1_a1+1)*x2_a1/(c_23*(x2_a1+(c_23/c_2)**self.s_2*(-x2_a1-x1_a1+1)+(c_21\
        /c_2)**self.s_2*x1_a1)**2),(c_32/c_3)**self.s_3*x2_a1/(c_32*((c_32/c_3)**self.s_3\
        *x2_a1-x2_a1+(c_31/c_3)**self.s_3*x1_a1-x1_a1+1))-(c_32/c_3)**self.s_3\
        *(-x2_a1-x1_a1+1)/(c_32*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/c_3\
        )**self.s_3*x1_a1-x1_a1+1))+(c_32/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*(-\
        x2_a1-x1_a1+1)*x2_a1/(c_32*((c_32/c_3)**self.s_3*x2_a1-x2_a1+(c_31/c_3\
        )**self.s_3*x1_a1-x1_a1+1)**2),0,-1,0,0,0,0]

        Row10 = [(c_12/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13/c_1)**self.s_1)*x1_b1*x2_b1/(c_12\
        *((c_12/c_1)**self.s_1*x2_b1+(c_13/c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1\
        )**2)-(c_12/c_1)**self.s_1*x1_b1/(c_12*((c_12/c_1)**self.s_1*x2_b1+(c_13\
        /c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1)),(c_21/c_2)**self.s_2*(1-(c_23/c_2\
        )**self.s_2)*x1_b1*x2_b1/(c_21*(x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1\
        +1)+(c_21/c_2)**self.s_2*x1_b1)**2)-(c_21/c_2)**self.s_2*x1_b1/(c_21*(\
        x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1+1)+(c_21/c_2)**self.s_2*x1_b1)),(c_13\
        /c_1)**self.s_1*x1_b1/(c_13*((c_12/c_1)**self.s_1*x2_b1+(c_13/c_1)**self.s_1\
        *(-x2_b1-x1_b1+1)+x1_b1))+(c_13/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13\
        /c_1)**self.s_1)*x1_b1*(-x2_b1-x1_b1+1)/(c_13*((c_12/c_1)**self.s_1*x2_b1\
        +(c_13/c_1)**self.s_1*(-x2_b1-x1_b1+1)+x1_b1)**2),(c_31/c_3)**self.s_3*\
        x1_b1/(c_31*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/c_3)**self.s_3*x1_b1-\
        x1_b1+1))+(c_31/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*x1_b1*(-x2_b1-x1_b1\
        +1)/(c_31*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/c_3)**self.s_3*x1_b1-\
        x1_b1+1)**2),(c_23/c_2)**self.s_2*x2_b1/(c_23*(x2_b1+(c_23/c_2)**self.s_2*\
        (-x2_b1-x1_b1+1)+(c_21/c_2)**self.s_2*x1_b1))-(c_23/c_2)**self.s_2*(-x2_b1\
        -x1_b1+1)/(c_23*(x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1+1)+(c_21/\
        c_2)**self.s_2*x1_b1))+(c_23/c_2)**self.s_2*(1-(c_23/c_2)**self.s_2)*(-x2_b1-\
        x1_b1+1)*x2_b1/(c_23*(x2_b1+(c_23/c_2)**self.s_2*(-x2_b1-x1_b1+1)+(c_21\
        /c_2)**self.s_2*x1_b1)**2),(c_32/c_3)**self.s_3*x2_b1/(c_32*((c_32/c_3)**self.s_3\
        *x2_b1-x2_b1+(c_31/c_3)**self.s_3*x1_b1-x1_b1+1))-(c_32/c_3)**self.s_3\
        *(-x2_b1-x1_b1+1)/(c_32*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/c_3\
        )**self.s_3*x1_b1-x1_b1+1))+(c_32/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*(\
        -x2_b1-x1_b1+1)*x2_b1/(c_32*((c_32/c_3)**self.s_3*x2_b1-x2_b1+(c_31/c_3\
        )**self.s_3*x1_b1-x1_b1+1)**2),0,-1,0,0,0,0]

        Row11 = [(c_12/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13/c_1)**self.s_1)*x1_a2*x2_a2/(c_12\
        *((c_12/c_1)**self.s_1*x2_a2+(c_13/c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2\
        )**2)-(c_12/c_1)**self.s_1*x1_a2/(c_12*((c_12/c_1)**self.s_1*x2_a2+(c_13\
        /c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2)),(c_21/c_2)**self.s_2*(1-(c_23/c_2\
        )**self.s_2)*x1_a2*x2_a2/(c_21*(x2_a2+(c_23/c_2)**self.s_2*(-x2_a2-x1_a2\
        +1)+(c_21/c_2)**self.s_2*x1_a2)**2)-(c_21/c_2)**self.s_2*x1_a2/(c_21*(x2_a2\
        +(c_23/c_2)**self.s_2*(-x2_a2-x1_a2+1)+(c_21/c_2)**self.s_2*x1_a2)),(c_13\
        /c_1)**self.s_1*x1_a2/(c_13*((c_12/c_1)**self.s_1*x2_a2+(c_13/c_1)**self.s_1\
        *(-x2_a2-x1_a2+1)+x1_a2))+(c_13/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13\
        /c_1)**self.s_1)*x1_a2*(-x2_a2-x1_a2+1)/(c_13*((c_12/c_1)**self.s_1*x2_a2\
        +(c_13/c_1)**self.s_1*(-x2_a2-x1_a2+1)+x1_a2)**2),(c_31/c_3)**self.s_3*\
        x1_a2/(c_31*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/c_3)**self.s_3*x1_a2-\
        x1_a2+1))+(c_31/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*x1_a2*(-x2_a2-x1_a2\
        +1)/(c_31*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/c_3)**self.s_3*x1_a2-\
        x1_a2+1)**2),(c_23/c_2)**self.s_2*x2_a2/(c_23*(x2_a2+(c_23/c_2)**self.s_2*\
        (-x2_a2-x1_a2+1)+(c_21/c_2)**self.s_2*x1_a2))-(c_23/c_2)**self.s_2*(-x2_a2\
        -x1_a2+1)/(c_23*(x2_a2+(c_23/c_2)**self.s_2*(-x2_a2-x1_a2+1)+(c_21/\
        c_2)**self.s_2*x1_a2))+(c_23/c_2)**self.s_2*(1-(c_23/c_2)**self.s_2)*(-x2_a2-\
        x1_a2+1)*x2_a2/(c_23*(x2_a2+(c_23/c_2)**self.s_2*(-x2_a2-x1_a2+1)+(c_21\
        /c_2)**self.s_2*x1_a2)**2),(c_32/c_3)**self.s_3*x2_a2/(c_32*((c_32/c_3)**self.s_3\
        *x2_a2-x2_a2+(c_31/c_3)**self.s_3*x1_a2-x1_a2+1))-(c_32/c_3)**self.s_3\
        *(-x2_a2-x1_a2+1)/(c_32*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/c_3\
        )**self.s_3*x1_a2-x1_a2+1))+(c_32/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*(\
        -x2_a2-x1_a2+1)*x2_a2/(c_32*((c_32/c_3)**self.s_3*x2_a2-x2_a2+(c_31/c_3\
        )**self.s_3*x1_a2-x1_a2+1)**2),0,0,0,0,-1,0]

        Row12 = [(c_12/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13/c_1)**self.s_1)*x1_b2*x2_b2/(c_12\
        *((c_12/c_1)**self.s_1*x2_b2+(c_13/c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2\
        )**2)-(c_12/c_1)**self.s_1*x1_b2/(c_12*((c_12/c_1)**self.s_1*x2_b2+(c_13\
        /c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2)),(c_21/c_2)**self.s_2*(1-(c_23/c_2\
        )**self.s_2)*x1_b2*x2_b2/(c_21*(x2_b2+(c_23/c_2)**self.s_2*(-x2_b2-x1_b2\
        +1)+(c_21/c_2)**self.s_2*x1_b2)**2)-(c_21/c_2)**self.s_2*x1_b2/(c_21*(x2_b2\
        +(c_23/c_2)**self.s_2*(-x2_b2-x1_b2+1)+(c_21/c_2)**self.s_2*x1_b2)),(c_13\
        /c_1)**self.s_1*x1_b2/(c_13*((c_12/c_1)**self.s_1*x2_b2+(c_13/c_1)**self.s_1\
        *(-x2_b2-x1_b2+1)+x1_b2))+(c_13/c_1)**self.s_1*((c_12/c_1)**self.s_1-(c_13\
        /c_1)**self.s_1)*x1_b2*(-x2_b2-x1_b2+1)/(c_13*((c_12/c_1)**self.s_1*x2_b2\
        +(c_13/c_1)**self.s_1*(-x2_b2-x1_b2+1)+x1_b2)**2),(c_31/c_3)**self.s_3*\
        x1_b2/(c_31*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/c_3)**self.s_3*x1_b2-\
        x1_b2+1))+(c_31/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*x1_b2*(-x2_b2-x1_b2\
        +1)/(c_31*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/c_3)**self.s_3*x1_b2-\
        x1_b2+1)**2),(c_23/c_2)**self.s_2*x2_b2/(c_23*(x2_b2+(c_23/c_2)**self.s_2*\
        (-x2_b2-x1_b2+1)+(c_21/c_2)**self.s_2*x1_b2))-(c_23/c_2)**self.s_2*(-x2_b2\
        -x1_b2+1)/(c_23*(x2_b2+(c_23/c_2)**self.s_2*(-x2_b2-x1_b2+1)+(c_21/\
        c_2)**self.s_2*x1_b2))+(c_23/c_2)**self.s_2*(1-(c_23/c_2)**self.s_2)*(-x2_b2-\
        x1_b2+1)*x2_b2/(c_23*(x2_b2+(c_23/c_2)**self.s_2*(-x2_b2-x1_b2+1)+(c_21\
        /c_2)**self.s_2*x1_b2)**2),(c_32/c_3)**self.s_3*x2_b2/(c_32*((c_32/c_3)**self.s_3\
        *x2_b2-x2_b2+(c_31/c_3)**self.s_3*x1_b2-x1_b2+1))-(c_32/c_3)**self.s_3\
        *(-x2_b2-x1_b2+1)/(c_32*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/c_3
        )**self.s_3*x1_b2-x1_b2+1))+(c_32/c_3)**self.s_3*((c_32/c_3)**self.s_3-1)*(\
        -x2_b2-x1_b2+1)*x2_b2/(c_32*((c_32/c_3)**self.s_3*x2_b2-x2_b2+(c_31/c_3\
        )**self.s_3*x1_b2-x1_b2+1)**2),0,0,0,0,-1,0]               

        return array([Row1, Row2, Row3, Row4, Row5, Row6, Row7, Row8, Row9, Row10, Row11, Row12])

