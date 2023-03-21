"""
Created on Tue Mar 21 08:43:24 2023

@author: hanna
"""
from sympy import *
import numpy as np

"""
a = 5.537                              #[bar * L^2 / mol^2]
a = a * 1000**2                        #[bar * mL^2 / mol^2]
b = 0.03049                            #[L / mol]
b = b * 1000                           


R = 8.314             #[J / mol K]
T = 647.096 
"""


def func1():
    return R*T*(1 / (V_g - b) - 1 / (V_v - b)) - a*(1 / V_g**2 - 1 / V_v**2)  # 11

def func2():
    return R*T*(log((V_g - b) / (V_v - b)) / (V_g - V_v) ) - a*(1 / (V_g*V_v) + 1 / (V_g**2)) # 12

a, b, R, T, V_g, V_v = symbols('a b R T V_g V_v') 

f1 = func1()
f2 = func2()

df1_dg = diff(f1, V_g)
df1_dv = diff(f1, V_v)

df2_dg = diff(f2, V_g)
df2_dv = diff(f2, V_v)

Jacobi_uneval = np.array([[df1_dg, df1_dv], [df2_dg, df2_dv]])


def Jacobi(J, x0, y0):
    
    """
    Parameters
    ----------
    J : Array
        Unevaluated Jacobi determinant.
    x0 : 
    
    y0: 
        
    Returns
    -------
    Jacobi : Array
             Evaluated Jacobi determinant at (x0, y0)
    """
    
    return J.subs([(a, a), (b, b), (R, R), (T, T), (V_g, x0), (V_v, y0)])
    
    
    

