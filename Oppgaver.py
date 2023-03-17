import numpy as np
import matplotlib.pyplot as plt

# PARAMETERE

R = 8.314             #[J / mol K]
T_c = 647.096         #[K]
p_c = 22.064 * 10e6   #[Pa] 
#n = 1                 #[mol] 


# OPPGAVE 1a)

a = (27 * R**2 * T_c**2) / (64 * p_c)
b = (R * T_c) / (8 * p_c)
V_c = 3*b

print(a, b, V_c)
#Sammenlign V_c med eksperimentelle verdier
#Eksperimentell V_c = 55.948mL


# OPPGAVE 1b)















































































def ex_eq(c):
    return (np.sinh(2*c/T_c))**2-1

def d_ex_eq(c):
    return 2*np.sinh(4*c/T_c)/T_c


def newton_one_var(f, df, c_0, tol):
    max_it =100000
    i = 0
    x = np.array([c_0])
    while(np.abs(f(c_0))>=tol):
        x_i = x[-1] - (f(x[-1])/df(x[-1]))
        x = np.append(x, x_i)
        i+=1
        if(i>max_it):
            raise Exception("Maximum iterations reached")
    print(f"The root was found to be at {x[-1]} after {i} iterations")
    return x

c_0 = 1
c_arr = newton_one_var(ex_eq, d_ex_eq, c_0, 10e-2)