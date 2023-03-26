"""
Created on Fri Mar 24 16:15:59 2023

@author: hanna
"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy import optimize

# Oppgave 2a

# Vi skal bruke flere sett med eksperimentelle verdier for å utvikle en modell som senker dette avviket

file1 = "L_verdier.txt"

# Leser datafilen
L_data = np.loadtxt(file1)
T1 = L_data[:, 0] + 273.15               # Temperaturer [K]
L = L_data[:, 2]                         # Heat of vaporization [J/mol]
L = L * 10                               # Heat of vaporization [bar * mL / mol]

file2 = "Vg_verdier.txt"

# Leser datafilen
Vg_data = np.loadtxt(file2)
T2 = Vg_data[:, 0]                    # Temperaturer [K]
Vg = (18.01528/Vg_data[:, 3])*1000   # Gassvolum [mL]

file3 = "Vv_verdier.txt"

# Leser datafilen
Vv_data = np.loadtxt(file3)
T3 = Vv_data[:, 0]                    # Temperaturer [K]
Vv = (18.01528/Vv_data[:, 3])*1000   # Væskevolum [mL]

# Har ulike temperaturverdier eksperimentelt. Hvordan slå disse sammen? 
# Alle datasettene er IKKE evaluert ved de samme verdiene av T. 



# 2b
# Henter fra oppgave 1:
T_c = 647.096         #[K]

def modell_L(T, a, b):
    return a*(abs(T-T_c))**(b)
    
def modell_Vg(T, a, b):
    return a*np.e**(b/T)

def modell_Vv(T, a, b, c, d, e):
    return a*T**4 + b*T**3 + c*T**2 + d*T + e


plt.title("L")
plt.plot(T1, L, label="Exp. L")
a_L = optimize.curve_fit(modell_L, T1, L)[0][0]
b_L = optimize.curve_fit(modell_L, T1, L)[0][1]
plt.plot(T1, modell_L(T1, a_L, b_L), label="Tilnærming") 
plt.legend()
plt.show()

plt.title("Vg")
plt.plot(T2, Vg, label="Exp. Vg")  # 1/x**2
a_Vg = optimize.curve_fit(modell_Vg, T2, Vg)[0][0]
b_Vg = optimize.curve_fit(modell_Vg, T2, Vg)[0][1]
plt.plot(T2, modell_Vg(T2, a_Vg, b_Vg), label="Tilnærming" )
plt.legend()
plt.show()

plt.title("Vv")
plt.plot(T3, Vv, label="Exp. Vv")  
a_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][0]
b_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][1]
c_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][2]
d_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][3]
e_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][4]
plt.plot(T3, modell_Vv(T3, a_Vv, b_Vv, c_Vv, d_Vv, e_Vv), label="Tilnærming") 
plt.legend()
plt.show()

#2c)

def clapeyrons(T):
    L = modell_L(T, a_L, b_L)
    Vg = modell_Vg(T, a_Vg, b_Vg)
    Vv =  modell_Vv(T, a_Vv, b_Vv, c_Vv, d_Vv, e_Vv)
    return L / (T * (Vg - Vv))

def simpson(T_start, T_slutt, n, f):
    
    '''
    Utfører Simpsons metode
    a = startverdi
    b = sluttverdi
    n = antall intervaller
    f = funksjonen
    '''
    
    h = (T_slutt - T_start) / n                                 #lengden på intervallene i x-aksen
    S = 0                                       
    for i in range(1, int(n/2)):
        x = T_start + (2*i-1)*h
        S = S + 4*f(x)
        x = T_start + 2*i*h
        S = S + 2*f(x)
    S = S + 4*f(T_slutt - h) + f(T_start) + f(T_slutt)
    S = h*S/3
    
    return S

# def simpson_med_np(T_start, T_slutt, n, f):
#     '''
#     Utfører Simpsons metode med numpy
#     a = startverdi
#     b = sluttverdi
#     n = antall intervaller, må være et partall
#     f = funksjonen
#     '''
#     h = (T_slutt - T_start)/n                                            #lengde på intervall
#     x_verdier = np.linspace(T_start, T_slutt, n+1)                       #(start, slutt, antall punkter)
#     x_verdier1 = x_verdier[1:-1:2]                           
#     x_verdier2 = x_verdier[2:-2:2]
#     f_x1 = 4*f(x_verdier1)
#     f_x2 = 2*f(x_verdier2)
#     S = (f(a) + f(b) + np.sum(f_x1) + np.sum(f_x2)) * (h/3)
#     return S

T_0 = np.arange(273.5, 646.5, 1)
T_1 = np.arange(274.5, 647.5, 1)

p = simpson(T_0, T_1, 100, clapeyrons)
plt.plot(((T_0+T_1)/2), p)