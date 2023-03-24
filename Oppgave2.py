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
T1 = L_data[:, 0] + 273               # Temperaturer [K]
L = L_data[:, 2]                     # Heat of vaporization [J/mol]

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

def modell1(T, a, b):
    return a*(abs(T-T_c))**(b)
    
def modell2(T, a, b, c):
    return a/np.log10(T) + b/T + c

def modell3(T, a, b, c, d):
    return a*T**3 + b*T**2 + c*T + d


plt.title("L")
plt.plot(T1, L, label="Exp. L")
a = optimize.curve_fit(modell1, T1, L)[0][0]
b = optimize.curve_fit(modell1, T1, L)[0][1]
plt.plot(T1, modell1(T1, a, b), label="Tilnærming") 
plt.legend()
plt.show()

plt.title("Vg")
plt.plot(T2, Vg, label="Exp. Vg")  # 1/x**2
a = optimize.curve_fit(modell2, T2, Vg)[0][0]
b = optimize.curve_fit(modell2, T2, Vg)[0][1]
c = optimize.curve_fit(modell2, T2, Vg)[0][2]
#d = optimize.curve_fit(modell2, T2, Vg)[0][3]
plt.plot(T2, modell2(T2, a, b, c), label="Tilnærming" )
plt.legend()
plt.show()

plt.title("Vv")
plt.plot(T3, Vv, label="Exp. Vv")  
a = optimize.curve_fit(modell3, T3, Vv)[0][0]
b = optimize.curve_fit(modell3, T3, Vv)[0][1]
c = optimize.curve_fit(modell3, T3, Vv)[0][2]
d = optimize.curve_fit(modell3, T3, Vv)[0][3]
plt.plot(T3, modell3(T3, a, b, c, d), label="Tilnærming") 
plt.legend()
plt.show()

