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

file2 = "Vv_verdier.txt"

# Leser datafilen
Vv_data = np.loadtxt(file2)
T2 = Vv_data[:, 0]                    # Temperaturer [K]
Vv = (18.01528/Vv_data[:, 3])*1000   # Væskevolum [mL]

file3 = "Vg_verdier.txt"

# Leser datafilen
Vg_data = np.loadtxt(file3)
T3 = Vg_data[:, 0]                    # Temperaturer [K]
Vg = (18.01528/Vg_data[:, 3])*1000   # Gassvolum [mL]

# Har ulike temperaturverdier eksperimentelt. Hvordan slå disse sammen? 
# Alle datasettene er IKKE evaluert ved de samme verdiene av T. 



# 2b

def modell1(T):
    return 
    
def modell2(T, a):
    return a*T**2

plt.title("L")
plt.plot(T1, L)
plt.show()

def modell2(T, a):
    return a*T**2

plt.title("Vv")
plt.plot(T2, Vv, label = "Vv")  
plt.legend()
plt.plot(T2, modell(T2, optimize.curve_fit(modell, T2, Vv)[0])) 
print(optimize.curve_fit(modell, T2, Vv)[0])
plt.show()

plt.title("Vg")
plt.plot(T3, Vg)  # 1/x**2
plt.show()


