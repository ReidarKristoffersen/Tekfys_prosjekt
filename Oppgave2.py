"""
Created on Fri Mar 24 16:15:59 2023

@author: hanna
"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy import optimize, interpolate

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
T2 = Vg_data[:, 0]                       # Temperaturer [K]
Vg = (18.01528/Vg_data[:, 3])*1000       # Gassvolum [mL]


file3 = "Vv_verdier.txt"

# Leser datafilen
Vv_data = np.loadtxt(file3)
T3 = Vv_data[:, 0]                       # Temperaturer [K]
Vv = (18.01528/Vv_data[:, 3])*1000       # Væskevolum [mL]

# Har ulike temperaturverdier eksperimentelt. Hvordan slå disse sammen? 
# Alle datasettene er IKKE evaluert ved de samme verdiene av T. 

# Molar masse (H2O) = 18.01528
# Omgjøring fra massetetthet til volum (antar 1 mol stoff)



# 2b
# Henter fra oppgave 1:
T_c = 647.096         #[K]

def modell_L(T, a, b):
    """
    

    Parameters
    ----------
    T : TYPE
        DESCRIPTION.
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return a*(abs(T-T_c))**(b)
    
def modell_Vg(T, a, b):
    return a*np.e**(b/T)

def modell_Vv(T, a, b, c, d, e):
    return a*T**4 + b*T**3 + c*T**2 + d*T + e


plt.title("Latent varme, L")
plt.plot(T1, L, label="Exp. L")
plt.xlabel("Varme [bar*mL / mol")
plt.ylabel("Volum [mL]") #[J/mol]?
a_L = optimize.curve_fit(modell_L, T1, L)[0][0]
b_L = optimize.curve_fit(modell_L, T1, L)[0][1]
plt.plot(T1, modell_L(T1, a_L, b_L), label="Tilnærming") 
plt.legend()
plt.show()

plt.title("Gassvolum, Vg")
plt.plot(T2, Vg, label="Exp. Vg")
plt.xlabel("Temperatur [K]")
plt.ylabel("Volum [mL]")
a_Vg = optimize.curve_fit(modell_Vg, T2, Vg)[0][0]
b_Vg = optimize.curve_fit(modell_Vg, T2, Vg)[0][1]
plt.plot(T2, modell_Vg(T2, a_Vg, b_Vg), label="Tilnærming" )
plt.legend()
plt.show()

plt.title("Væskevolum, Vv")
plt.plot(T3, Vv, label="Exp. Vv")  
plt.xlabel("Temperatur [K]")
plt.ylabel("Volum [mL]")
a_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][0]
b_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][1]
c_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][2]
d_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][3]
e_Vv = optimize.curve_fit(modell_Vv, T3, Vv)[0][4]
plt.plot(T3, modell_Vv(T3, a_Vv, b_Vv, c_Vv, d_Vv, e_Vv), label="Tilnærming") 
plt.legend()
plt.show()

#2c)

def clapeyrons_curve(T):
    """
    

    Parameters
    ----------
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    # likning (13), bruker modellfunksjonen
    
    L = modell_L(T, a_L, b_L)
    Vg = modell_Vg(T, a_Vg, b_Vg)
    Vv =  modell_Vv(T, a_Vv, b_Vv, c_Vv, d_Vv, e_Vv)
    return L / (T * (Vg - Vv))

def simpson(T_start, T_slutt, n, f):
    """
    

    Parameters
    ----------
    T_start : TYPE
        DESCRIPTION.
    T_slutt : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.
    f : TYPE
        DESCRIPTION.

    Returns
    -------
    S : TYPE
        DESCRIPTION.

    """
    h = (T_slutt - T_start) / n                 #lengden på intervallene i x-aksen
    S = 0                                       
    for i in range(1, int(n/2)):
        x = T_start + (2*i-1)*h
        S = S + 4*f(x)
        x = T_start + 2*i*h
        S = S + 2*f(x)
    S = S + 4*f(T_slutt - h) + f(T_start) + f(T_slutt)
    S = h*S/3
    
    return S

# for å sjekke at Simpson fungerer
def f1(x):
    return x**2
print(simpson(1, 2, 8, f1))
# eksakt svar er 7/3 = 2.333

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


T_0 = 274
T_1 = np.arange(275, 647, 1)

# Leser datafilen
Vg_data = np.loadtxt("Vg_verdier.txt")
T_exp = Vg_data[:, 0]                      # Temperaturer [K]
p_g_exp = (Vg_data[:, 2])                  # Trykk [bar]

p_curve = simpson(T_0, T_1, 4, clapeyrons_curve)       #er denne i bar?
plt.plot(T_1, p_curve, label = "Simpson - curve fit")
plt.plot(T_exp,p_g_exp, label = "Eksperimentell")
plt.xlabel("Temperatur [K]")
plt.ylabel("Trykk [bar]")
plt.legend()
plt.show()

#2d)

def CubicSpline_interpol():
    """
    

    Returns
    -------
    None.

    """
    L_cub = interpolate.CubicSpline(T1, L)
    Vg_cub = interpolate.CubicSpline(T2, Vg)
    Vv_cub = interpolate.CubicSpline(T3, Vv)
    return L_cub, Vg_cub, Vv_cub

L_cubic, Vg_cubic, Vv_cubic = CubicSpline_interpol()

plt.plot(T1, L_cubic(T1), label = "L cubic", color = "green")
plt.plot(T1, modell_L(T1, a_L, b_L), label = "L curve fit", linestyle="--")
plt.plot(T1, L, label="Exp. L", linestyle = ":", color = "red")
plt.legend()
plt.show()

plt.plot(T2, Vg_cubic(T2), label = "Vg cubic", color = "green")
plt.plot(T2, modell_Vg(T2, a_Vg, b_Vg), label = "Vg curve fit", linestyle="--")
plt.plot(T2, Vg, label="Exp. Vg", linestyle = ":", color = "red")
plt.legend()
plt.show()

plt.plot(T3, Vv_cubic(T3), label = "Vv cubic", color = "green")
plt.plot(T3, modell_Vv(T3, a_Vv, b_Vv, c_Vv, d_Vv, e_Vv), label = "Vv curve fit", linestyle="--")
plt.plot(T3, Vv, label="Exp. Vv", linestyle = ":", color = "red")  
plt.legend()
plt.show()    

#2e)

def clapeyrons_cubic(T):
    return L_cubic(T) / (T * (Vg_cubic(T) - Vv_cubic(T)))

p_cubic = simpson(T_0, T_1, 4, clapeyrons_cubic)       #er denne i bar?
plt.plot(T_1, p_cubic, label = "Simpson - cubic")
plt.plot(T_exp, p_g_exp, label = "Eksperimentell")
plt.xlabel("Temperatur [K]")
plt.ylabel("Trykk [bar]")
plt.legend()
plt.show()

print(T_exp[0])
print(p_g_exp[0])

#2f)

p_an_g = vdW(T_arr, V_of_T.T[0])                  #vdW
plt.plot(T_arr, p_an_g, label = "van der Waals")
plt.plot(T_1, p_curve, label = "Simpson - curve fit")
plt.plot(T_1, p_cubic, label = "Simpson - cubic")

#konstant L som vi må finne
def analytisk(L, p_0, T_0, T): 
    return p_0 * np.e((L / R) * ((1/T_0) - (1/T)))

plt.plot(T, analytisk(L, p_0, T_0, T), label = "Analytisk")

plt.plot(T_exp, p_g_exp, label = "Eksperimentell")
plt.legend()
plt.show()
    
