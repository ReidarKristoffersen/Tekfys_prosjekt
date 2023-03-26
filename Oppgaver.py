import numpy as np
import matplotlib.pyplot as plt
import scipy
import sympy 

# PARAMETERE

R = 8.314             #[J / mol K]
R = 8314.46261815324  #[L Pa/ mol K]
R = R*1000            #[mL Pa/ mol K]
R = R*10**(-5)        #[mL bar/ mol K]
T_c = 647.096         #[K]
p_c = 22.064 * 10e6   #[Pa] 
p_c = p_c*10**(-5)    #[bar] 
#n = 1                #[mol] 


# OPPGAVE 1a)

#Fra nettet (engineeringtoolbox)
a = 5.537                              #[bar * L^2 / mol^2]
a = a * 1000**2                        #[bar * mL^2 / mol^2]
b = 0.03049                            #[L / mol]
b = b * 1000                           #[mL / mol]
V_c = 3 * b


print(a, b, V_c)
#Eksperimentell V_c = 55.948mL.
#Her ser vi at volumet er 91.47 mL. Man kan bare få samsvar med maks to termodynamiske variable. Her: p_c og T_c.


# OPPGAVE 1b)

V_val = np.linspace(75, 300, 1000)     #[mL]

def vdW(T,V):
    """
    Uses Van der Waals tilstandsligning for å 

    Parameters
    ----------
    T : float
        Temperatur
    V : float
        Volum

    Returns
    -------
    p : float
        Trykk ved gitt temperatur og volum

    """
    p= ((R * T) / (V - b)) - (a / V**2)
    return p
vdW_tilst_likn = vdW(T_c, V_val)

plt.figure(figsize=(15,8))
plt.plot(V_val, vdW_tilst_likn)
plt.xlabel("V [ml]")
plt.ylabel("p [bar]")
plt.title("1 mol $H_2 O$, van der Waals isoterm, T = 647.096 K")
plt.grid()
plt.show() 

#Det som kan observeres her som er ulikt fra plottet i oppgaveteksten er at her så er den avtagende hele veien. I oppgaveteksten øker den først til kritisk punkt før den blir avtagende.
#Det vi også kan se er at trykket er høyere for alle verdier av temperaturen.







































































#1c)

c= 1
def ex_eq(T):
    """
    Eksempel ligningen vi vil løse, for å teste implimentering av Newtons metode

    Parameters
    ----------
    c : float
        Variabel

    Returns
    -------
        float
        Funksjonsverdi i gitt c-verdi

    """
    return (np.sinh(2*c/T))**2-1

def d_ex_eq(T):
    """
    Den deriverte av eksempel ligningen vi vil løse, for å teste implimentering av Newtons metode

    Parameters
    ----------
    c : float
        Variabel

    Returns
    -------
        float
        derivert funksjonsverdi i gitt c-verdi
    """
    return (-2)*c*np.sinh(4*c/T)/T**2


def newton_one_var(f, df, T_0, tol=10e-12):
    """
    Gjennomfører newton's metode på funksjonen tatt inn som paramter. Returnerer alle x_i i en følge

    Parameters
    ----------
    f : Function
        Ligning som det tilnærmes en løsning av
    df : Function
        Deriverte av f
    c_0 : float
        Startverdi for f
    tol : TYPE
        Toleranse for løsning

    Returns
    -------
    x : array
        Liste med alle x_i regnet, hvor siste vil være den endelige løsningen innnenfor toleransen.

    """
    max_it =100000
    i = 0
    x = np.array([T_0])
    while(np.abs(f(x[-1]))>=tol):
        x_i = x[-1] - (f(x[-1])/df(x[-1]))
        #print(f"f: {f(x[-1])}")
        #print(f"df: {df(x[-1])}")
        x = np.append(x, x_i)
        i+=1
        if(i>max_it):
            break
            raise Exception("Maximum iterations reached")
    print(f"The root was found to be at {x[-1]} after {i+1} iterations")
    return x

T_0 =1
T_arr = newton_one_var(ex_eq, d_ex_eq, T_0, 10e-6)
r = T_arr[-1]
print(T_arr)


T_0 = 3 #Vet nå at løsningen er ca. 285. Sjekker hvor fortere metoden konverger hvis vi velger c_0 = 300
T_arr_300 = newton_one_var(ex_eq, d_ex_eq, T_0, 10e-6)
print(r-2*c/np.log(1+np.sqrt(2))) #Sjekker at løsningen stemmer med ligningen gitt i oppgaveteksten


#1d)

e_i = np.abs(T_arr-r) #definer array med feil for hver iterasjon
plt.plot(e_i)
plt.show()
p_i = np.log(e_i[2:]/e_i[1:-1])/np.log(e_i[1:-1]/e_i[0:-2]) 
plt.plot(p_i)
plt.title("$p_i$, konvergerer mot q")
plt.xlabel("i")
plt.ylabel("$p_i$")
plt.show()

#Ser at p_i konvergerer mot 2. Dette er q-verdien samsvarer med at Newtons metode konvergerer kvadratisk


#1e)
"""
SYMPY GREIER, IKKE I BRUK PER NÅ

def func1():
    return R_sym*T_sym*(1 / (V_g_sym - b_sym) - 1 / (V_v_sym - b_sym)) - a_sym*(1 / V_g_sym**2 - 1 / V_v_sym**2)  # 11

def func2():
    return R_sym*T_sym*(sympy.log((V_g_sym - b_sym) / (V_v_sym - b_sym)) / (V_g_sym - V_v_sym) ) - a_sym*(1 / (V_g_sym*V_v_sym) - 1 / (V_g_sym**2)) - R_sym*T_sym/(V_g_sym-b_sym) # 12

a_sym, b_sym, R_sym, T_sym, V_g_sym, V_v_sym = sympy.symbols('a b R T V_g V_v') 

f1 = func1()
f2 = func2()

df1_dg = sympy.diff(f1, V_g_sym)
df1_dv = sympy.diff(f1, V_v_sym)

df2_dg = sympy.diff(f2, V_g_sym)
df2_dv = sympy.diff(f2, V_v_sym)

Jacobi_uneval = sympy.Matrix([[df1_dg, df1_dv], [df2_dg, df2_dv]])

def Jacobi(J, x0, y0,T):
    
   
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
    
    return np.array(J.subs([(a_sym, a), (b_sym, b), (R_sym, R), (T_sym, T), (V_g_sym, x0), (V_v_sym, y0)]), dtype=np.float64)
"""
def f1(V_g, V_v, T):
    return R*T*(1 / (V_g - b) - 1 / (V_v - b)) - a*(1 / V_g**2 - 1 / V_v**2)  # 11

def f2(V_g, V_v, T):
    return R*T*(np.log((V_g - b) / (V_v - b)) / (V_g - V_v) ) - a*(1 / (V_g*V_v) - 1 / (V_g**2)) - R*T/(V_g-b) # 12

def f_tot(V_arr, T=273):
    return np.array([f1(V_arr[0], V_arr[1], T), f2(V_arr[0], V_arr[1], T)]).T


def Jacobi_numpy(V_g,V_v,T):
    """
    Returrnerer jacobimatrisa ved gitt V_g og V_v

    Parameters
    ----------
    V_g : float
        Gassvolum
    V_v : float
        væskevolum
    T : float
        Temperature.

    Returns
    -------
    J : Array 2x2
        Jacobi

    """
    df1_dg = -R*T/(V_g - b)**2 + 2*a/V_g**3
    df1_dv=R*T/(V_v - b)**2 - 2*a/V_v**3
    df2_dg = R*T/(V_g - b)**2 + R*T/((V_g - V_v)*(V_g - b)) - R*T*np.log((V_g - b)/(V_v - b))/(V_g - V_v)**2 - a*(-1/(V_g**2*V_v) + 2/V_g**3)
    df2_dv = -R*T/((V_g - V_v)*(V_v - b)) + R*T*np.log((V_g - b)/(V_v - b))/(V_g - V_v)**2 + a/(V_g*V_v**2)
    J = np.array([[df1_dg, df1_dv], [df2_dg, df2_dv]])
    return J


def newton_two_var(f, V_arr, T, tol = 10e-12):
    """
    

    Parameters
    ----------
    f : TYPE
        DESCRIPTION.
    Jacobi_uneval : TYPE
        DESCRIPTION.
    V_array : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    tol: 

    Returns
    -------
    None.

    """
    max_it = 100
    i = 0
    while np.linalg.norm(f(V_arr, T)) >= tol:
        #print(V_array)
        J = Jacobi_numpy(V_arr[0], V_arr[1], T)
        #print("\n")
        #print(J)
        #print("\n")
        f_ans = f(V_arr, T)
        #print(f"V:{V_arr}")
        #print(f"f:{f_ans}")
        #print(J.shape, f_ans.T.shape)
        
        #print(f_ans.reshape(2,1))
        #delta = np.linalg.solve(J, -f_ans.reshape(2,1))
        delta = np.linalg.solve(J, -f_ans) #ny
        V_arr += delta #[V_g, V_v]
        
        i += 1
        if(i>max_it):
            print("Maximum iterations reached")
            break
    return V_arr
            
V_array = np.array([12.6e+3, 35.7])
#fSolve = scipy.optimize.fsolve(f_tot, x0 = V_array)

print(newton_two_var(f_tot, V_array, 274))

T_arr = np.linspace(274,647,747)
V_of_T = np.zeros((len(T_arr),2))

for i in range(len(T_arr)):
    newton_two_var(f_tot, V_array, T_arr[i])
    V_of_T[i]= V_array

#BRUKER LOGARITMISK PLOTT FOR Å VISE AT DE MØTES, USIKKER PÅ HVA VI SKAL SI OM ENHETER
plt.semilogy(T_arr,V_of_T, label = ["$V_g$","$V_v$"]) 
plt.title("Logaritmisk plot av $V_g$ og $V_v$")
plt.xlabel("T [K]")
plt.vlines(T_c,20,10**4, color="gray", label = "$T_c$")
plt.legend()
plt.show()

#1f) SKJØNNER IKKE



# Leser datafilen
Vg_data = np.loadtxt("Vg_verdier.txt")
T_exp = Vg_data[:, 0]                    # Temperaturer [K]
p_g_exp = (Vg_data[:, 2])                  # Trykk [bar]

#Trippelpinkt ved T = 273.16, kritisk ved T = 647.096 (All data mellom dette)
p_an_g = vdW(T_arr, V_of_T.T[0])
p_an_v = vdW(T_arr,V_of_T.T[1])

plt.semilogy(T_arr, p_an_g, label = "Numerisk gasss")
plt.semilogy(T_arr, p_an_v, label = "Numerisk væske")
plt.semilogy(T_exp,p_g_exp, label = "Eksperimentell")
plt.xlabel("T [K]")
plt.ylabel("log of p [log(bar)]")
plt.legend()
plt.show()
#SER AT AVVIKET TYDELIG MINKER NÅR TEMPERATUREN ØKER 


#1g)

rel_index = 602
T_rel  = T_arr[rel_index]
V_rel = V_of_T[rel_index]
print(V_rel)

rel_V_arr = np.linspace(V_rel[1]-3,V_rel[0]+100,1000)

p_rel = vdW(T_rel, rel_V_arr)

p_V_g = vdW(T_rel, V_rel[0])

plt.plot(rel_V_arr, p_rel, label = "Damptrykk")
plt.scatter(V_rel[1], p_V_g, color = "red", label = "$V_v$")
plt.scatter(V_rel[0], p_V_g, color = "orange", label = "$V_g$")
plt.hlines(p_V_g, V_rel[1], V_rel[0], color = "blue", linestyles="--", label = "Faseovergang")
plt.title(f"Damptrykk ved T = {T_rel}K")
plt.xlabel("T [K]")
plt.ylabel("p [bar]")
plt.legend()
plt.show()

#1e)

p_V_v = vdW(T_c, V_of_T.T[1])
p_V_g = vdW(T_c, V_of_T.T[0])


#plt.plot(V_of_T.T[1][400:],p_V_v)
#plt.plot(V_of_T.T[0][:100],p_V_g)
plt.plot(V_of_T.T[1],p_V_v)
plt.plot(V_of_T.T[0] ,p_V_g)

