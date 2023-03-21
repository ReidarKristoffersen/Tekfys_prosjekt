import numpy as np
import matplotlib.pyplot as plt

# PARAMETERE

R = 8.314             #[J / mol K]
T_c = 647.096         #[K]
p_c = 22.064 * 10e6   #[Pa] 
#n = 1                 #[mol] 


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
vdW_tilst_likn = ((R * T_c) * 10 / (V_val - b)) - (a / V_val**2)

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

T_0 = 1
T_arr = newton_one_var(ex_eq, d_ex_eq, T_0, 10e-6)
r = T_arr[-1]
print(T_arr)


T_0 = 3 # Må velge en hensiktsmessig verdi. Kan ikke være for høy, da metoden ikke vil konvergere. 
T_arr_300 = newton_one_var(ex_eq, d_ex_eq, T_0, 10e-6)
print(r-2*c/np.log(1+np.sqrt(2))) # Sjekker at løsningen stemmer med ligningen gitt i oppgaveteksten


#1d)

e_i = np.abs(T_arr-r) # definer array med feil for hver iterasjon
plt.plot(e_i)
plt.show()
p_i = np.log(e_i[2:]/e_i[1:-1])/np.log(e_i[1:-1]/e_i[0:-2]) # Tar ikke med første element ettersom vi her fikk en rar verdi i plottet?????
plt.plot(p_i)
plt.title("$p_i$, konvergerer mot q")
plt.xlabel("i")
plt.ylabel("$p_i$")
plt.show()

#Ser at p_i konvergerer mot 2. Dette er q-verdien samsvarer med at Newtons metode konvergerer kvadratisk
