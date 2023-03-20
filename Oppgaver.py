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











































































#1c)

def ex_eq(c):
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
    return (np.sinh(2*c/T_c))**2-1

def d_ex_eq(c):
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
    return 2*np.sinh(4*c/T_c)/T_c


def newton_one_var(f, df, c_0, tol):
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
    x = np.array([c_0])
    while(np.abs(f(x[-1]))>=tol):
        x_i = x[-1] - (f(x[-1])/df(x[-1]))
        x = np.append(x, x_i)
        i+=1
        if(i>max_it):
            raise Exception("Maximum iterations reached")
    print(f"The root was found to be at {x[-1]} after {i} iterations")
    return x

c_0 = 1 #Starter med c = 1 som gitt i oppgavetekst
c_arr = newton_one_var(ex_eq, d_ex_eq, c_0, 10e-12)
r = c_arr[-1]
print(r)


c_0 = 300 #Vet nå at løsningen er ca. 285. Sjekker hvor fortere metoden konverger hvis vi velger c_0 = 300
c_arr_300 = newton_one_var(ex_eq, d_ex_eq, c_0, 10e-12)
print(T_c-2*r/np.log(1+np.sqrt(2))) #Sjekker at løsningen stemmer med ligningen gitt i oppgaveteksten


#1d)

e_i = np.abs(c_arr-r) #definer array med feil foor hver iterasjon
plt.plot(e_i[1:])
plt.show()
p_i = np.log(e_i[3:]/e_i[2:-1])/np.log(e_i[2:-1]/e_i[1:-2]) #Tar ikke med første element ettersom vi her fikk en rar verdi i plottet?????
plt.plot(p_i)
plt.title("$p_i$, konvergerer mot q")
plt.xlabel("i")
plt.ylabel("$p_i$")
plt.show()

#Ser at p_i konvergerer mot 2. Dette er q-verdien samsvarer med at Newtons metode konvergerer kvadratisk