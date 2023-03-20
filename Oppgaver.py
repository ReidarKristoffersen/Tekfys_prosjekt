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