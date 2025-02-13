from model import Model
from minimum import error_table
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Creamos el objeto de modelo
# modelo = Model(Mtot=5.0, Rtot=11.0570, Ltot=75.9213, Tc=1.9554)

# # Método para la representación gráfica
# modelo.visualize(merge=True)

# # Tabla de resultados
# print(modelo)

# Tabla resumen
# error_table(1, Mtot=5.0, Rtot=11.0570, Ltot=75.9213, Tc=1.9554, tol=50)
# Para tol=0.01 el tiempo de ejecución es 1 min y 45 seg
# Para tol=1 el tiempo de ejecución es 1 min y 30 seg
# Para tol=50 el tiempo de ejecución es 1 min y 30 seg

# Encontrando el mínimo
# def f(x):
#     return Model(Mtot=5.0, Rtot=x[0], Ltot=x[1], Tc=x[2]).error()

# result = optimize.minimize(f, x0=[10.93, 73.52, 1.95])
# print(f"Minimum at \n   Rtot = {result.x[0]:.4f}\n   Ltot = {result.x[1]:.4f}\n   Tc = {result.x[2]:.4f}")
# print(f"Error: {result.fun:.4f} %")


####### Representación en el diagrama rho-T
# modelo = Model() 

# R = 8.31447e7
# a = 7.56578e-15
# mu = modelo.mu
# mu_e = 2/(1+modelo.X)

# K0 = R/mu
# K1 = (1.0036e13)/(mu_e**(5/3))
# K2 = (1.2435e15)/(mu_e**(4/3))
# K3 = a/3

# logT = np.array([6, 10], dtype=float)
# logrho_I_II = 1.5*np.log10(K0/K1) + 1.5*logT
# logrho_I_III = 3*np.log10(K0/K2) + 3*logT
# logrho_II_III = 3*np.log10(K2/K1)*np.ones((2,))
# logrho_I_IV = np.log10(K3/(10*K0)) + 3*logT

# plt.plot(logT, logrho_I_II)
# plt.plot(logT, logrho_I_III)
# plt.plot(logT, logrho_II_III)
# plt.plot(logT, logrho_I_IV)
# plt.xlim([6,10])
# plt.ylim([-1,11])

# # Traza de la estrella
# plt.plot(np.log10(modelo.get("T")*1e7), np.log10(modelo.model["rho"]))

# # Zona de fusión de hidrógeno (diapositiva 4, tema 8) 
# epsilon_min = 1e20*modelo.X**2
# epsilon0 = 10**(-17.1)
# n = 4

# logrho4 = np.log10(epsilon_min/epsilon0)-20*np.linspace(6,7.5)
# plt.plot(np.linspace(6,7.5), logrho4)
#
# plt.show()


# Minimizar el error según la temperatura

# def f(x):
#     return Model(Tc=x[0]).error()

# result = optimize.minimize(fun=f, x0=[1.5], tol=50)

# print(f"T = {result.x[0]}. Error = {result.fun}. Converged = {result.success}")






# dataset = pd.read_csv("./6 class csv.csv.xls")

# O = dataset[dataset["Spectral Class"] == 'O']
# B = dataset[dataset["Spectral Class"] == 'B']
# A = dataset[dataset["Spectral Class"] == 'A']
# F = dataset[dataset["Spectral Class"] == 'F']
# G = dataset[dataset["Spectral Class"] == 'G']
# K = dataset[dataset["Spectral Class"] == 'K']
# M = dataset[dataset["Spectral Class"] == 'M']

# plt.style.use('ggplot')

# plt.grid()
# plt.scatter(np.log(O["Temperature (K)"]), np.log(O["Luminosity(L/Lo)"]),10, label = 'O')
# plt.scatter(np.log(B["Temperature (K)"]), np.log(B["Luminosity(L/Lo)"]),8, label = 'B')
# plt.scatter(np.log(A["Temperature (K)"]), np.log(A["Luminosity(L/Lo)"]),6, label = 'A')
# plt.scatter(np.log(F["Temperature (K)"]), np.log(F["Luminosity(L/Lo)"]),4, label = 'F')
# plt.scatter(np.log(G["Temperature (K)"]), np.log(G["Luminosity(L/Lo)"]),3, label = 'G')
# plt.scatter(np.log(K["Temperature (K)"]), np.log(K["Luminosity(L/Lo)"]),2, label = 'K')
# plt.scatter(np.log(M["Temperature (K)"]), np.log(M["Luminosity(L/Lo)"]),1, label = 'M')
# plt.gca().invert_xaxis()
# plt.title("Hertzsprung-Russell Diagram")
# plt.ylabel("log Luminosity")
# plt.xlabel("log T")
# plt.legend()

# L = np.log(73.52 / 3.842)
# T = np.log(Model(Rtot=10.93, Ltot=73.52, Tc=1.95).get("T").iloc[0]*1e7)
# plt.scatter(T, L, label = 'model')

# plt.show()


# ESTÁ EN MAGNITUD ABSOLUTA, HABRÍA QUE PASARLO A LUMINOSIDADES SOLARES

data = pd.read_csv("data.csv")

# format the points on the plot
transparency = 0.2
size = 1

# draws a scatter plot
plt.scatter(np.log10(data.temp), data.absmag, s=size, edgecolors='none', alpha=transparency)
# plt.xlim(np.log10(15000),np.log10(2000))
# plt.ylim(20,-15)
plt.title("H-R Diagram (log)")
plt.ylabel("Absolute Magnitude")
plt.xlabel("Log T (log K)")

L = np.log(73.52 / 3.842) # Luminosidad en luminosidades solares
M = 4.83 - 2.5*np.log(L)  # Magnitud absoluta comparando con la solar
T = np.log(Model(Rtot=10.93, Ltot=73.52, Tc=1.95).get("T").iloc[0]*1e7)
print(Model(Rtot=10.93, Ltot=73.52, Tc=1.95).get("T"))
plt.scatter(T, M, label = 'model')

plt.show()