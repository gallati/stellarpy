import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from modelo import Modelo

Rtot = 11.5
Ltot = 70.0
dR = 0.5  # Incremento en el radio
dL = 5.0  # Incremento en la luminosidad
n = 3     # Tamaño del máximo incremento [Ltot - n*dL, ..., Ltot + 0*dL, ..., Ltot + n*dL] 
d = 2*n+1 # Dimensiones de la tabla 

# Inicializamos la matriz donde se almacenan los errores
matriz_error = np.zeros((d, d), dtype=float)


# Crear la figura y el eje
fig, ax = plt.subplots()


for i in range(d):
    for j in range(d):
        matriz_error[i, j] = Modelo(Rtot=Rtot + (j-n)*dR, Ltot=Ltot+(i-n)*dL).error()
        ax.text(j, i, f'{matriz_error[i, j]:.2f} %', ha='center', va='center', color='white')




cax = ax.imshow(matriz_error, cmap='coolwarm', vmax=100.0)

ax.set_xticks(range(d))
ax.set_xticklabels([f"R$_{{Tot}}$ + {j-n}·$\\delta$R" for j in range(d)])

ax.set_yticks(range(d))
ax.set_yticklabels([f"L$_{{Tot}}$ + {i-n}·$\\delta$L" for i in range(d)])


plt.colorbar(cax)
plt.show()





# plt.plot()

# def minimo_error_modelo(x):
#     return Modelo(Mtot=5.0, Rtot=x[0], Ltot=x[1], Tc=1.95).error()

# initial_guess = [11.5, 70.0]

# result = minimize(minimo_error_modelo, initial_guess)

# print(result)

# print(Modelo(Mtot=5.0, Rtot=10.93, Ltot=73.52, Tc=1.95).error())