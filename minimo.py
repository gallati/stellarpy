import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from modelo import Modelo

def grafica_errores(n, Mtot=5.0, Rtot=11.5, Ltot=70.0, Tc=2.0, dR=0.5, dL=5.0):
    """
    - Descripción:
    Genera una gráfica que muestra el error relativo total para variaciones de L y R.

    - Parámetros:
        n: Tamaño del máximo incremento. 
            Se calculan incrementos como [Xtot - n*dX, ..., Xtot + 0*dX, ..., Xtot + n*dX]
        Mtot: Masa total del modelo (default 5.0)
        Rtot: Radio total del modelo (default 11.5)
        Ltot: Luminosidad total del modelo (default 70.0)
        Tc: Temperatura central del modelo (default 2.0)
        dR: Variación del radio (default 0.5)
        dL: Variación de la luminosidad (default 5.0)
    """

    # Dimensiones de la tabla e inicialización de la matriz de errores
    d = 2*n+1
    matriz_error = np.zeros((d, d), dtype=float)

    # Creamos la figura
    fig, ax = plt.subplots()

    # Calculamos la matriz de errores
    for i in range(d):
        for j in range(d):
            matriz_error[i, j] = Modelo(Mtot=Mtot, Rtot=Rtot + (j-n)*dR, Ltot=Ltot+(i-n)*dL, Tc=Tc).error()
            ax.text(j, i, f'{matriz_error[i, j]:.2f} %', ha='center', va='center')

    # Generamos el imshow y lo representamos
    imshow = ax.imshow(matriz_error, cmap='coolwarm', vmax=100.0)
    ax.set_title(f"M$_{{Tot}}$={Mtot}   R$_{{Tot}}$={Rtot}   L$_{{Tot}}$={Ltot}   T$_{{c}}$={Tc}   $\\delta$R={dR}   $\\delta$L={dL}")
    ax.set_xticks(range(d))
    ax.set_xticklabels([f"{j-n}·$\\delta$R" for j in range(d)])
    ax.set_yticks(range(d))
    ax.set_yticklabels([f"{i-n}·$\\delta$L" for i in range(d)])

    plt.colorbar(imshow)
    plt.show()


grafica_errores(1, Mtot=5.0, Rtot=10.93, Ltot=73.57, Tc=1.95, dR=0.1, dL=1.0)

# plt.plot()

# def minimo_error_modelo(x):
#     return Modelo(Mtot=5.0, Rtot=x[0], Ltot=x[1], Tc=1.95).error()

# initial_guess = [11.5, 70.0]

# result = minimize(minimo_error_modelo, initial_guess)

# print(result)

# print(Modelo(Mtot=5.0, Rtot=10.93, Ltot=73.52, Tc=1.95).error())

