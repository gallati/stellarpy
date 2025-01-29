from model import Model
from minimum import error_table
from scipy import optimize

# Creamos el objeto de modelo
modelo = Model(Mtot=5.0, Rtot=11.0570, Ltot=75.9213, Tc=1.9554)

# Método para la representación gráfica
modelo.visualize(merge=True)

# Tabla de resultados
print(modelo)

# Tabla resumen
error_table(2, Mtot=5.0, Rtot=11.0570, Ltot=75.9213, Tc=1.9554)

# Encontrando el mínimo
# def f(x):
#     return Model(Mtot=5.0, Rtot=x[0], Ltot=x[1], Tc=x[2]).error()

# result = optimize.minimize(f, x0=[10.93, 73.52, 1.95])
# print(f"Minimum at \n   Rtot = {result.x[0]:.4f}\n   Ltot = {result.x[1]:.4f}\n   Tc = {result.x[2]:.4f}")
# print(f"Error: {result.fun:.4f} %")