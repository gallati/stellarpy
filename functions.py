import pandas as pd
import numpy as np

# En todas las funciones las unidades son aquellas dadas por la Tabla 2,
# lo que lleva a cambiar expresiones como la de generación de energía (9)

######################################################################
######################### Constantes físicas #########################
######################################################################


# Definición de parámetros constantes

Mtot = 5.0
Rtot = 11.5
Ltot = 70.0
Tc = 2.0
X = 0.75
Y = 0.22
Z = 1 - X - Y
mu = 1/(2*X + 3*Y/4 + Z/2)
R = 8.31447 * 10**7


# Definición de la Tabla 1 como DataFrame (valores de T están en unidades del modelo)

epsilon_df = pd.DataFrame(data=[("pp", 0.40, 0.60, -6.84, 6.0),
                                ("pp", 0.60, 0.95, -6.04, 5.0), 
                                ("pp", 0.95, 1.20, -5.56, 4.5), 
                                ("pp", 1.20, 1.65, -5.02, 4.0), 
                                ("pp", 1.65, 2.40, -4.40, 3.5), 
                                ("CN", 1.20, 1.60, -22.2, 20.0),
                                ("CN", 1.60, 2.25, -19.8, 18.0), 
                                ("CN", 2.25, 2.75, -17.1, 16.0), 
                                ("CN", 2.75, 3.60, -15.6, 15.0), 
                                ("CN", 3.60, 5.00, -12.5, 13.0) ], 
                            columns=["Ciclo", "Tmin", "Tmax", "log10(epsilon1)", "nu"])


# Definición de parámetros no constantes

def rho(P, T):
    """
    Ecuación (6)
    """
    return (mu/R)*(P/T)

def epsilon(P, T):
    """
    Ecuación (9)
    """

    # Calculamos los valores de epsilon1, X1, X2 y nu según la temperatura
    epsilon1, X1, X2, nu, ciclo = calculo_tabla1(T)

    # Devolvemos el valor correspondiente junto con el ciclo que lo produce (valores de T están en unidades del modelo)
    return (epsilon1*X1*X2*rho(P, T)*(T/10)**nu, ciclo)

def calculo_tabla1(T):
    """
    Dada una temperatura T, devuelve los valores de epsilon1, X1, X2 y nu
    junto con el ciclo que genera la energía.
    """

    # Aplicamos un filtro al DataFrame de la Tabla 1
    filtro = (epsilon_df["Tmin"] <= T) & (T < epsilon_df["Tmax"])
    epsilon_values = epsilon_df[filtro]

    # Si el DataFrame está vacío, devolvemos 0.0
    if epsilon_values.empty:
        return (0.0, 0.0, 0.0, 0.0, "--")
    
    # Si el DataFrame no está vacío, devolvemos para el mayor epsilon1
    else:
        # Tomamos el índice de la fila para el que log10(epsilon1) es mayor y extraemos sus elementos
        index = epsilon_values["log10(epsilon1)"].idxmax()
        ciclo, _, _, log_epsilon1, nu = epsilon_values.loc[index].to_list()

        # Calculamos X1 o X2 según el ciclo
        if ciclo == "pp":
            X1, X2 = X, X
        elif ciclo == "CN":
            X1, X2 = X, Z/3

        # Devolvemos los valores de epsilon1, X1, X2, nu y el ciclo que lo genera
        return 10**log_epsilon1, X1, X2, nu, ciclo



######################################################################
### Ecuaciones fundamentales en las unidades del modelo (Tabla 3)  ###
######################################################################


# Caso radiativo

def dMdr_rad(r, P, T):
    """
    Ecuación (18)
    """
    Cm = 0.01523*mu
    return Cm * (P*r**2) / T

def dPdr_rad(r, P, T, M):
    """
    Ecuación (19)
    """
    Cp = 8.084*mu
    return - Cp * P*M / (T*r**2)

def dLdr_rad(r, P, T):
    """
    Ecuación (20)
    """
    epsilon1, X1, X2, nu, ciclo = calculo_tabla1(T)
    Cl = 0.01845*epsilon1*X1*X2*(10**nu)*mu**2
    return Cl * ((P*r)**2) * (T**(nu-2)), ciclo

def dTdr_rad(r, P, T, L):
    """
    Ecuación (21)
    """
    Ct = 0.01679*Z*(1+X)*mu*mu
    return -Ct * (L*(P**2)) / ((T**8.5)*(r**2))


# Caso convectivo

def dMdr_conv(r, T):
    """
    Ecuación (22)
    """
    Cm = 0.01523*mu
    return Cm * K*(T**1.5)*r**2

def dPdr_conv(r, T, M):
    """
    Ecuación (23)
    """
    Cp = 8.084*mu
    return -Cp * K*(T**1.5)*M / (r**2)

def dLdr_conv(r, T):
    """
    Ecuación (24)
    """
    Cl = 0.01845*epsilon1*X1*X2*(10**nu)*mu**2
    return Cl * (K**2)*(T**(3+nu))*r**2

def dTdr_conv(r, M):
    """
    Ecuación (25)
    """
    Ct = 3.234*mu
    return -Ct * M / (r**2)



######################################################################
############ Ecuaciones para encontrar valores iniciales #############
######################################################################


# Valores iniciales en la superficie

# Envoltura radiativa
def T_inicial_superficie(r):
    """
    Ecuación (35)
    """
    A1 = 1.9022*mu*Mtot
    return A1*(1/r - 1/Rtot)

# Envoltura convectiva
def P_inicial_superficie(r, T):
    """
    Ecuación (36)
    """
    A2 = 10.645*((Mtot/Ltot)/(mu*Z*(1+X)))**0.5
    return A2*T**4.25


# Valores iniciales en el centro

def M_inicial_centro(r, T):
    return 0.005077*mu*K*(T**1.5)*r**3

def L_inicial_centro(r):
    return 0.006150*epsilon1*X1*X2*(10**nu)*(mu**2)*(K**2)*(Tc**(3+nu))*(r**3)

def T_inicial_centro(r):
    return Tc - 0.008207*(mu**2)*K*(Tc**1.5)*(r**2)

def P_inicial_centro(r, T):
    return K*T**2.5


# Constante del polítropo

def politropo(P_dado, T_dado)
    return P_dado/(T_dado**2.5)



######################################################################
###################### Definición de algoritmos ######################
######################################################################


############################### Paso 2 ###############################

def paso2(modelo, derivadas, h, i):

    # Calculamos la estimación de P y T en la capa i+1

    P = modelo["P"][i]          # P en la capa i
    T = modelo["T"][i]          # T en la capa i
    fP = derivadas["fP"][i]     # fP en la capa i
    fT = derivadas["fT"][i]     # fT en la capa i
    AP1 = h * (derivadas["fP"][i] -  derivadas["fP"][i])                            # AP1 en la capa i
    AP2 = h * (derivadas["fP"][i] - 2*derivadas["fP"][i-1] + derivadas["fP"][i-2])  # AP2 en la capa i
    AT1 = h * (derivadas["fT"][i] -  derivadas["fT"][i])                            # AT1 en la capa i

    P_est = P + h*fP + AP1/2 + 5*AP2/12     # P estimado en la capa i+1
    T_est = T + h*fT + AT1/2                # T estimado en la capa i+1

    return P_est, T_est


############################### Paso 3 ###############################

def paso3(modelo, derivadas, P_dado, T_dado, M_dado, h, i):

    # Calculamos fM en la capa i+1 (P, T y M se dan en la capa i+1)

    r = modelo["r"][i] + h     # r en la capa i+1
    fM = dMdr_rad(r, P_dado, T_dado)     # fM en la capa i+1


    # Calculamos M en la capa i+1 teniendo

    AM1 = h * (fM - derivadas["fM"][i])       # AM1 en la capa i+1

    M_cal = M_dado + h*fM - AM1/2             # M calculado en la capa i+1

    return M_cal, fM


############################### Paso 4 ###############################

def paso4(modelo, derivadas, P_dado, T_dado, M_dado, h, i):

    # Calculamos fP en la capa i+1 (P, T y M se dan en la capa i+1)

    r = modelo["r"][i] + h     # r en la capa i+1
    fP = dPdr_rad(r, P_dado, T_dado, M_dado)  # fP en la capa i+1


    # Calculamos P en la capa i+1

    P = modelo["P"][i]                   # P en la capa i
    AP1 = h * (fP - derivadas["fP"][i])  # AP1 en la capa i+1

    P_cal = P + h*fP - AP1/2             # P calculado en la capa i+1

    return P_cal, fP


############################### Paso 5 ###############################

def paso5(P_cal, P_est):

    # Comprobamos si el error relativo es menor al error relativo máximo

    Error_relativo_maximo = 0.0001
    
    return abs(P_cal-P_est)/abs(P_cal) < Error_relativo_maximo


############################### Paso 6 ###############################

def paso6(modelo, derivadas, P_dado, T_dado, L_dado, h, i):

    # Calculamos fL en la capa i+1 (P, T y L se dan en la capa i+1)

    r = modelo["r"][i] + h     # r en la capa i+1
    fL, ciclo = dLdr_rad(r,P_dado,T_dado)       # fL en la capa i+1

    # Calculamos L en la capa i+1
    AL1 = h * (fL - derivadas["fL"][i])                            # AL1 en la capa i+1
    AL2 = h * (fL - 2*derivadas["fL"][i] + derivadas["fL"][i-1])   # AL2 en la capa i+1
    L_cal = L_dado + h*fL - AL1/2  - AL2/12                        # L calculado en la capa i+1

    return L_cal, fL, ciclo


############################### Paso 7 ###############################

def paso7(modelo, derivadas, P_dado, T_dado, L_dado, h, i):

    # Calculamos fT en la capa i+1 (P, T y L se dan en la capa i+1)

    r = modelo["r"][i] + h     # r en la capa i+1
    fT = dTdr_rad(r, P_dado, T_dado, L_dado)  # fT en la capa i+1


    # Calculamos T en la capa i+1

    T = modelo["T"][i]                   # T en la capa i
    AT1 = h * (fT - derivadas["fT"][i])  # AT1 en la capa i+1

    T_cal = T + h*fT - AT1/2            # T calculado en la capa i+1

    return T_cal, fT


############################### Paso 8 ###############################

def paso8(T_cal, T_est):

    # Comprobamos si el error relativo es menor al error relativo máximo

    Error_relativo_maximo = 0.0001
    
    return abs(T_cal-T_est)/abs(T_cal) < Error_relativo_maximo


############################### Paso 9 ###############################

def paso9(P_dado, T_dado, fP_dado, fT_dado):
    return T_dado*fP_dado/(P_dado*fT_dado)


############################### Paso 10 ###############################

def paso10(n1):
    return n1 <= 2.5


############################### Paso X ###############################

def pasoX(x_cal, x_est):

    # Comprobamos si el error relativo es menor al error relativo máximo

    Error_relativo_maximo = 0.0001
    
    return abs(x_cal-x_est)/abs(x_cal) < Error_relativo_maximo

