import numpy as np
import pandas as pd
from functions import *

# Configuración de los DataFrame

pd.set_option('colheader_justify', 'center', "display.precision", 7, 'display.max_rows', 100)
modelo = pd.DataFrame(data = [], columns=["E", "fase", "r", "P", "T", "L", "M", "n+1"])
derivadas = pd.DataFrame(data = [], columns=["fP", "fT", "fL", "fM"])


######################################################################
######################## Primeras tres capas #########################
######################################################################

Rini = 0.9*Rtot
h = -Rini/100
r = Rini

for i in range(3):

    # Calculamos y almacenamos los valores del modelo
    T = T_inicial_superficie(r)
    P = P_inicial_superficie(r, T)
    M = Mtot
    L = Ltot
    modelo.loc[len(modelo)] = {"E":"--", "fase":"INICIO", "r":r, "P":P, "T":T, "L":L, "M":M, "n+1":"-"}

    # Calculamos y almacenamos los valores de las f_i (derivadas)
    fT = dTdr_rad(r, P, T, L)
    fP = dPdr_rad(r, P, T, M)
    fM = 0.0
    fL = 0.0
    derivadas.loc[len(derivadas)] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

    # Aumentamos el valor del radio
    r += h


######################################################################
############################ Fase A.1.1. #############################
######################################################################

fase = "A.1.1."         # Fase actual (A.1.1.)

loop1 = True

while loop1:

    # Aplicamos el paso 2
    P_est, T_est = paso2(modelo, derivadas, h, i)

    loop2 = True

    while loop2:

        loop3 = True

        while loop3:
            
            # Aplicamos el paso 4
            P_cal, fP = paso4(modelo, derivadas, P_est, T_est, Mtot, h, i)

            # Aplicamos el paso 5
            if pasoX(P_cal, P_est):
                loop3 = False
            else:
                P_est = P_cal

        # Aplicamos el paso 7
        T_cal, fT = paso7(modelo, derivadas, P_cal, T_est, Ltot, h, i)

        # Aplicamos el paso 8
        if pasoX(T_cal, T_est):
            loop2 = False
        else:
            T_est = T_cal

    # Aplicamos el paso 3
    M_cal, fM = paso3(modelo, derivadas, P_cal, T_cal, Mtot, h, i)

    # Comparamos la masa calculada con la masa total
    if not pasoX(M_cal, Mtot):
        loop1 = False
    else:
        # Añadimos las variables calculadas al modelo (asumimos que M y L permanecen constantes)
        r = modelo["r"][i] + h
        modelo.loc[len(modelo)] = {"E":"--", "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":Ltot, "M":Mtot, "n+1":"-"}
        derivadas.loc[len(derivadas)] = {"fP":fP, "fT":fT, "fL":0.0, "fM":fM}

        # Calculamos la siguiente capa
        i += 1


######################################################################
############################ Fase A.1.2. #############################
######################################################################

fase = "A.1.2."         # Fase actual ("A.1.2.")

loop1 = True

while loop1:

    # Aplicamos el paso 2
    P_est, T_est = paso2(modelo, derivadas, h, i)

    loop2 = True

    while loop2:

        loop3 = True

        while loop3:
            
            # Aplicamos el paso 3
            M_cal, fM = paso3(modelo, derivadas, P_est, T_est, modelo["M"][i], h, i)

            # Aplicamos el paso 4
            P_cal, fP = paso4(modelo, derivadas, P_est, T_est, M_cal, h, i)

            # Aplicamos el paso 5
            if pasoX(P_cal, P_est):
                loop3 = False
            else:
                P_est = P_cal

        # Aplicamos el paso 7
        T_cal, fT = paso7(modelo, derivadas, P_cal, T_est, Ltot, h, i)

        # Aplicamos el paso 8
        if pasoX(T_cal, T_est):
            loop2 = False
        else:
            T_est = T_cal

    # Aplicamos el paso 6
    L_cal, fL, ciclo = paso6(modelo, derivadas, P_cal, T_cal, Ltot, h, i)

    # Comparamos la luminosidad calculada con la luminosidad total
    if not pasoX(L_cal, Ltot):
        loop1 = False
    else:
        # Añadimos las variables calculadas al modelo (asumimos que L permanece constante)
        r = modelo["r"][i] + h
        modelo.loc[len(modelo)] = {"E":"--", "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":Ltot, "M":M_cal, "n+1":"-"}
        derivadas.loc[len(derivadas)] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

        # Calculamos la siguiente capa
        i += 1


######################################################################
############################ Fase A.1.3. #############################
######################################################################

fase = "A.1.3."         # Fase actual ("A.1.3.")

loop1 = True

while loop1:

    # Aplicamos el paso 2
    P_est, T_est = paso2(modelo, derivadas, h, i)

    loop2 = True

    while loop2:

        loop3 = True

        while loop3:
            
            # Aplicamos el paso 3
            M_cal, fM = paso3(modelo, derivadas, P_est, T_est, modelo["M"][i], h, i)

            # Aplicamos el paso 4
            P_cal, fP = paso4(modelo, derivadas, P_est, T_est, M_cal, h, i)

            # Aplicamos el paso 5
            if pasoX(P_cal, P_est):
                loop3 = False
            else:
                P_est = P_cal

        # Aplicamos el paso 6
        L_cal, fL, ciclo = paso6(modelo, derivadas, P_cal, T_est, modelo["L"][i], h, i)

        # Aplicamos el paso 7
        T_cal, fT = paso7(modelo, derivadas, P_cal, T_est, L_cal, h, i)

        # Aplicamos el paso 8
        if pasoX(T_cal, T_est):
            loop2 = False
        else:
            T_est = T_cal

    # Aplicamos el paso 9
    n1 = paso9(P_cal, T_cal, fP, fT)

    # Comparamos n+1 con 2.5
    if paso10(n1):
        loop1 = False
    else:
        # Añadimos las variables calculadas al modelo
        r = modelo["r"][i] + h
        modelo.loc[len(modelo)] = {"E":ciclo, "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":L_cal, "M":M_cal, "n+1":n1}
        derivadas.loc[len(derivadas)] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

        # Calculamos la siguiente capa
        i += 1

print(modelo)