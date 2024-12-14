import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class Modelo:
    
    # Definición del inicializador del modelo
    def __init__(self, Mtot=5.0, Rtot=11.5, Ltot=70.0, Tc=2.0):

        # Definición de parámetros variables
        self.Mtot = Mtot
        self.Rtot = Rtot
        self.Ltot = Ltot
        self.Tc = Tc

        # Definición de parámetros constantes
        self.X = 0.75
        self.Y = 0.22
        self.Z = 1 - self.X - self.Y
        self.mu = 1/(2*self.X + 3*self.Y/4 + self.Z/2)
        self.R = 8.31447 * 10**7

        # Definición de la Tabla 1 (ritmo de generación de energía)
        self.epsilon_df = pd.DataFrame(data=[("pp", 0.40, 0.60, -6.84, 6.0),
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

        # Ejecución del cálculo
        self._calcular()

    # Definición de la instrucción print sobre el modelo
    def __repr__(self):
        return self.modelo.__repr__()
    def __str__(self):
        return self.modelo.to_string()
    
    # Definición de la función para obtener el error relativo total
    def error(self):
        return self.error_relativo_total
    
    # Definición de la función para la representación gráfica
    def grafica(self, x_axis="r", which=["P", "T", "L", "M", "rho"], merge=False):
        """
        Función para graficar las variables calculadas por el modelo.

        Parámetros:
            - x_axis (string, default = 'r'): 
                String que indica qué variable se utiliza como variable independiente. Admite los valores 
                siguientes: 'r', 'P', 'T', 'L', 'M' y 'rho'.

            - which (iterable, default = ["P", "T", "L", "M", "rho"]): 
                Lista (o iterable) que contiene las variables dependientes que se quieren graficar en formato 
                string. Admite los mismos valores que x_axis.

            - merge (bool, default = False):
                Si es True, representa todas las variables especificadas en which en una misma figura y
                normalizadas.
                Si es False, representa todas las variables especificadas en which en diferentes figuras
                sin normalizar.
        """
        
        # Definimos los títulos y etiquetas para cada variable
        plots = pd.DataFrame(data=[("Radio", "r / $10^{10}$ cm"),
                                ("Presión", "P / $10^{15}$ din cm$^{-2}$"), 
                                ("Temperatura", "T / $10^{7}$ K"),
                                ("Luminosidad", "L / $10^{33}$ erg s$^{-1}$"),
                                ("Masa", "M / $10^{33}$ g"),
                                ("Densiadad","$\\rho$ / g$cm^{-3}$")],
                            columns=["title", "label"],
                            index=["r", "P", "T", "L", "M", "rho"])

        # Si queremos las gráficas en la misma figura (curvas normalizadas)
        if merge:
            plt.figure(figsize=(8, 6))
            for variable in which:
                plt.plot(self.modelo[x_axis], self.modelo[variable]/self.modelo[variable].max(), label=plots.loc[variable]["title"])
                plt.xlabel(plots.loc[x_axis]["label"])
                plt.legend()
                plt.grid(visible=True)

        # Si queremos utilizar diferentes figuras (curvas sin normalizar)
        else:
            for variable in which:
                plt.figure(figsize=(8, 6))
                plt.plot(self.modelo[x_axis], self.modelo[variable])
                plt.title(plots.loc[variable]["title"])
                plt.xlabel(plots.loc[x_axis]["label"])
                plt.ylabel(plots.loc[variable]["label"])
                plt.grid(visible=True)

        plt.show()


    def _calcular(self):

        ################################################################################
        ################################################################################
        ##################### ----- Definición de ecuaciones ----- #####################
        ################################################################################
        ################################################################################

        # Densidad
        def rho(self, P, T):
            """
            Ecuación (6)
            Convertimos la presión y la temperatura a las unidades estándar
            """
            return (self.mu/self.R)*((P*1e15)/(T*1e7))
        
        # Constante del polítropo
        def K_pol(self, P_dado, T_dado):
            return P_dado/(T_dado**2.5)
        
        def politropo(self, K, T_dado):
            return K*(T_dado**2.5)

        # Generación de energía
        def epsilon(self, P, T, row):
            """
            Ecuación (9) para DataFrames
            """
            ciclo, _, _, log10epsilon1, nu = row.to_list()

            if ciclo == "pp":
                    X1, X2 = self.X, self.X
            elif ciclo == "CN":
                    X1, X2 = self.X, self.Z/3

            return (10**log10epsilon1)*X1*X2*rho(self, P, T)*(T*10)**nu, X1, X2

        def calculo_tabla1(self, P, T):
            """
            Dada una presión P y temperatura T  devuelve los valores de 
            epsilon1, X1, X2 y nu junto con el ciclo que genera la energía.
            """
            # Aplicamos un filtro al DataFrame de la Tabla 1
            filtro = (self.epsilon_df["Tmin"] <= T) & (T < self.epsilon_df["Tmax"])
            epsilon_values = self.epsilon_df[filtro].reset_index(drop=True)

            # Si el DataFrame está vacío, devolvemos 0.0
            if epsilon_values.empty:
                return (0.0, 0.0, 0.0, 0.0, "--")
            # En caso contrario, tomamos el ciclo que genere mayor energía
            else:
                # Calculamos el valor de epsilon, X1 y X2 para cada fila
                epsilon_values[["epsilon", "X1", "X2"]] = epsilon_values.apply(lambda row: epsilon(self, P, T, row), axis=1, result_type='expand')
                # Tomamos el índice de la fila para el que epsilon es mayor y extraemos sus elementos
                index = epsilon_values["epsilon"].idxmax()
                ciclo, _, _, log10epsilon1, nu, _, X1, X2 = epsilon_values.loc[index].to_list()
                # Devolvemos los valores de epsilon1, X1, X2, nu y el ciclo que lo genera
                return 10**log10epsilon1, X1, X2, nu, ciclo


        ################################################################################
        # ------ Ecuaciones fundamentales en las unidades del modelo (Tabla 3) ------- #
        ################################################################################
        
        # ------------------------------ Caso radiativo ------------------------------ #

        def dMdr_rad(self, r, P, T):
            """
            Ecuación (18)
            """
            Cm = 0.01523*self.mu
            return Cm * (P*r**2) / T

        def dPdr_rad(self, r, P, T, M):
            """
            Ecuación (19)
            """
            Cp = 8.084*self.mu
            return - Cp * P*M / (T*r**2)

        def dLdr_rad(self, r, P, T):
            """
            Ecuación (20)
            """
            epsilon1, X1, X2, nu, ciclo = calculo_tabla1(self, P, T)
            Cl = 0.01845*epsilon1*X1*X2*(10**nu)*self.mu**2
            return Cl * ((P*r)**2) * (T**(nu-2)), ciclo

        def dTdr_rad(self, r, P, T, L):
            """
            Ecuación (21)
            """
            Ct = 0.01679*self.Z*(1+self.X)*self.mu*self.mu
            return -Ct * (L*(P**2)) / ((T**8.5)*(r**2))


        # ------------------------------ Caso convectivo ----------------------------- #

        def dMdr_conv(self, r, T, K):
            """
            Ecuación (22)
            """
            Cm = 0.01523*self.mu
            return Cm * K*(T**1.5)*r**2

        def dPdr_conv(self, r, T, M, K):
            """
            Ecuación (23)
            """
            if r == 0.0:
                return 0.0
            else:
                Cp = 8.084*self.mu
                return -Cp * K*(T**1.5)*M / (r**2)

        def dLdr_conv(self, r, P, T, K):
            """
            Ecuación (24)
            """
            epsilon1, X1, X2, nu, ciclo = calculo_tabla1(self, P, T)
            Cl = 0.01845*epsilon1*X1*X2*(10**nu)*self.mu**2
            return Cl * (K**2)*(T**(3+nu))*r**2, ciclo

        def dTdr_conv(self, r, M):
            """
            Ecuación (25)
            """
            if r == 0.0:
                return 0.0
            else:
                Ct = 3.234*self.mu
                return -Ct * M / (r**2)


        ################################################################################
        # --------------- Ecuaciones para encontrar valores iniciales ---------------- #
        ################################################################################

        # -------------------- Valores iniciales en la superficie -------------------- #

        # Envoltura radiativa
        def T_inicial_superficie(self, r):
            """
            Ecuación (35)
            """
            A1 = 1.9022*self.mu*self.Mtot
            return A1*(1/r - 1/self.Rtot)

        # Envoltura convectiva
        def P_inicial_superficie(self, r, T):
            """
            Ecuación (36)
            """
            A2 = 10.645*((self.Mtot/self.Ltot)/(self.mu*self.Z*(1+self.X)))**0.5
            return A2*T**4.25


        # ---------------------- Valores iniciales en el centro ---------------------- #

        def M_inicial_centro(self, r, T, K):
            """
            Ecuación (43)
            """
            return 0.005077*self.mu*K*(T**1.5)*r**3

        def L_inicial_centro(self, r, P, T, K):
            """
            Ecuación (44)
            """
            epsilon1, X1, X2, nu, ciclo = calculo_tabla1(self, P, T)
            return 0.006150*epsilon1*X1*X2*(10**nu)*(self.mu**2)*(K**2)*(self.Tc**(3+nu))*(r**3), ciclo

        def T_inicial_centro(self, r, K):
            """
            Ecuación (45)
            """
            return self.Tc - 0.008207*(self.mu**2)*K*(self.Tc**1.5)*(r**2)

        def P_inicial_centro(self, r, T, K):
            """
            Ecuación (46)
            """
            return K*T**2.5


        ################################################################################
        ################################################################################
        ####################### ----- Definición de pasos ----- ########################
        ################################################################################
        ################################################################################

        # ---------------------------------- Paso 2 ---------------------------------- #

        def paso2(self, capa, derivadas, h, i):

            # Calculamos la estimación de P y T en la capa i+1
            P = capa["P"]          # P en la capa i
            T = capa["T"]          # T en la capa i
            fP = derivadas["fP"][i]     # fP en la capa i
            fT = derivadas["fT"][i]     # fT en la capa i
            AP1 = h * (derivadas["fP"][i] -  derivadas["fP"][i])                            # AP1 en la capa i
            AP2 = h * (derivadas["fP"][i] - 2*derivadas["fP"][i-1] + derivadas["fP"][i-2])  # AP2 en la capa i
            AT1 = h * (derivadas["fT"][i] -  derivadas["fT"][i])                            # AT1 en la capa i

            P_est = P + h*fP + AP1/2 + 5*AP2/12     # P estimado en la capa i+1
            T_est = T + h*fT + AT1/2                # T estimado en la capa i+1

            return P_est, T_est


        # ---------------------------------- Paso 3 ---------------------------------- #

        def paso3(self, capa, derivadas, P_dado, T_dado, M_dado, h, i):

            # Calculamos fM en la capa i+1 (P, T y M se dan en la capa i+1)

            r = capa["r"] + h     # r en la capa i+1
            fM = dMdr_rad(self, r, P_dado, T_dado)     # fM en la capa i+1


            # Calculamos M en la capa i+1 teniendo

            AM1 = h * (fM - derivadas["fM"][i])       # AM1 en la capa i+1

            M_cal = M_dado + h*fM - AM1/2             # M calculado en la capa i+1

            return M_cal, fM


        # ---------------------------------- Paso 4 ---------------------------------- #

        def paso4(self, capa, derivadas, P_dado, T_dado, M_dado, h, i):

            # Calculamos fP en la capa i+1 (P, T y M se dan en la capa i+1)

            r = capa["r"] + h                         # r en la capa i+1
            fP = dPdr_rad(self, r, P_dado, T_dado, M_dado)  # fP en la capa i+1


            # Calculamos P en la capa i+1

            P = capa["P"]                        # P en la capa i
            AP1 = h * (fP - derivadas["fP"][i])  # AP1 en la capa i+1

            P_cal = P + h*fP - AP1/2             # P calculado en la capa i+1

            return P_cal, fP


        # ---------------------------------- Paso 5 ---------------------------------- #

        def paso5(self, P_cal, P_est):

            # Comprobamos si el error relativo es menor al error relativo máximo

            Error_relativo_maximo = 0.0001
            
            return abs(P_cal-P_est)/abs(P_cal) < Error_relativo_maximo


        # ---------------------------------- Paso 6 ---------------------------------- #

        def paso6(self, capa, derivadas, P_dado, T_dado, L_dado, h, i):

            # Calculamos fL en la capa i+1 (P, T y L se dan en la capa i+1)

            r = capa["r"] + h     # r en la capa i+1
            fL, ciclo = dLdr_rad(self, r,P_dado,T_dado)       # fL en la capa i+1

            # Calculamos L en la capa i+1
            AL1 = h * (fL - derivadas["fL"][i])                            # AL1 en la capa i+1
            AL2 = h * (fL - 2*derivadas["fL"][i] + derivadas["fL"][i-1])        # AL2 en la capa i+1
            L_cal = L_dado + h*fL - AL1/2  - AL2/12                   # L calculado en la capa i+1

            return L_cal, fL, ciclo


        # ---------------------------------- Paso 7 ---------------------------------- #

        def paso7(self, capa, derivadas, P_dado, T_dado, L_dado, h, i):

            # Calculamos fT en la capa i+1 (P, T y L se dan en la capa i+1)

            r = capa["r"] + h                         # r en la capa i+1
            fT = dTdr_rad(self, r, P_dado, T_dado, L_dado)  # fT en la capa i+1


            # Calculamos T en la capa i+1

            T = capa["T"]                   # T en la capa i
            AT1 = h * (fT - derivadas["fT"][i])  # AT1 en la capa i+1

            T_cal = T + h*fT - AT1/2            # T calculado en la capa i+1

            return T_cal, fT


        # -------------------------------- Paso 7 bis -------------------------------- #

        def paso7bis(self, capa, derivadas, M_dado, h, i):
            
            # Calculamos fT en la capa i+1 (M se da en la capa i+1)

            r = capa["r"] + h          # r en la capa i+1
            fT = dTdr_conv(self, r, M_dado)  # fT en la capa i+1


            # Calculamos T en la capa i+1

            T = capa["T"]                   # T en la capa i
            AT1 = h * (fT - derivadas["fT"][i])     # AT1 en la capa i+1

            T_cal = T + h*fT - AT1/2            # T calculado en la capa i+1

            return T_cal, fT


        # ---------------------------------- Paso 8 ---------------------------------- #

        def paso8(self, T_cal, T_est):

            # Comprobamos si el error relativo es menor al error relativo máximo

            Error_relativo_maximo = 0.0001
            
            return abs(T_cal-T_est)/abs(T_cal) < Error_relativo_maximo


        # ---------------------------------- Paso 9 ---------------------------------- #

        def paso9(self, P_dado, T_dado, fP_dado, fT_dado):
            return T_dado*fP_dado/(P_dado*fT_dado)


        # ---------------------------------- Paso 10 --------------------------------- #

        def paso10(self, n1):
            return n1 <= 2.5


        # ---------------------------------- Paso X ---------------------------------- #

        def pasoX(self, x_cal, x_est):

            # Comprobamos si el error relativo es menor al error relativo máximo

            Error_relativo_maximo = 0.0001
            
            return abs(x_cal-x_est)/abs(x_cal) < Error_relativo_maximo


        ################################################################################
        # --------------------------- Error relativo total --------------------------- #
        ################################################################################

        def err_rel_total(self, X_rad, X_conv):
            """
            Devuelve el error relativo total para dos arrays (X_rad, X_conv)
            Cada array contiene el valor de P, T, L, M
            """
            return (sum(((X_rad - X_conv)/X_rad)**2))**0.5


        ################################################################################
        ################################################################################
        ###################### ----- Ejecución del cálculo ----- #######################
        ################################################################################
        ################################################################################

        pd.set_option('colheader_justify', 'center', "display.precision", 9, 'display.max_rows', 1000)
        
        exterior = pd.DataFrame(data = [], columns=["E", "fase", "r", "P", "T", "L", "M", "rho", "n+1"])
        derivadas_exterior = pd.DataFrame(data = [], columns=["fP", "fT", "fL", "fM"])
        interior = pd.DataFrame(data = [], columns=["E", "fase", "r", "P", "T", "L", "M", "rho", "n+1"])
        derivadas_interior = pd.DataFrame(data = [], columns=["fP", "fT", "fL", "fM"])


        ################################################################################
        # ---------------------- Primeras tres capas (superficie) -------------------- #
        ################################################################################

        Rini = 0.9*self.Rtot
        paso = 100
        h = - Rini/paso
        r = Rini

        for i in range(3):

            # Calculamos y almacenamos los valores del modelo
            T = T_inicial_superficie(self, r)
            P = P_inicial_superficie(self, r, T)
            M = self.Mtot
            L = self.Ltot
            exterior.loc[i] = {"E":"--", "fase":"INICIO", "r":r, "P":P, "T":T, "L":L, "M":M, "rho":rho(self, P, T), "n+1":"-"}

            # Calculamos y almacenamos los valores de las f_i (derivadas)
            fT = dTdr_rad(self, r, P, T, L)
            fP = dPdr_rad(self, r, P, T, M)
            fM = 0.0
            fL = 0.0
            derivadas_exterior.loc[i] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Aumentamos el valor del radio
            r += h


        ################################################################################
        # -------------------------------- Fase A.1.1. ------------------------------- #
        ################################################################################

        fase = "A.1.1."         # Fase actual (A.1.1.)

        loop1 = True

        while loop1:

            # Capa actual y derivadas actuales y una y dos capas atrás
            capa = exterior.loc[i]
            derivadas = derivadas_exterior.loc[i-2:i]

            # Aplicamos el paso 2
            P_est, T_est = paso2(self, capa, derivadas, h, i)

            loop2 = True

            while loop2:

                loop3 = True

                while loop3:
                    
                    # Aplicamos el paso 4
                    P_cal, fP = paso4(self, capa, derivadas, P_est, T_est, self.Mtot, h, i)

                    # Aplicamos el paso 5
                    if pasoX(self, P_cal, P_est):
                        loop3 = False
                    else:
                        P_est = P_cal

                # Aplicamos el paso 7
                T_cal, fT = paso7(self, capa, derivadas, P_cal, T_est, self.Ltot, h, i)

                # Aplicamos el paso 8
                if pasoX(self, T_cal, T_est):
                    loop2 = False
                else:
                    T_est = T_cal

            # Aplicamos el paso 3
            M_cal, fM = paso3(self, capa, derivadas, P_cal, T_cal, self.Mtot, h, i)

            # Comparamos la masa calculada con la masa total
            if not pasoX(self, M_cal, self.Mtot):
                loop1 = False
            else:
                # Añadimos las variables calculadas al modelo (asumimos que M y L permanecen constantes)
                r = capa["r"] + h
                exterior.loc[i+1] = {"E":"--", "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":self.Ltot, "M":self.Mtot, "rho":rho(self, P_cal, T_cal), "n+1":"-"}
                derivadas_exterior.loc[i+1] = {"fP":fP, "fT":fT, "fL":0.0, "fM":fM}

                # Calculamos la siguiente capa
                i += 1


        ################################################################################
        # -------------------------------- Fase A.1.2. ------------------------------- #
        ################################################################################

        fase = "A.1.2."         # Fase actual ("A.1.2.")

        loop1 = True

        while loop1:

            # Capa actual y derivadas actuales y una y dos capas atrás
            capa = exterior.loc[i]
            derivadas = derivadas_exterior.loc[i-2:i]

            # Aplicamos el paso 2
            P_est, T_est = paso2(self, capa, derivadas, h, i)

            loop2 = True

            while loop2:

                loop3 = True

                while loop3:
                    
                    # Aplicamos el paso 3
                    M_cal, fM = paso3(self, capa, derivadas, P_est, T_est, capa["M"], h, i)

                    # Aplicamos el paso 4
                    P_cal, fP = paso4(self, capa, derivadas, P_est, T_est, M_cal, h, i)

                    # Aplicamos el paso 5
                    if pasoX(self, P_cal, P_est):
                        loop3 = False
                    else:
                        P_est = P_cal

                # Aplicamos el paso 7
                T_cal, fT = paso7(self, capa, derivadas, P_cal, T_est, self.Ltot, h, i)

                # Aplicamos el paso 8
                if pasoX(self, T_cal, T_est):
                    loop2 = False
                else:
                    T_est = T_cal

            # Aplicamos el paso 6
            L_cal, fL, ciclo = paso6(self, capa, derivadas, P_cal, T_cal, self.Ltot, h, i)

            # Comparamos la luminosidad calculada con la luminosidad total
            if not pasoX(self, L_cal, self.Ltot):
                loop1 = False
            else:
                # Añadimos las variables calculadas al modelo (asumimos que L permanece constante)
                r = capa["r"] + h
                exterior.loc[i+1] = {"E":"--", "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":self.Ltot, "M":M_cal, "rho":rho(self, P_cal,T_cal), "n+1":"-"}
                derivadas_exterior.loc[i+1] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

                # Calculamos la siguiente capa
                i += 1


        ################################################################################
        # -------------------------------- Fase A.1.3. ------------------------------- #
        ################################################################################

        fase = "A.1.3."         # Fase actual ("A.1.3.")

        loop1 = True

        while loop1:

            # Capa actual y derivadas actuales y una y dos capas atrás
            capa = exterior.loc[i]
            derivadas = derivadas_exterior.loc[i-2:i]

            # Aplicamos el paso 2
            P_est, T_est = paso2(self, capa, derivadas, h, i)

            loop2 = True

            while loop2:

                loop3 = True

                while loop3:
                    
                    # Aplicamos el paso 3
                    M_cal, fM = paso3(self, capa, derivadas, P_est, T_est, capa["M"], h, i)

                    # Aplicamos el paso 4
                    P_cal, fP = paso4(self, capa, derivadas, P_est, T_est, M_cal, h, i)

                    # Aplicamos el paso 5
                    if pasoX(self, P_cal, P_est):
                        loop3 = False
                    else:
                        P_est = P_cal

                # Aplicamos el paso 6
                L_cal, fL, ciclo = paso6(self, capa, derivadas, P_cal, T_est, capa["L"], h, i)

                # Aplicamos el paso 7
                T_cal, fT = paso7(self, capa, derivadas, P_cal, T_est, L_cal, h, i)

                # Aplicamos el paso 8
                if pasoX(self, T_cal, T_est):
                    loop2 = False
                else:
                    T_est = T_cal

            # Aplicamos el paso 9
            n1 = paso9(self, P_cal, T_cal, fP, fT)

            # Comparamos n+1 con 2.5
            if paso10(self, n1):
                loop1 = False
                K = K_pol(self, P_cal, T_cal)          # Cálculo de K en la capa i+1 (Fase A.2.)
                
            else:
                # Añadimos las variables calculadas al modelo
                r = capa["r"] + h
                exterior.loc[i+1] = {"E":ciclo, "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":L_cal, "M":M_cal, "rho":rho(self, P_cal, T_cal), "n+1":n1}
                derivadas_exterior.loc[i+1] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

                # Calculamos la siguiente capa
                i += 1

        # Almacenamos el número de capas que quedan por calcular

        convec_index = 101 - len(exterior)


        ################################################################################
        # --------------------------------- Fase A.2. -------------------------------- #
        ################################################################################

        # Cálculos innecesarios


        ################################################################################
        # ------------------------ Primeras tres capas (centro) ---------------------- #
        ################################################################################

        h = Rini/paso   # Redefinimos el paso
        r = 0.0         # Integramos desde r = 0.0

        # Iteramos las primeras tres capas
        for i in range(3):

            # Calculamos y almacenamos los valores del modelo
            T = T_inicial_centro(self, r, K)
            P = P_inicial_centro(self, r, T, K)
            M = M_inicial_centro(self, r, T, K)
            L, _ = L_inicial_centro(self, r, P, T, K)
            interior.loc[i] = {"E":"--", "fase":"CENTRO", "r":r, "P":P, "T":T, "L":L, "M":M, "rho":rho(self, P ,T), "n+1":"-"}

            # Calculamos y almacenamos los valores de las f_i (derivadas)
            fT = dTdr_conv(self, r, M)
            fP = dPdr_conv(self, r, T, M, K)
            fM = dMdr_conv(self, r, T, K)
            fL, _ = dLdr_conv(self, r, P, T, K)
            derivadas_interior.loc[i] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Aumentamos el valor del radio
            r += h


        ################################################################################
        # ----------------------------- Capas posteriores ---------------------------- #
        ################################################################################

        fase = "CONVEC"     # Fase actual ("CONVEC")

        # Iteramos sobre las capas que faltan
        for i in range(2, convec_index):

            # Capa actual y derivadas actuales y una y dos capas atrás
            capa = interior.loc[i]
            derivadas = derivadas_interior.loc[i-2:i]

            # Aplicamos el paso 2 bis
            _, T_est = paso2(self, capa, derivadas, h, i)

            loop2 = True

            while loop2:

                # Aplicamos el cálculo de P con la constante del politropo K
                P_est = politropo(self, K, T_est)

                # Aplicamos el paso 3
                M_cal, fM = paso3(self, capa, derivadas, P_est, T_est, capa["M"], h, i)

                # Aplicamos el paso 7 bis
                if capa["r"] + h < 1e-6:  # Tomamos como 0 valores muy pequeños
                    T_cal = T_est
                elif capa["r"] + h > 0.0:
                    T_cal, fT = paso7bis(self, capa, derivadas, M_cal, h, i)

                # Aplicamos el paso 8
                if pasoX(self, T_cal, T_est):
                    loop2 = False
                else:
                    T_est = T_cal

            # Aplicamos el cálculo de P con la constante del politropo K
            P_cal = politropo(self, K, T_cal)

            # Aplicamos el paso 6 
            L_cal, fL, ciclo = paso6(self, capa, derivadas, P_cal, T_cal, capa["L"], h, i)

            # Calculamos el radio y fP
            r = capa["r"] + h
            fP = dPdr_conv(self, r, T_cal, M_cal, K)

            # Añadimos las variables calculadas al modelo
            interior.loc[i+1] = {"E":ciclo, "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":L_cal, "M":M_cal, "rho":rho(self, P_cal, T_cal), "n+1":"-"}
            derivadas_interior.loc[i+1] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Calculamos la siguiente capa
            i += 1


        # Almacenamos los datos que da la integración desde el centro y desde la superficie
        rad = np.array(exterior.loc[len(exterior)-1, "P":"M"].to_list(), dtype=float)
        conv = np.array(interior.loc[len(interior)-1, "P":"M"].to_list(), dtype=float)

        # Redefinimos los índices del interior para que el centro sea i = paso
        interior.index = range(paso, paso - convec_index - 1, -1)


        ################################################################################
        # ------------------------- Capas más superficiales -------------------------- #
        ################################################################################

        # Reutilizamos el código del primer apartado

        r = Rini + h
        i = -1

        # Calculamos hasta llegar a Rtot
        while self.Rtot - r > 0.0:

            # Calculamos y almacenamos los valores del modelo
            T = T_inicial_superficie(self, r)
            P = np.real(P_inicial_superficie(self, r, T))    # Al trabajar con r casi Rtot, tomamos la parte real
            M = self.Mtot
            L = self.Ltot
            exterior.loc[i] = {"E":"--", "fase":"^^^^^^", "r":r, "P":P, "T":T, "L":L, "M":M, "rho":rho(self, P, T), "n+1":"-"}

            # Calculamos y almacenamos los valores de las f_i (derivadas)
            fT = dTdr_rad(self, r, P, T, L)
            fP = dPdr_rad(self, r, P, T, M)
            fM = 0.0
            fL = 0.0
            derivadas_exterior.loc[i] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Aumentamos el valor del radio y subimos una capa
            r += h
            i -= 1


        ################################################################################
        # ---------------------- Modelo y error relativo total ----------------------- #
        ################################################################################

        # A partir de los datos que da la integración desde el centro y desde la superficie 
        # calculamos el error relativo total
        self.error_relativo_total = err_rel_total(self, rad, conv)*100 

        # Juntamos ambas partes tomando como punto intermedio el calculado desde el exterior
        self.modelo = pd.concat([exterior, interior.iloc[:-1]]).sort_index()
        # Juntamos ambas partes tomando como punto intermedio el calculado desde el interior
        # self.modelo = pd.concat([exterior.sort_index().iloc[:-1], interior]).sort_index()

Modelo(Rtot=10.93, Ltot=73.57, Tc=1.95).grafica(x_axis="M")

# print(Modelo(Rtot=10.93, Ltot=73.57, Tc=1.95).error())