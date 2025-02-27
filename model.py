import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
import numpy as np

class Model:
    """
    Stellar-interior numerical model. Object that once initializated estimates the
    numerical value of radius, pressure, temperature, mass, luminosity and density
    throughout a star with given initial parameters.
    
    ## Parameters
        * Mtot (float, default = 5.0) : total mass of the star.
        * Rtot (float, default = 11.5) : total radius of the star.
        * Ltot (float, default = 70.0) : total luminosity of the star.
        * Tc (float, default = 2.0) : central temperature of the star.

    ## Methods
        * error() : function that returns the total relative error of the numerical 
        calculation.
        * visualize(x_axis="r", which=["P", "T", "L", "M", "rho"], merge=False) : 
        function that enables the graphical representation of the calculated variables.

    ## Units
    In order to properly estimate the variables of the star, the unit system adopted 
    for parameter input and results interpretation varies with respect to CGS.

    * radius (r)                         ->   1e10 cm
    * pressure (P)                       ->   1e15 din cm^-2
    * temperature (T)                    ->   1e7 K
    * mass (M)                           ->   1e33 g
    * luminosity (L)                     ->   1e33 erg s^-1
    * density (rho)                      ->   1 g cm^-3
    * energy generation rate (epsilon)   ->   1 erg g^-1 s^-1
    * opacity (kappa)                    ->   1 cm^2 g^-1
    """

    # Model initializer
    def __init__(self, Mtot=5.0, Rtot=11.5, Ltot=70.0, Tc=2.0):

        # Variable parameters
        self.Mtot = Mtot
        self.Rtot = Rtot
        self.Ltot = Ltot
        self.Tc = Tc

        # Constant parameters
        self.X = 0.75
        self.Y = 0.22
        self.Z = 1 - self.X - self.Y
        self.mu = 1/(2*self.X + 3*self.Y/4 + self.Z/2)
        self.R = 8.31447 * 10**7

        # Defining the energy generation rate table using modified units (10e7 K)
        self.epsilon_df = pd.DataFrame(data=[("PP", 0.40, 0.60, -6.84, 6.0),
                                ("PP", 0.60, 0.95, -6.04, 5.0), 
                                ("PP", 0.95, 1.20, -5.56, 4.5), 
                                ("PP", 1.20, 1.65, -5.02, 4.0), 
                                ("PP", 1.65, 2.40, -4.40, 3.5), 
                                ("CN", 1.20, 1.60, -22.2, 20.0),
                                ("CN", 1.60, 2.25, -19.8, 18.0), 
                                ("CN", 2.25, 2.75, -17.1, 16.0), 
                                ("CN", 2.75, 3.60, -15.6, 15.0), 
                                ("CN", 3.60, 5.00, -12.5, 13.0) ], 
                            columns=["Cycle", "Tmin", "Tmax", "log10(epsilon1)", "nu"])

        # Calculating the model variables
        self._calculate()

    # Defining the print instruction
    def __repr__(self):
        return self.model.__repr__()
    def __str__(self):
        return self.model.to_string()
    
    # Defining the get method
    def get(self, variable="all"):
        """
        Function to access variables.
        
        ## Parameters
            * variable (string, default = 'all'): 

            If default ('all'), a Data Frame object is returned containing the calculated 
            values of the variables. For queries on specific variables you must enter one
            of the following strings: 'r', 'P', 'T', 'L', 'M' and 'rho'
        """
        if variable == "all":
            return self.model
        else:
            return self.model[variable]
            
    
    # Defining the error method
    def error(self):
        """
        Returns total relative error of the numerical calculation of the star-interior 
        model.
        """
        return self.totalRelativeError

    # Defining the plot method
    def visualize(self, x_axis="r", which=["P", "T", "L", "M", "rho"], merge=False):
        """
        Function to graph the calculated variables.

        ## Parameters
            * x_axis (string, default = 'r'): 

            String to select the independent variable of the plot. It can only be one of
            the variables calculated with the model: 'r', 'P', 'T', 'L', 'M' and 'rho'.

            * which (array-like, default = ['P', 'T', 'L', 'M', 'rho']): 

            Array-like containing the dependent variables desirable to plot in string format.
            Supports the same values as x_axis: 'r', 'P', 'T', 'L', 'M' and 'rho'.

            * merge (bool, default = False): 
            If True, it plots all variables specified in 'which' normalized in the same figure.
            If False, it plots all variables specified in 'which' without normalizing in
            different figures.
        """

        # Changing font, colors and figure size
        plt.rcParams["font.family"] = "serif"
        # plt.style.use("tableau-colorblind10")
        # 'classic' 'grayscale' 'seaborn-v0_8-colorblind' 'seaborn-v0_8-dark-palette' 'seaborn-v0_8-muted'
        figSize = (10, 7)
        
        # Defining titles and labels for each variable 
        plots = pd.DataFrame(data=[("Radius", "r / $10^{10}$ cm"),
                                ("Pressure", "P / $10^{15}$ din cm$^{-2}$"), 
                                ("Temperature", "T / $10^{7}$ K"),
                                ("Luminosity", "L / $10^{33}$ erg s$^{-1}$"),
                                ("Mass", "M / $10^{33}$ g"),
                                ("Density","$\\rho$ / g cm$^{-3}$")],
                            columns=["title", "label"],
                            index=["r", "P", "T", "L", "M", "rho"])

        # Calculating the x value for which the transition to the convective zone occurs
        transition = self.model[self.model["fase"] == "CONVEC"].iloc[0][x_axis] 

        # All curves in the same figure (normalized plots)
        if merge:
            plt.figure(figsize=figSize)
            for variable in which:
                plt.plot(self.model[x_axis], self.model[variable]/self.model[variable].max(), label=plots.loc[variable]["title"], linewidth=1.5)

            # Customizing the plot
            plt.title("Stellar-interior model", fontsize=20)          # Title
            plt.xlabel(plots.loc[x_axis]["label"], fontsize=16)       # x axis label
            plt.ylabel("Normalized magnitude", fontsize=16)           # y axis label
            plt.tick_params(axis="both", labelsize=14, length=10)     # Numbering size
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))   # Minor ticks (x axis)
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))   # Minor ticks (y axis)
            
            # Marking the convective zone of the star
            plt.axvspan(-1, transition, color="gray", alpha=0.2, label="Convective zone")
            plt.legend(fontsize=12)                                                         # Leyend
            plt.grid(which="major", linestyle="-", linewidth=1, visible=True)               # Major grid
            plt.grid(which="minor", linestyle=":", linewidth=0.5, visible=True, alpha=0.5)  # Minor grid

        # Curves in different figures (without normalization)
        else:
            for variable in which:
                plt.figure(figsize=figSize)
                plt.plot(self.model[x_axis], self.model[variable], color="k", linewidth=1.5)

                # Customizing the plots for each figure
                plt.title(plots.loc[variable]["title"], fontsize=20)      # Title
                plt.xlabel(plots.loc[x_axis]["label"], fontsize=16)       # x axis label
                plt.ylabel(plots.loc[variable]["label"], fontsize=16)     # y axis label
                plt.tick_params(axis="both", labelsize=14, length=10)   # Numbering size
                plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))   # Minor ticks (x axis)
                plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))   # Minor ticks (y axis)

                # Marking the convective zone of the star for each figure
                plt.axvspan(-1, transition, color="gray", alpha=0.2, label="Convective zone")
                plt.legend(fontsize=16)                                                         # Leyend
                plt.grid(which="major", linestyle="-", linewidth=1, visible=True)               # Major grid
                plt.grid(which="minor", linestyle=":", linewidth=0.5, visible=True, alpha=0.5)  # Minor grid

        plt.show()


    def _calculate(self):

        ################################################################################
        ################################################################################
        ######################## ----- Defining equations ----- ########################
        ################################################################################
        ################################################################################

        # Density
        def rho(self, P, T):
            """
            Equation (6)
            Both pressure and temperature are expressed in the modified unit system.
            """
            return (self.mu/self.R)*((P*1e15)/(T*1e7))
        
        # Polytrope constant
        def K_pol(self, P, T):
            return P/(T**2.5)
        
        def polytrope(self, K, T):
            return K*(T**2.5)

        # Energy generation rate
        def epsilon(self, P, T, row):
            """
            Equation (9) using DataFrames.
            Given a row from the energy generation rate table, a value for pressure (P) and a value for
            temperature (T) it returns the energy generation rate (epsilon) and the values for X1 and X2.
            Temperature is expressed in the modified unit system.
            """

            # Selecting the appropriate row parameters of the energy generation rate table
            cycle, _, _, log10epsilon1, nu = row.to_list()

            # Depending on whether the cycle is 'PP' or 'CN', X1 and X2 change
            if cycle == "PP":
                    X1, X2 = self.X, self.X
            elif cycle == "CN":
                    X1, X2 = self.X, self.Z/3

            return (10**log10epsilon1)*X1*X2*rho(self, P, T)*(T*10)**nu, X1, X2

        # Calculating the most efficient energy generation cycle for a given pressure and temperature
        def calculate_optimal_cycle(self, P, T):
            """
            For a given pressure and temperature it retruns the values of epsilon1, X1, X2 and nu
            using the most efficient energy generation cycle.
            """

            # If fusion temperature has not been reached yet
            if T < self.epsilon_df.iloc[0,1]:
                return (0.0, 0.0, 0.0, 0.0, "--")

            # Otherwise, only the cycle that generates more energy should be considered
            else:
                # Applying a temperature filter to the energy generation rate table
                filter = (self.epsilon_df["Tmin"] <= T) & (T < self.epsilon_df["Tmax"])
                epsilon_values = self.epsilon_df[filter].reset_index(drop=True)
                # Calcultaing epsilon, X1 and X2 for each row
                epsilon_values[["epsilon", "X1", "X2"]] = epsilon_values.apply(lambda row: epsilon(self, P, T, row), axis=1, result_type="expand")
                # Selecting the row for which the energy generation rate is greater
                index = epsilon_values["epsilon"].idxmax()
                cycle, _, _, log10epsilon1, nu, _, X1, X2 = epsilon_values.loc[index].to_list()
                # Returning epsilon1, X1, X2, nu and the cycle
                return 10**log10epsilon1, X1, X2, nu, cycle


        ################################################################################
        # ---------------------- Fundamental equations (Table 3) --------------------- #
        ################################################################################
        
        # ------------------------------ Radiative case ------------------------------ #

        def dMdr_rad(self, r, P, T):
            """
            Equation (18)
            """
            Cm = 0.01523*self.mu
            return Cm * (P*r**2) / T

        def dPdr_rad(self, r, P, T, M):
            """
            Equation (19)
            """
            Cp = 8.084*self.mu
            return - Cp * P * M / (T * (r**2))

        def dLdr_rad(self, r, P, T):
            """
            Equation (20)
            """
            epsilon1, X1, X2, nu, cycle = calculate_optimal_cycle(self, P, T)
            Cl = 0.01845*epsilon1*X1*X2*(10**nu)*self.mu**2
            return Cl * ((P*r)**2) * (T**(nu-2)), cycle

        def dTdr_rad(self, r, P, T, L):
            """
            Equation (21)
            """
            Ct = 0.01679*self.Z*(1+self.X)*self.mu*self.mu
            return -Ct * (L*(P**2)) / ((T**8.5)*(r**2))


        # ------------------------------ Convective case ----------------------------- #

        def dMdr_conv(self, r, T, K):
            """
            Equation (22)
            """
            Cm = 0.01523*self.mu
            return Cm * K*(T**1.5)*r**2

        def dPdr_conv(self, r, T, M, K):
            """
            Equation (23)
            """
            if r == 0.0:
                return 0.0
            else:
                Cp = 8.084*self.mu
                return -Cp * K*(T**1.5)*M / (r**2)

        def dLdr_conv(self, r, P, T, K):
            """
            Equation (24)
            """
            epsilon1, X1, X2, nu, cycle = calculate_optimal_cycle(self, P, T)
            Cl = 0.01845*epsilon1*X1*X2*(10**nu)*self.mu**2
            return Cl * (K**2)*(T**(3+nu))*r**2, cycle

        def dTdr_conv(self, r, M):
            """
            Equation (25)
            """
            if r == 0.0:
                return 0.0
            else:
                Ct = 3.234*self.mu
                return -Ct * M / (r**2)


        ################################################################################
        # --------------------- Equations to find initial values --------------------- #
        ################################################################################

        # ----------------------- Initial values ​​on the surface ---------------------- #

        # Radiative envelope
        def initial_surface_T(self, r):
            """
            Equation (35)
            """
            A1 = 1.9022*self.mu*self.Mtot
            return A1*(1/r - 1/self.Rtot)

        # Convective envelope
        def initial_surface_P(self, r, T):
            """
            Equation (36)
            """
            A2 = 10.645*((self.Mtot/self.Ltot)/(self.mu*self.Z*(1+self.X)))**0.5
            return A2*T**4.25


        # ----------------------- Initial values ​​at the center ----------------------- #

        def initial_center_M(self, r, T, K):
            """
            Equation (43)
            """
            return 0.005077*self.mu*K*(T**1.5)*r**3

        def initial_center_L(self, r, P, T, K):
            """
            Equation (44)
            """
            epsilon1, X1, X2, nu, cycle = calculate_optimal_cycle(self, P, T)
            return 0.006150*epsilon1*X1*X2*(10**nu)*(self.mu**2)*(K**2)*(self.Tc**(3+nu))*(r**3), cycle

        def initial_center_T(self, r, K):
            """
            Equation (45)
            """
            return self.Tc - 0.008207*(self.mu**2)*K*(self.Tc**1.5)*(r**2)

        def initial_center_P(self, r, T, K):
            """
            Equation (46)
            """
            return K*T**2.5


        ################################################################################
        ################################################################################
        ########################## ----- Defining steps ----- ##########################
        ################################################################################
        ################################################################################

        # ---------------------------------- Step 1 ---------------------------------- #

        def step1(self, shell, h):
            """
            Given the variables of the star for a shell it returns the value of the 
            radius for the next shell.
            """

            return shell["r"] + h


        # ---------------------------------- Step 2 ---------------------------------- #

        def step2(self, shell, derivatives, h, i):
            """
            Given the variables of the star for a shell and the derivatives of the current 
            and two previous shells it estimates the value of the pressure and temperature
            for the next shell.
            """

            P = shell["P"]                                                      # P for the i shell
            T = shell["T"]                                                      # T for the i shell
            fP = derivatives["fP"][i]                                           # fP for the i shell
            fT = derivatives["fT"][i]                                           # fT for the i shell
            AP1 = h * (fP -  derivatives["fP"][i-1])                            # AP1 for the i shell
            AP2 = h * (fP - 2*derivatives["fP"][i-1] + derivatives["fP"][i-2])  # AP2 for the i shell
            AT1 = h * (fT -  derivatives["fT"][i-1])                            # AT1 for the i shell

            P_est = P + h*fP + AP1/2 + 5*AP2/12    # Estimated pressure for the i+1 shell
            T_est = T + h*fT + AT1/2               # Estimated temperature for the i+1 shell

            return P_est, T_est


        # ---------------------------------- Step 3 ---------------------------------- #

        def step3(self, shell, derivatives, r, P, T, M, h, i):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            and temperature for the next shell and the derivatives of the current and two 
            previous shells it calculates the value of the mass and its derivative for the 
            next shell.
            """

            fM = dMdr_rad(self, r, P, T)             # fM for the i+1 shell
            AM1 = h * (fM - derivatives["fM"][i])    # AM1 for the i+1 shell

            M_cal = M + h*fM - AM1/2                 # Calculated mass for the i+1 shell

            return M_cal, fM


        # ---------------------------------- Step 4 ---------------------------------- #

        def step4(self, shell, derivatives, r, P, T, M, h, i):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            temperature and mass for the next shell and the derivatives of the current and 
            two previous shells it calculates the value of the pressure and its derivative 
            for the next shell.
            """

            fP = dPdr_rad(self, r, P, T, M)          # fP for the i+1 shell
            AP1 = h * (fP - derivatives["fP"][i])    # AP1 for the i+1 shell
            P = shell["P"]                           # P for the i shell

            P_cal = P + h*fP - AP1/2                 # Calculated pressure for the i+1 shell

            return P_cal, fP


        # ---------------------------------- Step 5 ---------------------------------- #

        def step5(self, P_cal, P_est):
            """
            Given two values for the pressure, it tests if its relative difference is
            smaller than the maximum relative error.
            """

            maximumRelativeError = 0.0001
            
            return abs(P_cal-P_est)/abs(P_cal) < maximumRelativeError


        # ---------------------------------- Step 6 ---------------------------------- #

        def step6(self, shell, derivatives, r, P, T, L, h, i):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure 
            and temperature for the next shell and the derivatives of the current and two 
            previous shells it calculates the value of the luminosity and its derivative 
            for the next shell. It also returns which energy generation cycle produces 
            the energy.
            """

            fL, cycle = dLdr_rad(self, r, P, T)                                 # fL for the i+1 shell
            AL1 = h * (fL - derivatives["fL"][i])                               # AL1 for the i+1 shell
            AL2 = h * (fL - 2*derivatives["fL"][i] + derivatives["fL"][i-1])    # AL2 for the i+1 shell
            
            L_cal = L + h*fL - AL1/2  - AL2/12        # Calculated luminosity for the i+1 shell

            return L_cal, fL, cycle


        # ---------------------------------- Step 7 ---------------------------------- #

        def step7(self, shell, derivatives, r, P, T, L, h, i):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            temperature and luminosity for the next shell and the derivatives of the current 
            and two previous shells it calculates the value of the temperature and its 
            derivative for the next shell. 
            """

            fT = dTdr_rad(self, r, P, T, L)          # fT for the i+1 shell
            AT1 = h * (fT - derivatives["fT"][i])    # AT1 for the i+1 shell
            T = shell["T"]                           # T for the i shell

            T_cal = T + h*fT - AT1/2                 # Calculated temperature for the i+1 shell

            return T_cal, fT


        # -------------------------------- Step 7 bis -------------------------------- #

        def step7bis(self, shell, derivatives, r, M, h, i):
            """
            Given the variables of the star for a shell, an estimation of radius and mass 
            for the next shell and the derivatives of the current and two previous shells it 
            calculates the value of the temperature and its derivative for the next shell. 
            """

            fT = dTdr_conv(self, r, M)               # fT for the i+1 shell
            AT1 = h * (fT - derivatives["fT"][i])    # AT1 for the i+1 shell
            T = shell["T"]                           # T for the i shell

            T_cal = T + h*fT - AT1/2                 # Calculated temperature for the i+1 shell

            return T_cal, fT


        # ---------------------------------- Step 8 ---------------------------------- #

        def step8(self, T_cal, T_est):
            """
            Given two values for the temperature, it tests if its relative difference 
            is smaller than the maximum relative error.
            """

            maximumRelativeError = 0.0001
            
            return abs(T_cal-T_est)/abs(T_cal) < maximumRelativeError


        # ---------------------------------- Step 9 ---------------------------------- #

        def step9(self, P, T, fP, fT):
            """
            Given calculated values for temperature and pressure and its derivatives 
            for the next shell it calculates the n+1 coefficient.
            """
            return T*fP/(P*fT)


        # ---------------------------------- Step 10 --------------------------------- #

        def step10(self, n1):
            """
            Given a value for the n+1 coeficient it tests if it is smaller than 2.5.
            """
            return n1 <= 2.5


        # ---------------------------------- Step X ---------------------------------- #

        def stepX(self, x_cal, x_est):
            """
            Given two values for a variable, it tests if its relative difference is 
            smaller than the maximum relative error.
            """

            maximumRelativeError = 0.0001   
            
            return abs(x_cal-x_est)/abs(x_cal) < maximumRelativeError


        ################################################################################
        # ---------------------- Computing total relative error ---------------------- #
        ################################################################################

        def calculate_total_relative_error(self, X_rad, X_conv):
            """
            Given two array-like objects computes the total relative error
            """
            return (sum(((X_rad - X_conv)/X_rad)**2))**0.5


        ################################################################################
        ################################################################################
        #################### ----- Performing the calculation ----- ####################
        ################################################################################
        ################################################################################

        # Customizing the DataFrame options
        pd.set_option("colheader_justify", "center", "display.max_rows", 1000)
        pd.options.display.float_format = '{:.7f}'.format

        # Defining the DataFrames that will store the data
        outer = pd.DataFrame(data = [], columns=["E", "fase", "r", "P", "T", "L", "M", "rho", "n+1"])
        outer_derivatives = pd.DataFrame(data = [], columns=["fP", "fT", "fL", "fM"])
        inside = pd.DataFrame(data = [], columns=["E", "fase", "r", "P", "T", "L", "M", "rho", "n+1"])
        inside_derivatives = pd.DataFrame(data = [], columns=["fP", "fT", "fL", "fM"])


        ################################################################################
        # ----------------------- First three shells (surface) ----------------------- #
        ################################################################################

        Rini = 0.9*self.Rtot            # Defining the initial radius
        shellBreaks = 100               # Setting the shell partition of Rini
        innerShells = shellBreaks + 1   # Setting the number of shells from Rini
        h = - Rini/shellBreaks          # Defining the radius step between shells (negative)
        r = Rini                        # Initializing the radius variable

        for i in range(3):

            # Calculating and storing the values for the first three shells
            T = initial_surface_T(self, r)
            P = initial_surface_P(self, r, T)
            M = self.Mtot
            L = self.Ltot
            outer.loc[i] = {"E":"--", "fase":"INICIO", "r":r, "P":P, "T":T, "L":L, "M":M, "rho":rho(self, P, T), "n+1":"-"}

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_rad(self, r, P, T, L)
            fP = dPdr_rad(self, r, P, T, M)
            fM = 0.0
            fL = 0.0
            outer_derivatives.loc[i] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Moving forward one step
            r += h


        ################################################################################
        # --------------------------------- Fase A.1. -------------------------------- #
        ################################################################################

        fase = "RADIAT"         # Current fase

        while True:             # First loop

            # Values for the current shell and derivatives of the current and two previous shells
            shell = outer.loc[i]
            derivatives = outer_derivatives.loc[i-2:i]

            # Performing the first step
            r = step1(self, shell, h)
            # Performing the second step
            P_est, T_est = step2(self, shell, derivatives, h, i)

            while True:         # Second loop

                while True:     # Third loop
                    
                    # Performing the third step
                    M_cal, fM = step3(self, shell, derivatives, r, P_est, T_est, shell["M"], h, i)
                    # Performing the fourth step
                    P_cal, fP = step4(self, shell, derivatives, r, P_est, T_est, M_cal, h, i)
                    # Performing the fifth step
                    if stepX(self, P_cal, P_est):
                        break   # Third loop break
                    else:
                        P_est = P_cal

                # Performing the sixth step
                L_cal, fL, cycle = step6(self, shell, derivatives, r, P_cal, T_est, shell["L"], h, i)
                # Performing the seventh step
                T_cal, fT = step7(self, shell, derivatives, r, P_cal, T_est, L_cal, h, i)
                # Performing the eighth step
                if stepX(self, T_cal, T_est):
                    break       # Second loop break
                else:
                    T_est = T_cal
                    P_est = P_cal

            # Performing the ninth step
            n1 = step9(self, P_cal, T_cal, fP, fT)
            # Performing the tenth step
            if step10(self, n1):
                break           # First loop break
            else:
                # Storing values and derivatives
                outer.loc[i+1] = {"E":cycle, "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":L_cal, "M":M_cal, "rho":rho(self, P_cal, T_cal), "n+1":n1}
                outer_derivatives.loc[i+1] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

                # Moving forward one step
                i += 1

        # Storing how many shells are left to compute
        convec_index = innerShells - len(outer)


        ################################################################################
        # --------------------------------- Fase A.2. -------------------------------- #
        ################################################################################

        # Unnecesary calculations

        ################################################################################
        # ------------------------- First three shells (center) ---------------------- #
        ################################################################################
        
        h = Rini/shellBreaks            # Defining the radius step between shells (positive)
        r = 0.0                         # Initializing the radius variable
        K = K_pol(self, P_cal, T_cal)   # Calculating the polytrope coefficient        
        
        for i in range(3):

            # Calculating and storing the values for the first three shells
            T = initial_center_T(self, r, K)
            P = initial_center_P(self, r, T, K)
            M = initial_center_M(self, r, self.Tc, K)
            L, _ = initial_center_L(self, r, P, self.Tc, K)
            inside.loc[i] = {"E":"--", "fase":"CENTRO", "r":r, "P":P, "T":T, "L":L, "M":M, "rho":rho(self, P ,T), "n+1":"-"}

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_conv(self, r, M)
            fP = dPdr_conv(self, r, T, M, K)
            fM = dMdr_conv(self, r, T, K)
            fL, _ = dLdr_conv(self, r, P, T, K)
            inside_derivatives.loc[i] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Moving forward one step
            r += h


        ################################################################################
        # -------------------------------- Next shells ------------------------------- #
        ################################################################################

        fase = "CONVEC"         # Current fase

        # Iterating over the final layers
        for i in range(2, convec_index):

            # Values for the current shell and derivatives of the current and two previous shells
            shell = inside.loc[i]
            derivatives = inside_derivatives.loc[i-2:i]

            # Performing the first step
            r = step1(self, shell, h)
            # Performing the second step bis
            _, T_est = step2(self, shell, derivatives, h, i)

            while True:         # Loop 
                
                # Estimating pressure using the polytrope constant
                P_est = polytrope(self, K, T_est)
                # Performing the third step
                M_cal, fM = step3(self, shell, derivatives, r, P_est, T_est, shell["M"], h, i)
                # Performing the seventh step bis
                if shell["r"] + h < 1e-6:           # Very small values are considered zero
                    T_cal = T_est
                elif shell["r"] + h > 0.0:
                    T_cal, fT = step7bis(self, shell, derivatives, r, M_cal, h, i)
                # Performing the eighth step
                if stepX(self, T_cal, T_est):
                    break       # Loop break
                else:
                    T_est = T_cal

            # Calculating pressure using the polytrope constant
            P_cal = polytrope(self, K, T_cal)
            # Performing the sixth step
            L_cal, fL, cycle = step6(self, shell, derivatives, r, P_cal, T_cal, shell["L"], h, i)
            # Calculating the derivative of the pressure for the i shell
            fP = dPdr_conv(self, r, T_cal, M_cal, K)
            # Storing values and derivatives
            inside.loc[i+1] = {"E":cycle, "fase":fase, "r":r, "P":P_cal, "T":T_cal, "L":L_cal, "M":M_cal, "rho":rho(self, P_cal, T_cal), "n+1":"-"}
            inside_derivatives.loc[i+1] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Moving forward one step
            i += 1

        # Storing data from both surface integration and center integration
        rad = np.array(outer.loc[len(outer)-1, "P":"M"].to_list(), dtype=float)
        conv = np.array(inside.loc[len(inside)-1, "P":"M"].to_list(), dtype=float)

        # Redefining inside dataframe index so r=0.0 corresponds to i=shellBreaks
        inside.index = range(shellBreaks, shellBreaks - convec_index - 1, -1)


        ################################################################################
        # ------------------------- More superficial shells -------------------------- #
        ################################################################################

        r = Rini + h            # Initializing the radius variable
        i = -1                      

        # Iterating until reaching the total radius
        while self.Rtot - r > 0.0:

            # Calculating and storing the values for the first three shells
            T = initial_surface_T(self, r)
            P = np.real(initial_surface_P(self, r, T))    # Since r ~ Rtot, the program crashes if we do not take the real part
            M = self.Mtot
            L = self.Ltot
            outer.loc[i] = {"E":"--", "fase":"^^^^^^", "r":r, "P":P, "T":T, "L":L, "M":M, "rho":rho(self, P, T), "n+1":"-"}

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_rad(self, r, P, T, L)
            fP = dPdr_rad(self, r, P, T, M)
            fM = 0.0
            fL = 0.0
            outer_derivatives.loc[i] = {"fP":fP, "fT":fT, "fL":fL, "fM":fM}

            # Moving backward one step
            r += h
            i -= 1


        ################################################################################
        # ---------------------- Model and total relative error ---------------------- #
        ################################################################################

        # Calculating total relative error using surface integration data and center integration data
        self.totalRelativeError = calculate_total_relative_error(self, rad, conv)*100 

        # Joning outer and inside DataFrames using outer middle point
        self.model = pd.concat([outer, inside.iloc[:-1]]).sort_index()
        # Joning outer and inside DataFrames using inside middle point
        # self.model = pd.concat([outer.sort_index().iloc[:-1], inside]).sort_index()