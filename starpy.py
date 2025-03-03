import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import pandas as pd
import numpy as np

class Star:
    """
    Stellar-interior numerical model. Object that once initializated estimates the
    numerical value of radius, pressure, temperature, mass, luminosity and density
    throughout a star with given initial parameters.
    
    ## Parameters
        * Mtot (float, default = 5.0) : total mass of the star.
        * Rtot (float, default = 11.5) : total radius of the star.
        * Ltot (float, default = 70.0) : total luminosity of the star.
        * Tc (float, default = 2.0) : central temperature of the star.
        * X (float, default = 0.75) :  fraction of star mass in H.
        * Y (float, default = 0.22 : fraction of mass in He.

    ## Methods
        * error() : function that returns the total relative error of the numerical 
        calculation.
        * visualize(x_axis="r", which=["P", "T", "L", "M", "rho"], merge=False, solar_units=True, figsize=(10, 7)) : 
        function that enables the graphical representation of the calculated variables.

    ## Units
    In order to properly estimate the variables of the star, the unit system adopted 
    for parameter input and results interpretation varies with respect to CGS.

    * radius (r)                         ->   1e10 cm
    * pressure (P)                       ->   1e15 dyn cm^-2
    * temperature (T)                    ->   1e7 K
    * mass (M)                           ->   1e33 g
    * luminosity (L)                     ->   1e33 erg s^-1
    * density (rho)                      ->   1 g cm^-3
    * energy generation rate (epsilon)   ->   1 erg g^-1 s^-1
    * opacity (kappa)                    ->   1 cm^2 g^-1
    """

    # Model initializer
    def __init__(self, Mtot=5.0, Rtot=11.5, Ltot=70.0, Tc=2.0, X=0.75, Y=0.22):

        # Variable parameters
        self.Mtot = Mtot                                # Total mass of the star
        self.Rtot = Rtot                                # Total radius of the star
        self.Ltot = Ltot                                # Total luminosity of the star
        self.Tc = Tc                                    # Central temperature of the star

        # Constant parameters
        self.X = X                                      # Fraction of star mass in H
        self.Y = Y                                      # Fraction of star mass in He
        self.Z = 1 - self.X - self.Y                    # Metalicity
        self.mu = 1/(2*self.X + 3*self.Y/4 + self.Z/2)  # Mean molecular weight

        # Calculating the model variables
        self.model, self.totalRelativeError = self._calculate()

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

    # Defining the plot method for star variables
    def visualize(self, x_axis="r", which=["P", "T", "L", "M", "rho"], merge=False, solar_units=True, figsize=(10, 7)):
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


            * solar_units (bool, default = True)

            If True, all plots will be graphed using solar units.
            If False, all plots will be graphed using the model units.


            * figsize (two-dimensional array-like, default = (10, 7))

            Two-dimensional array-like for a better customization on the figures size.
        """

        # Changing font
        plt.rcParams["font.family"] = "serif"

        # Solar units are used
        if solar_units:
            # Defining sun parameters
            Msun = 1.9884    # 1e33 g
            Rsun = 6.957     # 1e10 cm
            Lsun = 3.842     # 1e33 erg s^-1
            # Redefining model units
            modelData = self.model
            modelData["r"] = modelData["r"] / Rsun  # cm
            modelData["P"] = modelData["P"] * 1e15  # dyn cm^-2
            modelData["T"] = modelData["T"] * 1e7   # K
            modelData["M"] = modelData["M"] / Msun  # g
            modelData["L"] = modelData["L"] / Lsun  # erg s^-1

            # Defining titles and labels for each variable 
            plots = pd.DataFrame(data=[("Radius", "r / R$_{\\odot}$"),
                                    ("Pressure", "P / dyn cm$^{-2}$"), 
                                    ("Temperature", "T / K"),
                                    ("Luminosity", "L / L$_{\\odot}$"),
                                    ("Mass", "M / M$_{\\odot}$"),
                                    ("Density", "$\\rho$ / g cm$^{-3}$")],
                                columns=["title", "label"],
                                index=["r", "P", "T", "L", "M", "rho"])
        
        # Model units are used
        else:
            # Using model units
            modelData = self.model
            # Defining titles and labels for each variable 
            plots = pd.DataFrame(data=[("Radius", "r / $10^{10}$ cm"),
                                    ("Pressure", "P / $10^{15}$ dyn cm$^{-2}$"), 
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
            plt.figure(figsize=figsize)
            for variable in which:
                plt.plot(self.model[x_axis], self.model[variable]/self.model[variable].max(), label=plots.loc[variable]["title"], linewidth=1.5)

            # Customizing the plot
            plt.title("Stellar-interior model", fontsize=20)                    # Title
            plt.xlabel(plots.loc[x_axis]["label"], fontsize=16)                 # x axis label
            plt.ylabel("Normalized magnitude", fontsize=16)                     # y axis label
            plt.tick_params(axis="both", labelsize=14)                          # Numbering size
            plt.gca().tick_params(direction="in", which="major", length=8)      # Major ticks size and orientation
            plt.gca().tick_params(direction="in", which="minor", length=3)      # Minor ticks size and orientation
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (x axis)
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (y axis)
            
            # Marking the convective zone of the star
            plt.axvspan(-1, transition, color="gray", alpha=0.2, label="Convective zone")
            plt.legend(fontsize=12)                                                         # Leyend
            plt.grid(which="major", linestyle="-", linewidth=1, visible=True)               # Major grid
            plt.grid(which="minor", linestyle=":", linewidth=0.5, visible=True, alpha=0.5)  # Minor grid

        # Curves in different figures (without normalization)
        else:
            for variable in which:
                plt.figure(figsize=figsize)
                plt.plot(self.model[x_axis], self.model[variable], color="k", linewidth=1.5)

                # Customizing the plots for each figure
                plt.title(plots.loc[variable]["title"], fontsize=20)                # Title
                plt.xlabel(plots.loc[x_axis]["label"], fontsize=16)                 # x axis label
                plt.ylabel(plots.loc[variable]["label"], fontsize=16)               # y axis label
                plt.tick_params(axis="both", labelsize=14)                          # Numbering size
                plt.gca().tick_params(direction="in", which="major", length=8)      # Major ticks size and orientation
                plt.gca().tick_params(direction="in", which="minor", length=3)      # Minor ticks size and orientation
                plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (x axis)
                plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (y axis)

                # Marking the convective zone of the star for each figure
                plt.axvspan(-1, transition, color="gray", alpha=0.2, label="Convective zone")
                plt.legend(fontsize=16)                                                         # Leyend
                plt.grid(which="major", linestyle="-", linewidth=1, visible=True)               # Major grid
                plt.grid(which="minor", linestyle=":", linewidth=0.5, visible=True, alpha=0.5)  # Minor grid

        plt.show()

    # Defining the Temperature-Density Diagram method
    def TDD(self, figsize=(8,7)):
        """
        Function to plot the Temperature-Density Diagram, i.e. the values throughout the star for
        temperature and density. Several regions are distinguished depending on the dominant pressure.

        I: ideal gas. II: degeneracy. III: relativistic degeneracy. IV: radiation pressure.

        ## Parameters
            * figsize (two-dimensional array-like, default = (10, 7))

            Two-dimensional array-like for a better customization on the figures size.
        """

        # Defining fundamental constants and star constants
        R = 8.31447 * 10**7             # Universal gas constant
        a = 7.56578e-15                 # Radiation density constant
        mu_e = 2/(1+self.X)             # Mean molecular electron weight

        # Selecting T and rho
        T = np.log10(self.get("T")*1e7)       # Star log10 temperature in K
        rho = np.log10(self.get("rho"))       # Star log10 density in g/cm^3

        # Changing font and figure size. Setting the plot limits
        plt.rcParams["font.family"] = "serif"
        plt.figure(figsize=figsize)
        xlims = [min(T) / 1.1, 10]
        ylims = [min(rho) * 1.1, 10.5]

        # Defining polytropic constants
        K0 = R/self.mu                      # Ideal gas polytrope constant
        K1 = (1.0036e13)/(mu_e**(5/3))      # Degeneracy polytrope constant
        K2 = (1.2435e15)/(mu_e**(4/3))      # Relativistic degeneracy polytrope constant
        K3 = a/3                            # Radiation pressure polytrope constant

        # Defining transition zones
        logT_1_2 = np.array([xlims[0], 9], dtype=float)             # Ideal gas - Degeneracy (x axis)
        logT_1_3 = np.array([9, xlims[1]], dtype=float)             # Ideal gas - Relativistic degeneracy (x axis)
        logT_2_3 = np.array([xlims[0], xlims[1]], dtype=float)      # Degeneracy - Relativistic degeneracy (x axis)

        logrho_I_II = 1.5*np.log10(K0/K1) + 1.5*logT_1_2            # Ideal gas - Degeneracy (y axis)
        logrho_I_III = 3*np.log10(K0/K2) + 3*logT_1_3               # Ideal gas - Relativistic degeneracy (y axis)
        logrho_II_III = 3*np.log10(K2/K1)*np.ones((2,))             # Degeneracy - Relativistic degeneracy (y axis)
        logrho_I_IV = np.log10(K3/(10*K0)) + 3*logT_2_3             # Ideal gas - Radiation pressure (y axis)

        # Graphing the transition zones
        plt.plot(logT_1_2, logrho_I_II, color = "k", linestyle="solid")             # Ideal gas - Degeneracy
        plt.plot(logT_1_3, logrho_I_III, color = "k", linestyle="solid")            # Ideal gas - Relativistic degeneracy
        plt.plot(logT_1_2, logrho_II_III, color = "k", linestyle="solid")           # Ideal gas - Radiation pressure
        plt.plot(logT_2_3, logrho_I_IV, color = "k", linestyle="solid")             # Degeneracy - Relativistic degeneracy
        # Fillin between the transition zones
        plt.fill_between(logT_2_3, logrho_I_IV, ylims[1], color="#d9d9d9")          # Ideal gas
        plt.fill_between(logT_2_3, logrho_I_IV, ylims[0], color="#f2f2f2")          # Radiation pressure
        plt.fill_between(logT_1_2, logrho_I_II, logrho_II_III, color="#a6a6a6")     # Degeneracy
        plt.fill_between(logT_1_2, logrho_II_III, ylims[1], color="#5c5c5c")        # Relativistic degeneracy
        plt.fill_between(logT_1_3, logrho_I_III, ylims[1], color="#5c5c5c")         # Relativistic degeneracy

        # Adding text to the plot
        plt.text(0.10, 0.30, "Ideal gas", transform=plt.gca().transAxes, fontsize=12, color="brown")
        plt.text(0.15, 0.65, "Degeneracy", transform=plt.gca().transAxes, fontsize=12, color="brown")
        plt.text(0.25, 0.90, "Relativistic degeneracy", transform=plt.gca().transAxes, fontsize=12, color="brown")
        plt.text(0.65, 0.2, "Radiation pressure", transform=plt.gca().transAxes, fontsize=12, color="brown")

        # Graphing the star variables in the diagram
        plt.plot(T, rho, color="#006769")
        plt.scatter(T.iloc[len(T)-1], rho.iloc[len(rho)-1], s=60, color="#006769", marker="o", label="Star center")
        plt.scatter(T.iloc[0], rho.iloc[0], s=60, color="#006769", marker="d", label="Star surface")

        # Adding title and labels. Setting the ticks parameters
        plt.title("Temperature-Density Diagram", fontsize=20)               # Title
        plt.xlabel("Log [T(K)]", fontsize=16)                               # x axis label
        plt.ylabel("Log [$\\rho$(g cm$^{-3})$]", fontsize=16)               # y axis label
        plt.tick_params(axis="both", labelsize=14)                          # Numbering size
        plt.gca().tick_params(direction="in", which="major", length=8)      # Major ticks size and orientation
        plt.gca().tick_params(direction="in", which="minor", length=3)      # Minor ticks size and orientation
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))             # Setting minor ticks (x axis)
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))             # Setting minor ticks (y axis)
        plt.xlim(xlims)                                                     # x limits
        plt.ylim(ylims)                                                     # y limits
        plt.legend(loc="upper right", fontsize=14)                          # Legend

        plt.show()

    # Defining the Hertzsprung–Russell Diagram method
    def HRD(self):
        ...

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
            R = 8.31447 * 10**7
            return (self.mu/R)*((P*1e15)/(T*1e7))
        
        # Polytrope constant
        def K_pol(self, P, T):
            return P/(T**2.5)
        
        def polytrope(self, K, T):
            return K*(T**2.5)

        # Calculating the most efficient energy generation cycle for a given pressure and temperature
        def calculate_optimal_cycle(self, P, T):
            """
            For a given pressure and temperature it retruns the values of epsilon1, X1, X2 and nu
            using the most efficient energy generation cycle.
            """

            # Calculating PP cycle parameters
            if T < 0.40:
                epsilon1_PP, nu_PP = 0.0, 0.0
            elif T < 0.60:
                epsilon1_PP, nu_PP = 10**-6.84, 6.0
            elif T < 0.95:
                epsilon1_PP, nu_PP = 10**-6.04, 5.0
            elif T < 1.20:
                epsilon1_PP, nu_PP = 10**-5.56, 4.5
            elif T < 1.65:
                epsilon1_PP, nu_PP = 10**-5.02, 4.0
            elif T < 2.40:
                epsilon1_PP, nu_PP = 10**-4.40, 3.5
            else:
                epsilon1_PP, nu_PP = 0.0, 0.0

            # Calculating CN cycle parameters
            if T < 1.20:
                epsilon1_CN, nu_CN = 0.0, 0.0
            elif T < 1.60:
                epsilon1_CN, nu_CN = 10**-22.2, 20.0
            elif T < 2.25:
                epsilon1_CN, nu_CN = 10**-19.8, 18.0
            elif T < 2.75:
                epsilon1_CN, nu_CN = 10**-17.1, 16.0
            elif T < 3.60:
                epsilon1_CN, nu_CN = 10**-15.6, 15.0
            elif T < 5.00:
                epsilon1_CN, nu_CN = 10**-12.5, 13.0
            else:
                epsilon1_CN, nu_CN = 0.0, 0.0

            # Selecting the cycle for which the energy generation rate is greater
            epsilon_PP = epsilon1_PP * (self.X*self.X) * rho(self, P, T) * (T*10)**nu_PP
            epsilon_CN = epsilon1_CN * (self.X*self.Z/3) * rho(self, P, T) * (T*10)**nu_CN
            
            # If there is no energy generation
            if epsilon_PP == 0.0 and epsilon_CN == 0.0:
                return 0.0, 0.0, 0.0, 0.0, "--"
            # If the energy generation rate given by PP cycle is greater
            if epsilon_PP > epsilon_CN:
                return epsilon1_PP, self.X, self.X, nu_PP, "PP"
            # If the energy generation rate given by CN cycle is greater
            else:
                return epsilon1_CN, self.X, self.Z/3, nu_CN, "CN"


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

            return shell[2] + h


        # ---------------------------------- Step 2 ---------------------------------- #

        def step2(self, shell, derivatives, h):
            """
            Given the variables of the star for a shell and the derivatives of the current 
            and two previous shells it estimates the value of the pressure and temperature
            for the next shell.
            """

            P = shell[3]                                                    # P for the i shell
            T = shell[4]                                                    # T for the i shell
            fP = derivatives[2][0]                                          # fP for the i shell
            fT = derivatives[2][1]                                          # fT for the i shell
            AP1 = h * (fP -  derivatives[1][0])                             # AP1 for the i shell
            AP2 = h * (fP - 2*derivatives[1][0] + derivatives[0][0])        # AP2 for the i shell
            AT1 = h * (fT -  derivatives[1][1])                             # AT1 for the i shell

            P_est = P + h*fP + AP1/2 + 5*AP2/12    # Estimated pressure for the i+1 shell
            T_est = T + h*fT + AT1/2               # Estimated temperature for the i+1 shell

            return P_est, T_est


        # ---------------------------------- Step 3 ---------------------------------- #

        def step3(self, shell, derivatives, r, P, T, M, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            and temperature for the next shell and the derivatives of the current and two 
            previous shells it calculates the value of the mass and its derivative for the 
            next shell.
            """

            fM = dMdr_rad(self, r, P, T)             # fM for the i+1 shell
            AM1 = h * (fM - derivatives[2][3])       # AM1 for the i+1 shell

            M_cal = M + h*fM - AM1/2                 # Calculated mass for the i+1 shell

            return M_cal, fM


        # ---------------------------------- Step 4 ---------------------------------- #

        def step4(self, shell, derivatives, r, P, T, M, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            temperature and mass for the next shell and the derivatives of the current and 
            two previous shells it calculates the value of the pressure and its derivative 
            for the next shell.
            """

            fP = dPdr_rad(self, r, P, T, M)          # fP for the i+1 shell
            AP1 = h * (fP - derivatives[2][0])       # AP1 for the i+1 shell
            P = shell[3]                             # P for the i shell

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

        def step6(self, shell, derivatives, r, P, T, L, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure 
            and temperature for the next shell and the derivatives of the current and two 
            previous shells it calculates the value of the luminosity and its derivative 
            for the next shell. It also returns which energy generation cycle produces 
            the energy.
            """

            fL, cycle = dLdr_rad(self, r, P, T)                              # fL for the i+1 shell
            AL1 = h * (fL - derivatives[2][2])                               # AL1 for the i+1 shell
            AL2 = h * (fL - 2*derivatives[2][2] + derivatives[1][2])         # AL2 for the i+1 shell
            
            L_cal = L + h*fL - AL1/2  - AL2/12        # Calculated luminosity for the i+1 shell

            return L_cal, fL, cycle


        # ---------------------------------- Step 7 ---------------------------------- #

        def step7(self, shell, derivatives, r, P, T, L, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            temperature and luminosity for the next shell and the derivatives of the current 
            and two previous shells it calculates the value of the temperature and its 
            derivative for the next shell. 
            """

            fT = dTdr_rad(self, r, P, T, L)          # fT for the i+1 shell
            AT1 = h * (fT - derivatives[2][1])       # AT1 for the i+1 shell
            T = shell[4]                             # T for the i shell

            T_cal = T + h*fT - AT1/2                 # Calculated temperature for the i+1 shell

            return T_cal, fT


        # -------------------------------- Step 7 bis -------------------------------- #

        def step7bis(self, shell, derivatives, r, M, h):
            """
            Given the variables of the star for a shell, an estimation of radius and mass 
            for the next shell and the derivatives of the current and two previous shells it 
            calculates the value of the temperature and its derivative for the next shell. 
            """

            fT = dTdr_conv(self, r, M)               # fT for the i+1 shell
            AT1 = h * (fT - derivatives[2][1])       # AT1 for the i+1 shell
            T = shell[4]                             # T for the i shell

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

        # # Customizing the DataFrame options
        # pd.set_option("colheader_justify", "center", "display.max_rows", None)
        # pd.options.display.float_format = '{:.7f}'.format

        # # Defining the DataFrames that will store the data
        # outer = pd.DataFrame(data = [], columns=["E", "fase", "r", "P", "T", "L", "M", "rho", "n+1"])
        # outer_derivatives = pd.DataFrame(data = [], columns=["fP", "fT", "fL", "fM"])
        # inside = pd.DataFrame(data = [], columns=["E", "fase", "r", "P", "T", "L", "M", "rho", "n+1"])
        # inside_derivatives = pd.DataFrame(data = [], columns=["fP", "fT", "fL", "fM"])


        ################################################################################
        # ----------------------- First three shells (surface) ----------------------- #
        ################################################################################

        Rini = 0.9*self.Rtot            # Defining the initial radius
        shellBreaks = 100               # Setting the shell partition of Rini
        innerShells = shellBreaks + 1   # Setting the number of shells from Rini
        h = - Rini/shellBreaks          # Defining the radius step between shells (negative)
        r = Rini                        # Initializing the radius variable

        outer_data = []                      # Initializing the list that will store the data (surface integration)
        outer_derivatives = []          # Initializing the list that will store derivatives (surface integration)

        for i in range(3):

            # Calculating and storing the values for the first three shells
            T = initial_surface_T(self, r)
            P = initial_surface_P(self, r, T)
            M = self.Mtot
            L = self.Ltot
            outer_data.append(["--", "INICIO", r, P, T, L, M, rho(self, P, T), 0.0])

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_rad(self, r, P, T, L)
            fP = dPdr_rad(self, r, P, T, M)
            fM = 0.0
            fL = 0.0
            outer_derivatives.append([fP, fT, fL, fM])

            # Moving forward one step
            r += h


        ################################################################################
        # --------------------------------- Fase A.1. -------------------------------- #
        ################################################################################

        fase = "RADIAT"         # Current fase

        while True:             # First loop

            # Values for the current shell and derivatives of the current and two previous shells
            shell = outer_data[i]
            derivatives = outer_derivatives[i-2:i+1]

            # Performing the first step
            r = step1(self, shell, h)
            # Performing the second step
            P_est, T_est = step2(self, shell, derivatives, h)

            while True:         # Second loop

                while True:     # Third loop
                    
                    # Performing the third step
                    M_cal, fM = step3(self, shell, derivatives, r, P_est, T_est, shell[6], h)
                    # Performing the fourth step
                    P_cal, fP = step4(self, shell, derivatives, r, P_est, T_est, M_cal, h)
                    # Performing the fifth step
                    if stepX(self, P_cal, P_est):
                        break   # Third loop break
                    else:
                        P_est = P_cal

                # Performing the sixth step
                L_cal, fL, cycle = step6(self, shell, derivatives, r, P_cal, T_est, shell[5], h)
                # Performing the seventh step
                T_cal, fT = step7(self, shell, derivatives, r, P_cal, T_est, L_cal, h)
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
                outer_data.append([cycle, fase, r, P_cal, T_cal, L_cal, M_cal, rho(self, P_cal, T_cal), n1])
                outer_derivatives.append([fP, fT, fL, fM])

                # Moving forward one step
                i += 1

        # Storing how many shells are left to compute
        convec_index = innerShells - len(outer_data)


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
        
        inside_data = []                 # Initializing the list that will store the data (center integration)
        inside_derivatives = []          # Initializing the list that will store derivatives (center integration)

        for i in range(3):

            # Calculating and storing the values for the first three shells
            T = initial_center_T(self, r, K)
            P = initial_center_P(self, r, T, K)
            M = initial_center_M(self, r, self.Tc, K)
            L, _ = initial_center_L(self, r, P, self.Tc, K)
            inside_data.append(["--", "CENTER", r, P, T, L, M, rho(self, P, T), 0.0])

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_conv(self, r, M)
            fP = dPdr_conv(self, r, T, M, K)
            fM = dMdr_conv(self, r, T, K)
            fL, _ = dLdr_conv(self, r, P, T, K)
            inside_derivatives.append([fP, fT, fL, fM])

            # Moving forward one step
            r += h


        ################################################################################
        # -------------------------------- Next shells ------------------------------- #
        ################################################################################

        fase = "CONVEC"         # Current fase

        # Iterating over the final layers
        # for i in range(2, convec_index):
        for i in range(2, convec_index):

            # Values for the current shell and derivatives of the current and two previous shells
            shell = inside_data[i]
            derivatives = inside_derivatives[i-2:i+1]

            # Performing the first step
            r = step1(self, shell, h)
            # Performing the second step bis
            _, T_est = step2(self, shell, derivatives, h)

            while True:         # Loop 
                
                # Estimating pressure using the polytrope constant
                P_est = polytrope(self, K, T_est)
                # Performing the third step
                M_cal, fM = step3(self, shell, derivatives, r, P_est, T_est, shell[6], h)
                # Performing the seventh step bis
                if shell[2] + h < 1e-4:           # Very small values of "r" are considered zero
                    T_cal = T_est
                elif shell[2] + h > 0.0:
                    T_cal, fT = step7bis(self, shell, derivatives, r, M_cal, h)
                # Performing the eighth step
                if stepX(self, T_cal, T_est):
                    break       # Loop break
                else:
                    T_est = T_cal

            # Calculating pressure using the polytrope constant
            P_cal = polytrope(self, K, T_cal)
            # Performing the sixth step
            L_cal, fL, cycle = step6(self, shell, derivatives, r, P_cal, T_cal, shell[5], h)
            # Calculating the derivative of the pressure for the i shell
            fP = dPdr_conv(self, r, T_cal, M_cal, K)
            # Storing values and derivatives
            inside_data.append([cycle, fase, r, P_cal, T_cal, L_cal, M_cal, rho(self, P_cal, T_cal), 0.0])
            inside_derivatives.append([fP, fT, fL, fM])

            # Moving forward one step
            i += 1


        ################################################################################
        # ------------------------- More superficial shells -------------------------- #
        ################################################################################

        r = Rini + h            # Initializing the radius variable
        i = -1
        n = 0

        # Iterating until reaching the total radius
        while self.Rtot - r > 0.0:

            # Calculating and storing the values for the first three shells
            T = initial_surface_T(self, r)
            P = np.real(initial_surface_P(self, r, T))    # Since r ~ Rtot, the program crashes if we do not take the real part
            M = self.Mtot
            L = self.Ltot
            outer_data.insert(0, ["--", "^^^^^^", r, P, T, L, M, rho(self, P, T), 0.0])

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_rad(self, r, P, T, L)
            fP = dPdr_rad(self, r, P, T, M)
            fM = 0.0
            fL = 0.0
            outer_derivatives.insert(0, [fP, fT, fL, fM])

            # Moving backward one step
            r += h
            i -= 1
            n += 1


        ################################################################################
        # ---------------------- Model and total relative error ---------------------- #
        ################################################################################

        # Calculating total relative error using surface integration data and center integration data
        totalRelativeError = calculate_total_relative_error(self, np.array(outer_data[-1][2:7]), np.array(inside_data[-1][2:7]))*100 

        # Joning outer and inside data using outer middle point
        model = outer_data + list(reversed(inside_data[:-1]))
        # Joning outer and inside DataFrames using inside middle point
        # model = outer_data[:-1] + list(reversed(inside_data))

        # Returning a DataFrame object
        model = pd.DataFrame(model, columns=["E", "fase", "r", "P", "T", "L", "M", "rho", "n+1"])
        model.index -= n
        
        return model, totalRelativeError