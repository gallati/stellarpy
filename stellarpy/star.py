# Universidad Complutense de Madrid
# Facultad de Ciencias Físicas 
# Course 2024 / 2025
#
# Nuno Cerviño Luridiana
# ncervino@ucm.es 
#
# 
#
# Relation between variables and shell[i] index:
# i        :   0      1   2   3   4   5   6    7      8       9
# variable : cycle  fase  r   P   T   l   m   rho  epsilon  kappa

import matplotlib.pyplot as plt
from astropy import units as u
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patheffects as path_effects
import pandas as pd
import numpy as np
import os

class Star:
    """
    Stellar-interior numerical model. 
    
    Represents a star with a given mass, radius, luminosity, central temperature and chemical composition.

    ## Atributes:
        * Mtot (QuantityLike): Total mass of the star.
        * Rtot (QuantityLike): Total radius of the star.
        * Ltot (QuantityLike): Total luminosity of the star.
        * Tc (QuantityLike): Central temperature of the star.
        * X (QuantityLike): Fraction of star mass in H.
        * Y (QuantityLike): Fraction of mass in He.

    ## Methods:
        * get(): Returns the requested Star instance data.
        * parameters(): Returns Star instance atributes as a list following the order: [Mtot, Rtot, Ltot, Tc, X, Y]
        * redefine(): Redefines Star instance atributes.
        * error(): Returns the total relative error of the numerical calculation.
        * visualize(): Graphical representation of the calculated variables throughout the star.
        * TDD(): Graphical representation of the star in the Temperature-Density Diagram.
        * HR(): Graphical representation of the star in the Hertzsprung–Russell Diagram.

    ## Units:
        In order to properly estimate the variables of the star, the unit system adopted for internal calculations 
        of the model varies with respect to CGS. However, both input and output values of the model can be expressed
        in any unit system using Quantity objects from astropy.

        * radius (r)                         ->   1e10 cm
        * pressure (P)                       ->   1e15 dyn cm^-2
        * temperature (T)                    ->   1e7 K
        * mass (m)                           ->   1e33 g
        * luminosity (l)                     ->   1e33 erg s^-1
        * density (rho)                      ->   1 g cm^-3
        * energy generation rate (epsilon)   ->   1 erg g^-1 s^-1
        * opacity (kappa)                    ->   1 cm^2 g^-1
    """

    # Defining the model unit system
    UNIT_SYSTEM = {
        "radius": 1e10 * u.cm,
        "temperature": 1e7 * u.K,
        "mass": 1e33 * u.g,
        "luminosity": 1e33 * u.erg / u.s
    }

    # Model initializer
    def __init__(self, Mtot, Rtot, Ltot, Tc, X, Y):
        """
        Initializes a new instance of Star.

        ## Arguments:
            * Mtot (QuantityLike): Total mass of the star.
            * Rtot (QuantityLike): Total radius of the star.
            * Ltot (QuantityLike): Total luminosity of the star.
            * Tc (QuantityLike): Central temperature of the star.
            * X (QuantityLike): Fraction of star mass in H.
            * Y (QuantityLike): Fraction of mass in He.
        """

        # Checking if Mtot has units assigned
        if isinstance(Mtot, u.Quantity):
            self.massInputUnits = Mtot.unit                              # Storing input units
            self.Mtot = (Mtot.to(self.UNIT_SYSTEM["mass"])).value        # Changing to model units
        else:
            self.massInputUnits = u.dimensionless_unscaled               # Storing 'empty' input units
            self.Mtot = Mtot                                             # Preserving model units

        # Checking if Rtot has units assigned
        if isinstance(Rtot, u.Quantity):
            self.radiusInputUnits = Rtot.unit                            # Storing input units
            self.Rtot = (Rtot.to(self.UNIT_SYSTEM["radius"])).value      # Changing to model units
        else:
            self.radiusInputUnits = u.dimensionless_unscaled             # Storing 'empty' input units
            self.Rtot = Rtot                                             # Preserving model units

        # Checking if Ltot has units assigned
        if isinstance(Ltot, u.Quantity):
            self.luminosityInputUnits = Ltot.unit                        # Storing input units
            self.Ltot = (Ltot.to(self.UNIT_SYSTEM["luminosity"])).value  # Changing to model units
        else:
            self.luminosityInputUnits = u.dimensionless_unscaled         # Storing 'empty' input units
            self.Ltot = Ltot                                             # Preserving model units

        # Checking if Tc has units assigned
        if isinstance(Tc, u.Quantity):
            self.temperatureInputUnits = Tc.unit                         # Storing input units
            self.Tc = (Tc.to(self.UNIT_SYSTEM["temperature"])).value     # Changing to model units
        else:
            self.temperatureInputUnits = u.dimensionless_unscaled        # Storing 'empty' input units
            self.Tc = Tc                                                 # Preserving model units

        # Checking if X has units assigned
        if isinstance(X, u.Quantity):
            self.X = X.value
        else:
            self.X = X
        
        # Checking if Y has units assigned
        if isinstance(Y, u.Quantity):
            self.Y = Y.value
        else:
            self.Y = Y

        # Defining other parameters
        self.Z = 1 - self.X - self.Y                        # Metalicity
        self.mu = 1/(2*self.X + 3*self.Y/4 + self.Z/2)      # Mean molecular weight

        # Calculating the model variables
        self.model, self.totalRelativeError = self.__calculate()

    # Defining the __repr__ method
    def __repr__(self):
        return f"{self.__class__.__name__} instance. Atributes: < Mtot={self.Mtot!r}, Rtot={self.Rtot!r}, Ltot={self.Ltot!r}, Tc={self.Tc!r}, X={self.X!r}, Y={self.Y!r}; model units >"
    
    # Defining the __str__ method
    def __str__(self):
        return f"{self.__class__.__name__} instance. Atributes: < Mtot={self.Mtot:.4f}, Rtot={self.Rtot:.4f}, Ltot={self.Ltot:.4f}, Tc={self.Tc:.4f}, X={self.X:.4f}, Y={self.Y:.4f}; model units >"

    # Defining the get method
    def get(self, variable="all", input_units=True, to_csv=False, name=None):
        """
        Returns the requested Star instance data.
        
        ## Arguments:
            * variable (string, default = 'all'):
                If default ('all'), a Data Frame object is returned containing the calculated values of the variables. 
                For queries on specific variables you must enter one of the following strings: 'r', 'P', 'T', 'l', 'm', 'rho', 'epsilon' or 'kappa'.
            
            * input_units (bool, default = True):
                If True, requested data will be expressed using the same units as those used to initialize the Star instance.
                If False, model internal units will be used to express the requested data.

            * to_csv (bool, default = False):
                If True, requested data will be stored in a csv file in the current directory.
                If False, no csv file will be created.

            * name (string, default = None):
                String containing the name of the csv file, if created.
        """

        # Using model data
        modelData = self.model.copy(deep=True)

        # Unsing input units as output 
        if input_units:
            if self.massInputUnits != u.dimensionless_unscaled:
                modelData["m"] = modelData["m"].apply(lambda element: (element*self.UNIT_SYSTEM["mass"]).to(self.massInputUnits).value)
            if self.radiusInputUnits != u.dimensionless_unscaled:
                modelData["r"] = modelData["r"].apply(lambda element: (element*self.UNIT_SYSTEM["radius"]).to(self.radiusInputUnits).value)
            if self.luminosityInputUnits != u.dimensionless_unscaled:
                modelData["l"] = modelData["l"].apply(lambda element: (element*self.UNIT_SYSTEM["luminosity"]).to(self.luminosityInputUnits).value)
            if self.temperatureInputUnits != u.dimensionless_unscaled:
                modelData["T"] = modelData["T"].apply(lambda element: (element*self.UNIT_SYSTEM["temperature"]).to(self.temperatureInputUnits).value)

        # Checking if csv file has a name assigned
        if name is None:
            name = "stellar-interior-numerical-model"

        # Returning all variables
        if variable == "all":
            if to_csv:
                modelData.to_csv(f"{name}.csv")
            return modelData
        # Returning requested variable
        else:
            if to_csv:
                modelData[variable].to_csv(f"{name}.csv")
            return modelData[variable]

    # Defining the parameters method
    def parameters(self):
        """
        Returns Star instance atributes using model units as a list following the order: [Mtot, Rtot, Ltot, Tc, X, Y]
        """
        return [self.Mtot, self.Rtot, self.Ltot, self.Tc, self.X, self.Y]
    
    # Defining the redefining method
    def redefine(self, Mtot=None, Rtot=None, Ltot=None, Tc=None, X=None, Y=None):
        """
        Redefines Star instance atributes. Only given arguments will update.

        ## Arguments:
            * Mtot (QuantityLike, default = None): Total mass of the star.
            * Rtot (QuantityLike, default = None): Total radius of the star.
            * Ltot (QuantityLike, default = None): Total luminosity of the star.
            * Tc (QuantityLike, default = None): Central temperature of the star.
            * X (QuantityLike, default = None): Fraction of star mass in H.
            * Y (QuantityLike, default = None): Fraction of mass in He.
        """

        if Mtot is not None:
            # Checking if Mtot has units assigned
            if isinstance(Mtot, u.Quantity):
                self.massInputUnits = Mtot.unit                              # Storing input units
                self.Mtot = (Mtot.to(self.UNIT_SYSTEM["mass"])).value        # Changing to model units
            else:
                self.massInputUnits = u.dimensionless_unscaled               # Storing 'empty' input units
                self.Mtot = Mtot                                             # Preserving model units

        if Rtot is not None:
            # Checking if Rtot has units assigned
            if isinstance(Rtot, u.Quantity):
                self.radiusInputUnits = Rtot.unit                            # Storing input units
                self.Rtot = (Rtot.to(self.UNIT_SYSTEM["radius"])).value      # Changing to model units
            else:
                self.radiusInputUnits = u.dimensionless_unscaled             # Storing 'empty' input units
                self.Rtot = Rtot                                             # Preserving model units

        if Ltot is not None:
            # Checking if Ltot has units assigned
            if isinstance(Ltot, u.Quantity):
                self.luminosityInputUnits = Ltot.unit                        # Storing input units
                self.Ltot = (Ltot.to(self.UNIT_SYSTEM["luminosity"])).value  # Changing to model units
            else:
                self.luminosityInputUnits = u.dimensionless_unscaled         # Storing 'empty' input units
                self.Ltot = Ltot                                             # Preserving model units

        if Tc is not None:
            # Checking if Tc has units assigned
            if isinstance(Tc, u.Quantity):
                self.temperatureInputUnits = Tc.unit                         # Storing input units
                self.Tc = (Tc.to(self.UNIT_SYSTEM["temperature"])).value     # Changing to model units
            else:
                self.temperatureInputUnits = u.dimensionless_unscaled        # Storing 'empty' input units
                self.Tc = Tc                                                 # Preserving model units

        if X is not None:
            # Checking if X has units assigned
            if isinstance(X, u.Quantity):
                self.X = X.value
            else:
                self.X = X

        if Y is not None:
            # Checking if Y has units assigned
            if isinstance(Y, u.Quantity):
                self.Y = Y.value
            else:
                self.Y = Y
        
        # Calculating the model variables
        self.model, self.totalRelativeError = self.__calculate()

    # Defining the error method
    def error(self):
        """
        Returns the percentage of total relative error of the numerical calculation of the star-interior model.
        """
        return self.totalRelativeError

    # Defining the plot method for star variables
    def visualize(self, x_axis="r", which=["P", "T", "l", "m", "rho"], merge=False, normalize=True, figsize=(8, 6)):
        """
        Graphical representation of the calculated variables throughout the star.

        ## Arguments:
            * x_axis (string, default = 'r'): 
                String to select the independent variable of the plot from the following: 
                'r', 'P', 'T', 'l', 'm', 'rho', 'epsilon' and 'kappa'.

            * which (array-like, default = ['P', 'T', 'l', 'm', 'rho']): 
                Array-like containing the dependent variables desirable to plot in string format.
                Supports the same values as x_axis: 'r', 'P', 'T', 'l', 'm', 'rho', 'epsilon' and 'kappa'.

            * merge (bool, default = False):
                If True, all variables specified in 'which' are graphed in the same figure.
                If False, all variables specified in 'which' are graphed in different figures.

            * normalize (bool, default = True):
                If True, all plots will be graphed using normalized units.
                If False, all plots will be graphed using a mix between cgs and solar units.

            * figsize (two-dimensional array-like, default = (8, 6)):
                Two-dimensional array-like for a better customization on the figures size.
        """

        # Changing font
        plt.rcParams["font.family"] = "serif"

        # Using calculated values
        modelData = self.get(variable="all", input_units=False)

        # Normalized units are used
        if normalize:
            # Redefining model units
            modelData["r"] = modelData["r"] / max(modelData["r"])
            modelData["P"] = modelData["P"] / max(modelData["P"])
            modelData["T"] = modelData["T"] / max(modelData["T"])
            modelData["m"] = modelData["m"] / max(modelData["m"])
            modelData["l"] = modelData["l"] / max(modelData["l"])
            modelData["rho"] = modelData["rho"] / max(modelData["rho"])
            modelData["kappa"] = modelData["kappa"] / max(modelData["kappa"])
            modelData["epsilon"] = modelData["epsilon"] / max(modelData["epsilon"])

            # Defining titles and labels for each variable 
            plots = pd.DataFrame(data=[("Radius", "r / R$_{\\text{total}}$"),
                                    ("Pressure", "P / P$_\\text{c}$"), 
                                    ("Temperature", "T / T$_\\text{c}$"),
                                    ("Luminosity", "$\\ell$ / L$_{\\text{total}}$"),
                                    ("Mass", "m / M$_{\\text{total}}$"),
                                    ("Density", "$\\rho$ / $\\rho_\\text{c}$"),
                                    ("Energy generation rate", "$\\varepsilon$ / $\\varepsilon_{\\text{max}}$"),
                                    ("Opacity", "$\\kappa$ / $\\kappa_{\\text{max}}$")],
                                columns=["title", "label"],
                                index=["r", "P", "T", "l", "m", "rho", "epsilon", "kappa"])

        # Solar units are used
        else:
            # Defining sun parameters
            Msun = 1.9884    # 1e33 g
            Rsun = 6.957     # 1e10 cm
            Lsun = 3.828     # 1e33 erg s^-1
            # Redefining model units
            modelData["r"] = modelData["r"] / Rsun  # cm
            modelData["P"] = modelData["P"] * 1e15  # dyn cm^-2
            modelData["T"] = modelData["T"] * 1e7   # K
            modelData["m"] = modelData["m"] / Msun  # g
            modelData["l"] = modelData["l"] / Lsun  # erg s^-1

            # Defining titles and labels for each variable 
            plots = pd.DataFrame(data=[("Radius", "r / R$_{\\odot}$"),
                                    ("Pressure", "P / dyn cm$^{-2}$"), 
                                    ("Temperature", "T / K"),
                                    ("Luminosity", "$\\ell$ / L$_{\\odot}$"),
                                    ("Mass", "m / M$_{\\odot}$"),
                                    ("Density", "$\\rho$ / g cm$^{-3}$"),
                                    ("Energy generation rate", "$\\varepsilon$ / erg$\\,\\,$g$^{-1}\\,$s$^{-1}$"),
                                    ("Opacity", "$\\kappa$ / cm$^2$ g$^{-1}$")],
                                columns=["title", "label"],
                                index=["r", "P", "T", "l", "m", "rho", "epsilon", "kappa"])

        # Calculating the x value for which the transition to the convective zone occurs
        transition = modelData[modelData["fase"] == "CONVEC"].iloc[0][x_axis] 

        # All curves in the same figure (normalized plots)
        if merge:
            plt.figure(figsize=figsize)
            for variable in which:
                plt.plot(modelData[x_axis], modelData[variable]/modelData[variable].max(), label=plots.loc[variable]["title"], linewidth=2.5, alpha=0.8)

            # Customizing the plot
            plt.title("Stellar-interior model", fontsize=20, weight="bold")     # Title
            plt.xlabel(plots.loc[x_axis]["label"], fontsize=18)                 # x axis label
            # plt.ylabel("Normalized magnitude", fontsize=14)                   # y axis label
            plt.xlim((min(modelData[x_axis]), max(modelData[x_axis])))        # x limits
            plt.ylim((-0.05, 1.05))                                             # y limits
            plt.tick_params(axis="both", labelsize=16)                          # Numbering size
            plt.gca().tick_params(direction="in", which="major", length=8)      # Major ticks size and orientation
            plt.gca().tick_params(direction="in", which="minor", length=3)      # Minor ticks size and orientation
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (x axis)
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (y axis)
            
            if x_axis in ["r", "l", "m", "kappa"]:
                plt.axvspan(-1, transition, color="gray", alpha=0.5, label="Convective zone")                               # Marking the convective zone of the star
                plt.axvspan(transition, modelData[x_axis].iloc[0]*1.5, color="gray", alpha=0.1, label="Radiative zone")    # Marking the radiative zone of the star
            else:
                plt.axvspan(-1, transition, color="gray", alpha=0.1, label="Radiative zone")                                                # Marking the convective zone of the star
                plt.axvspan(transition, modelData[x_axis].iloc[len(modelData)-1]*1.5, color="gray", alpha=0.5, label="Convective zone")   # Marking the radiative zone of the star

            plt.legend(fontsize=15)                                                                                     # Leyend
            plt.grid(which="major", linestyle="-", color = "black", linewidth=0.5, alpha=0.4, visible=True)             # Major grid
            # plt.grid(which="minor", linestyle=":", linewidth=0.5, visible=True, alpha=0.5)                            # Minor grid

        # Curves in different figures
        else:
            for variable in which:
                plt.figure(figsize=figsize)
                plt.plot(modelData[x_axis], modelData[variable], color="red", linewidth=2.5, alpha=0.8)

                # Customizing the plots for each figure
                plt.title(plots.loc[variable]["title"], fontsize=20, weight="bold") # Title
                plt.xlabel(plots.loc[x_axis]["label"], fontsize=18)                 # x axis label
                plt.ylabel(plots.loc[variable]["label"], fontsize=18)               # y axis label
                plt.xlim((min(modelData[x_axis]), max(modelData[x_axis])))        # x limits
                plt.ylim((min(modelData[variable])-0.05*max(modelData[variable]), max(modelData[variable])*1.05))  # y limits
                plt.tick_params(axis="both", labelsize=16)                          # Numbering size
                plt.gca().tick_params(direction="in", which="major", length=8)      # Major ticks size and orientation
                plt.gca().tick_params(direction="in", which="minor", length=3)      # Minor ticks size and orientation
                plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (x axis)
                plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))             # Minor ticks (y axis)

                if x_axis in ["r", "l", "m", "kappa"]:
                    plt.axvspan(-1, transition, color="gray", alpha=0.5, label="Convective zone")                               # Marking the convective zone of the star
                    plt.axvspan(transition, modelData[x_axis].iloc[0]*1.5, color="gray", alpha=0.1, label="Radiative zone")    # Marking the radiative zone of the star
                else:
                    plt.axvspan(-1, transition, color="gray", alpha=0.1, label="Radiative zone")                                                # Marking the convective zone of the star
                    plt.axvspan(transition, modelData[x_axis].iloc[len(modelData)-1]*1.5, color="gray", alpha=0.5, label="Convective zone")   # Marking the radiative zone of the star

                plt.legend(fontsize=15)                                                                                     # Leyend
                plt.grid(which="major", linestyle="-", color = "black", linewidth=0.5, alpha=0.4, visible=True)             # Major grid
                # plt.grid(which="minor", linestyle=":", linewidth=0.5, visible=True)                                       # Minor grid

        plt.show()

    # Defining the Temperature-Density Diagram method
    def TDD(self):
        """
        Graphical representation of the star variables in the Temperature-Density Diagram. 
        Several regions are distinguished depending on the dominant pressure.
        I: ideal gas. II: degeneracy. III: relativistic degeneracy. IV: radiation pressure.
        """

        # Defining fundamental constants and star constants
        R = 8.31447 * 10**7             # Universal gas constant
        a = 7.56578e-15                 # Radiation density constant
        mu_e = 2/(1+self.X)             # Mean molecular electron weight

        # Selecting T and rho
        T = np.log10(self.get("T", input_units=False)*1e7)       # Star log10 temperature in K
        rho = np.log10(self.get("rho", input_units=False))       # Star log10 density in g/cm^3

        # Changing font and figure size. Setting the plot limits
        plt.rcParams["font.family"] = "serif"
        plt.figure(figsize=(9,7))
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
        # plt.text(0.07, 0.30, "Ideal gas", transform=plt.gca().transAxes, fontsize=17, color="white").set_path_effects([path_effects.withStroke(linewidth=2, foreground="black")])
        # plt.text(0.10, 0.65, "Degeneracy", transform=plt.gca().transAxes, fontsize=18, color="white").set_path_effects([path_effects.withStroke(linewidth=2, foreground="black")])
        # plt.text(0.12, 0.885, "Relativistic degeneracy", transform=plt.gca().transAxes, fontsize=17, color="white").set_path_effects([path_effects.withStroke(linewidth=2, foreground="black")])
        # plt.text(0.60, 0.2, "Radiation pressure", transform=plt.gca().transAxes, fontsize=19, color="white").set_path_effects([path_effects.withStroke(linewidth=2, foreground="black")])
        plt.text(0.07, 0.30, "Ideal gas", transform=plt.gca().transAxes, fontsize=20, color="black").set_path_effects([path_effects.withStroke(linewidth=2.5, foreground="white")])
        plt.text(0.10, 0.65, "Degeneracy", transform=plt.gca().transAxes, fontsize=20, color="white").set_path_effects([path_effects.withStroke(linewidth=2.5, foreground="black")])
        plt.text(0.12, 0.885, "Relativistic degeneracy", transform=plt.gca().transAxes, fontsize=20, color="white").set_path_effects([path_effects.withStroke(linewidth=2.5, foreground="black")])
        plt.text(0.55, 0.2, "Radiation pressure", transform=plt.gca().transAxes, fontsize=20, color="black").set_path_effects([path_effects.withStroke(linewidth=2.5, foreground="white")])

        # Graphing the star variables in the diagram
        plt.plot(T, rho, color="black", linewidth=3.5, zorder=1)
        plt.plot(T, rho, color="orange", linewidth=2.5, zorder=1)
        plt.scatter(T.iloc[len(T)-1], rho.iloc[len(rho)-1], s=280, color="orange", marker="*", edgecolors="black", linewidths=0.6, label="Star center", zorder=2)
        plt.scatter(T.iloc[0], rho.iloc[0], s=120, color="orange", marker="s", edgecolors="black", linewidths=0.6, label="Star surface", zorder=2)

        # Adding title and labels. Setting the ticks parameters
        plt.title("Temperature-Density Diagram", fontsize=20, weight="bold")# Title
        plt.xlabel("Log [T(K)]", fontsize=18)                               # x axis label
        plt.ylabel("Log [$\\rho$(g cm$^{-3})$]", labelpad=-7.5, fontsize=18)# y axis label
        plt.tick_params(axis="both", labelsize=16)                          # Numbering size
        plt.gca().tick_params(direction="in", which="major", length=8)      # Major ticks size and orientation
        plt.gca().tick_params(direction="in", which="minor", length=3)      # Minor ticks size and orientation
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))             # Setting minor ticks (x axis)
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))             # Setting minor ticks (y axis)
        plt.xlim(xlims)                                                     # x limits
        plt.ylim(ylims)                                                     # y limits
        plt.legend(frameon=True, edgecolor="black",loc="upper right", fontsize=16) # Legend

        plt.show()

    # Defining the Hertzsprung–Russell Diagram method
    def HR(self):
        """
        Graphical representation of the star in the Hertzsprung–Russell Diagram.
        """

        # Loading data
        current_directory = os.path.dirname(__file__)
        data_route = os.path.join(current_directory, "hertzsprung-russell-data.csv")
        stars_data = pd.read_csv(data_route)

        # Customizing the figure
        plt.rcParams["font.family"] = "serif"           # Font
        plt.figure(facecolor="black", figsize=(6, 7))   # Background color
        plt.gca().set_facecolor("black")
        for spine in plt.gca().spines.values():         # Axes color
            spine.set_edgecolor("white")

        # Plotting stars
        plt.scatter(stars_data["log(Teff)"], stars_data["log(L/Lo)"], 
                    s=stars_data["Plot Size"], c=stars_data["Color"], 
                    edgecolors="black", linewidth=0.5)
        
        # Defining Stefan-Boltzmann constant and sun parameters
        Lsun = 3.828        # 1e33 erg s^-1
        sigma = 5.67040e-5  # erg cm^-2 s^-1 K^-4
        # Defining log(L) and log(Teff)
        star_L = np.log10( self.Ltot/Lsun )
        star_Teff = np.log10( ((self.Ltot*1e33) / (4*np.pi*((self.Rtot*1e10)**2)*sigma ) )**0.25 )
        # Plotting model star
        plt.scatter(star_Teff, star_L, s=700, c="gold", alpha=0.2)
        plt.scatter(star_Teff, star_L, marker="*", s=400, c="gold", label="Star", edgecolors="black", linewidth=0.5)

        # Customizing the graph
        plt.title("Hertzsprung-Russell Diagram", color="white", weight="bold", fontsize=18)         # Title
        plt.xlabel("T$_\\text{eff}$ / K ", color="white", fontsize=16)                                              # x label
        plt.ylabel("Log (L / L$_{\\!\\odot}$)", color="white", fontsize=16)                           # y label
        plt.gca().tick_params(which="major", length=9, direction="in", colors="white", labelsize=13)  # Major ticks length and orientation
        plt.tick_params(axis="x", which="minor", length=4, direction="in", color="white")             # Minor x ticks
        plt.tick_params(axis="y", which="minor", length=4, direction="in", color="white")             # Minor y ticks
        plt.xticks(np.log10(np.array([31000, 9850, 5000, 2500, 1000]), dtype=float),                  # Major x ticks
            labels=["30,000", "10,000", "5,000", "2,500", "1,000"], color="white")

        # Showing the plot
        plt.gca().invert_xaxis()
        plt.minorticks_on()
        plt.legend(fontsize=16)

        # Printing the effective temperature
        print(f"Effective temperature: Teff = {10**star_Teff:.0f} K.")

        plt.show()

    # Defining the function to calculate 
    def __calculate(self):
        """
        Performs the calculation of the stellar-interior numerical model. 
        Returns both the calculated variables and the total relative error of the calculation.
        """

        ################################################################################
        ################################################################################
        ######################## ----- Defining equations ----- ########################
        ################################################################################
        ################################################################################

        # Density
        def perfect_gas_rho(self, P, T):
            """
            Equation (6)
            Both pressure and temperature are expected in the modified unit system.
            """
            R = 8.31447 * 10**7
            return (self.mu/R)*((P*1e15)/(T*1e7))
        
        # Polytrope constant
        def K_pol(self, P, T):
            return P/(T**2.5)
        
        def polytrope(self, K, T):
            return K*(T**2.5)
        
        # Opacity
        def opacity(self, rho, T):
            """
            Equation (8)
            Both density and temperature are expected in the modified unit system.
            """
            g_bf = 3.162
            return 4.34e25 * g_bf * self.Z * (1+self.X) * rho / ((T*1e7)**3.5)

        # Calculating the most efficient energy generation cycle for a given pressure and temperature
        def calculate_optimal_cycle(self, P, T):
            """
            For a given pressure and temperature it retruns the values of epsilon, epsilon1, X1, X2 and nu using the most efficient energy generation cycle.
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
            epsilon_PP = epsilon1_PP * (self.X*self.X) * perfect_gas_rho(self, P, T) * (T*10)**nu_PP
            epsilon_CN = epsilon1_CN * (self.X*self.Z/3) * perfect_gas_rho(self, P, T) * (T*10)**nu_CN
            
            # If there is no energy generation
            if epsilon_PP == 0.0 and epsilon_CN == 0.0:
                return 0.0, 0.0, 0.0, 0.0, 0.0, "--"
            # If the energy generation rate given by PP cycle is greater
            if epsilon_PP >= epsilon_CN:
                return epsilon_PP, epsilon1_PP, self.X, self.X, nu_PP, "PP"
            # If the energy generation rate given by CN cycle is greater
            else:
                return epsilon_CN, epsilon1_CN, self.X, self.Z/3, nu_CN, "CN"


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

        def dPdr_rad(self, r, P, T, m):
            """
            Equation (19)
            """
            Cp = 8.084*self.mu
            return - Cp * P * m / (T * (r**2))

        def dLdr_rad(self, r, P, T):
            """
            Equation (20)
            """
            epsilon, epsilon1, X1, X2, nu, cycle = calculate_optimal_cycle(self, P, T)
            Cl = 0.01845*epsilon1*X1*X2*(10**nu)*self.mu**2
            return Cl * ((P*r)**2) * (T**(nu-2)), epsilon, cycle

        def dTdr_rad(self, r, P, T, l):
            """
            Equation (21)
            """
            Ct = 0.01679*self.Z*(1+self.X)*self.mu*self.mu
            return -Ct * (l*(P**2)) / ((T**8.5)*(r**2))


        # ------------------------------ Convective case ----------------------------- #

        def dMdr_conv(self, r, T, K):
            """
            Equation (22)
            """
            Cm = 0.01523*self.mu
            return Cm * K*(T**1.5)*r**2

        def dPdr_conv(self, r, T, m, K):
            """
            Equation (23)
            """
            if r == 0.0:
                return 0.0
            else:
                Cp = 8.084*self.mu
                return -Cp * K*(T**1.5)*m / (r**2)

        def dLdr_conv(self, r, P, T, K):
            """
            Equation (24)
            """
            _, epsilon1, X1, X2, nu, cycle = calculate_optimal_cycle(self, P, T)
            Cl = 0.01845*epsilon1*X1*X2*(10**nu)*self.mu**2
            return Cl * (K**2)*(T**(3+nu))*r**2, cycle

        def dTdr_conv(self, r, m):
            """
            Equation (25)
            """
            if r == 0.0:
                return 0.0
            else:
                Ct = 3.234*self.mu
                return -Ct * m / (r**2)


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
            epsilon, epsilon1, X1, X2, nu, cycle = calculate_optimal_cycle(self, P, T)
            return 0.006150*epsilon1*X1*X2*(10**nu)*(self.mu**2)*(K**2)*(self.Tc**(3+nu))*(r**3), epsilon, cycle

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

        def step3(self, shell, derivatives, r, P, T, m, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            and temperature for the next shell and the derivatives of the current and two 
            previous shells it calculates the value of the mass and its derivative for the 
            next shell.
            """

            fm = dMdr_rad(self, r, P, T)             # fm for the i+1 shell
            AM1 = h * (fm - derivatives[2][3])       # AM1 for the i+1 shell

            m_cal = m + h*fm - AM1/2                 # Calculated mass for the i+1 shell

            return m_cal, fm


        # ---------------------------------- Step 4 ---------------------------------- #

        def step4(self, shell, derivatives, r, P, T, m, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            temperature and mass for the next shell and the derivatives of the current and 
            two previous shells it calculates the value of the pressure and its derivative 
            for the next shell.
            """

            fP = dPdr_rad(self, r, P, T, m)          # fP for the i+1 shell
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

        def step6(self, shell, derivatives, r, P, T, l, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure 
            and temperature for the next shell and the derivatives of the current and two 
            previous shells it calculates the value of the luminosity and its derivative 
            for the next shell. It also returns which energy generation cycle produces 
            the energy.
            """

            fl, epsilon, cycle = dLdr_rad(self, r, P, T)                              # fl for the i+1 shell
            AL1 = h * (fl - derivatives[2][2])                               # AL1 for the i+1 shell
            AL2 = h * (fl - 2*derivatives[2][2] + derivatives[1][2])         # AL2 for the i+1 shell
            
            l_cal = l + h*fl - AL1/2  - AL2/12        # Calculated luminosity for the i+1 shell

            return l_cal, fl, epsilon, cycle


        # ---------------------------------- Step 7 ---------------------------------- #

        def step7(self, shell, derivatives, r, P, T, l, h):
            """
            Given the variables of the star for a shell, an estimation of radius, pressure,
            temperature and luminosity for the next shell and the derivatives of the current 
            and two previous shells it calculates the value of the temperature and its 
            derivative for the next shell. 
            """

            fT = dTdr_rad(self, r, P, T, l)          # fT for the i+1 shell
            AT1 = h * (fT - derivatives[2][1])       # AT1 for the i+1 shell
            T = shell[4]                             # T for the i shell

            T_cal = T + h*fT - AT1/2                 # Calculated temperature for the i+1 shell

            return T_cal, fT


        # -------------------------------- Step 7 bis -------------------------------- #

        def step7bis(self, shell, derivatives, r, m, h):
            """
            Given the variables of the star for a shell, an estimation of radius and mass 
            for the next shell and the derivatives of the current and two previous shells it 
            calculates the value of the temperature and its derivative for the next shell. 
            """

            fT = dTdr_conv(self, r, m)               # fT for the i+1 shell
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

        ################################################################################
        # ----------------------- First three shells (surface) ----------------------- #
        ################################################################################

        Rini = 0.9*self.Rtot            # Defining the initial radius
        shellBreaks = 100               # Setting the shell partition of Rini
        innerShells = shellBreaks + 1   # Setting the number of shells from Rini
        h = - Rini/shellBreaks          # Defining the radius step between shells (negative)
        r = Rini                        # Initializing the radius variable

        outer_data = []                 # Initializing the list that will store the data (surface integration)
        outer_derivatives = []          # Initializing the list that will store derivatives (surface integration)

        for i in range(3):

            # Calculating and storing the values for the first three shells
            T = initial_surface_T(self, r)
            P = initial_surface_P(self, r, T)
            m = self.Mtot
            l = self.Ltot
            rho = perfect_gas_rho(self, P, T)
            kappa = opacity(self, rho, T)
            outer_data.append(["--", "START", r, P, T, l, m, rho, 0.0, kappa])

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_rad(self, r, P, T, l)
            fP = dPdr_rad(self, r, P, T, m)
            fm = 0.0
            fl = 0.0
            outer_derivatives.append([fP, fT, fl, fm])

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
                    m_cal, fm = step3(self, shell, derivatives, r, P_est, T_est, shell[6], h)
                    # Performing the fourth step
                    P_cal, fP = step4(self, shell, derivatives, r, P_est, T_est, m_cal, h)
                    # Performing the fifth step
                    if stepX(self, P_cal, P_est):
                        break   # Third loop break
                    else:
                        P_est = P_cal

                # Performing the sixth step
                l_cal, fl, epsilon, cycle = step6(self, shell, derivatives, r, P_cal, T_est, shell[5], h)
                # Performing the seventh step
                T_cal, fT = step7(self, shell, derivatives, r, P_cal, T_est, l_cal, h)
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
                # Computing density and opacity
                rho = perfect_gas_rho(self, P_cal, T_cal)
                kappa = opacity(self, rho, T_cal)
                # Storing values and derivatives
                outer_data.append([cycle, fase, r, P_cal, T_cal, l_cal, m_cal, rho, epsilon, kappa])
                outer_derivatives.append([fP, fT, fl, fm])

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
            m = initial_center_M(self, r, self.Tc, K)
            l, epsilon, _ = initial_center_L(self, r, P, self.Tc, K)
            rho = perfect_gas_rho(self, P, T)
            kappa = opacity(self, rho, T)
            inside_data.append(["--", "CENTER", r, P, T, l, m, rho, epsilon, kappa])

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_conv(self, r, m)
            fP = dPdr_conv(self, r, T, m, K)
            fm = dMdr_conv(self, r, T, K)
            fl, _ = dLdr_conv(self, r, P, T, K)
            inside_derivatives.append([fP, fT, fl, fm])

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
                m_cal, fm = step3(self, shell, derivatives, r, P_est, T_est, shell[6], h)
                # Performing the seventh step bis
                T_cal, fT = step7bis(self, shell, derivatives, r, m_cal, h)
                # Performing the eighth step
                if stepX(self, T_cal, T_est):
                    break       # Loop break
                else:
                    T_est = T_cal

            # Calculating pressure using the polytrope constant
            P_cal = polytrope(self, K, T_cal)
            # Calculating the derivative of the pressure for the i shell
            fP = dPdr_conv(self, r, T_cal, m_cal, K)
            # Performing the sixth step
            l_cal, fl, epsilon, cycle = step6(self, shell, derivatives, r, P_cal, T_cal, shell[5], h)
            # Computing density and opacity
            rho = perfect_gas_rho(self, P_cal, T_cal)
            kappa = opacity(self, rho, T_cal)
            # Storing values and derivatives
            inside_data.append([cycle, fase, r, P_cal, T_cal, l_cal, m_cal, rho, epsilon, kappa])
            inside_derivatives.append([fP, fT, fl, fm])

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
            m = self.Mtot
            l = self.Ltot
            rho = perfect_gas_rho(self, P, T)
            kappa = opacity(self, rho, T)
            outer_data.insert(0, ["--", "^^^^^^", r, P, T, l, m, rho, 0.0, kappa])

            # Calculating and storing the derivatives for the first three shells
            fT = dTdr_rad(self, r, P, T, l)
            fP = dPdr_rad(self, r, P, T, m)
            fm = 0.0
            fl = 0.0
            outer_derivatives.insert(0, [fP, fT, fl, fm])

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
        model = pd.DataFrame(model, columns=["E", "fase", "r", "P", "T", "l", "m", "rho", "epsilon", "kappa"])
        model.index -= n

        # Customizing the DataFrame options
        pd.set_option("colheader_justify", "center", "display.max_rows", None)
        pd.options.display.float_format = '{:.7f}'.format
        
        return model, totalRelativeError