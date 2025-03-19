<p align="center">
  <img src="images/StellarPyLogo.png" alt="StarPyLogo" width="300px" height="auto">
</p>

---

Python module for solving the stellar-interior equations.

 - Author: **Nuno Cerviño Luridiana**
 - Last update: March, 2025

It provides:

* A Star object to model massive stars.
* Two different functions to optimize modeling.
* A 351 stars data set for Hertzsprung-Russell diagram representation.


# Installation

Using the terminal, clone the repository to your local machine.

```sh
git clone https://github.com/gallati/stellar-interior-numerical-model
```

Enter the cloned directory and install the package.

```sh
cd stellarpy/pip install .
```

Now you can access all StellarPy functionalities!


# Star class

Represents a star with a given mass and chemical composition. To preform the stellar-interior numerical calculation, initial values for radius, luminosity and central temperature are required. How ever, StellarPy provides the functions `error_table` and `find_minimum` to optimize which values better depict the star.

## Initialization

The `Star` object takes the following initial parameters:

* `Rtot` (float, default = 11.5): Total radius of the star.
* `Ltot` (float, default = 70.0): Total luminosity of the star.
* `Tc` (float, default = 2.0): Central temperature of the star.
* `X` (float, default = 0.75): Fraction of star mass in H.
* `Y` (float, default = 0.22): Fraction of mass in He.
* `solar_units` (bool, default = False): Specifies whether solar units are to be used as input.

Once the object is initialized, a numerical estimation for radius, pressure, temperature, mass, luminosity, density, energy generation rate and opacity throughout the star is performed.


## Methods

Several built-in methods are provided for the `Star` object.

* `get`

    Returns the requested Star instance data. Arguments:

    * `variable` (string, default = 'all'):
        If default ('all'), a Data Frame object is returned containing the calculated values of the variables. 
        For queries on specific variables you must enter one of the following strings: 'r', 'P', 'T', 'l', 'm', 'rho', 'epsilon' or 'kappa'.

    * `solar_units` (bool, default = False):
        If True, all data will be given using solar units.
        If False, all data will be given using the model units.


* `parameters`
    
    Returns Star instance atributes as a list following the order: [Mtot, Rtot, Ltot, Tc, X, Y]


* `redefine`

    Redefines Star instance atributes. Arguments:

    * `Mtot` (float, default = 5.0): Total mass of the star.
    * `Rtot` (float, default = 11.5): Total radius of the star.
    * `Ltot` (float, default = 70.0): Total luminosity of the star.
    * `Tc` (float, default = 2.0): Central temperature of the star.
    * `X` (float, default = 0.75): Fraction of star mass in H.
    * `Y` (float, default = 0.22): Fraction of mass in He.


* `error`

    Returns the percentage of total relative error of the numerical calculation of the star-interior model.


* `visualize`
    
    Graphical representation of the calculated variables throughout the star. Arguments:
        
    * `x_axis` (string, default = 'r'): 
        String to select the independent variable of the plot from the following: 
        'r', 'P', 'T', 'l', 'm', 'rho', 'epsilon' and 'kappa'.

    * `which` (array-like, default = ['P', 'T', 'l', 'm', 'rho']): 
        Array-like containing the dependent variables desirable to plot in string format.
        Supports the same values as x_axis: 'r', 'P', 'T', 'l', 'm', 'rho', 'epsilon' and 'kappa'.

    * `merge` (bool, default = False):
        If True, all variables specified in 'which' are graphed in the same figure.
        If False, all variables specified in 'which' are graphed in different figures.

    * `normalize` (bool, default = True):
        If True, all plots will be graphed using normalized units.
        If False, all plots will be graphed using a mix between cgs and solar units.

    * `figsize` (two-dimensional array-like, default = (8, 6)):
        Two-dimensional array-like for a better customization on the figures size.


* `TDD`

    Graphical representation of the star variables in the Temperature-Density Diagram. Several regions are distinguished depending on the dominant pressure. I: ideal gas. II: degeneracy. III: relativistic degeneracy. IV: radiation pressure.


* `HR`

    Graphical representation of the star in the Hertzsprung–Russell Diagram.


## Units

In order to properly estimate the variables of the star, the unit system adopted for internal calculations of the model varies with respect to CGS. However, both input and output values of the model can be converted to solar units.

    radius (r)                         ->   1e10 cm
    pressure (P)                       ->   1e15 dyn cm^-2
    temperature (T)                    ->   1e7 K
    mass (M)                           ->   1e33 g
    luminosity (L)                     ->   1e33 erg s^-1
    density (rho)                      ->   1 g cm^-3
    energy generation rate (epsilon)   ->   1 erg g^-1 s^-1
    opacity (kappa)                    ->   1 cm^2 g^-1


# Optimization functions

* `error_table`

    Table containing the total relative error for total luminosity and total radius variations. Given a Star object, total relative error for variations of Ltot and Rtot is computed. Arguments:

    * `star` (Star): 
        Star object for which the minimum must be found.

    * `n` (float):
        Size of the maximum variation. The output table length will be (2*n+1).

    * `dR` (float, default = 0.5): 
        Total radius variation.

    * `dL` (float, default = 5.0): 
        Total luminosity variation.
    
    * `numbering` (bool, default = False):
        Enables table numbering.


* `find_minimum`
    
    Total relative error minimum finder for Star objects. Arguments:

    * `star` (Star):

        Star object for which the minimum must be found.

    * `x0` (list, default = None): 

        List containing specific initial parameters required for the minimum search, listed as [Rtot, Ltot, Tc] in model units. If no list is provided, current parameters of the Star object will be used.

# Usage example

Let us compute the value of the stellar-interior variables for a star with the follwing parameters:

$$M = 2.51 M_\odot \quad R = 1.59 R_\odot \quad L = 19.76 L_\odot \quad T_c = 1.9554\cdot10^7 \text{K}$$

To do so, the `Model` object is initialized as shown.

```sh
model = Model(Mtot = 5.0, Rtot = 11.0570, Ltot = 75.9213, Tc = 1.9554)
```

The `get()` method is required to access, for example, the effective temperature $T_{\text{eff}}$ of the star.

```sh
T = model.get(variable = 'T')
Teff = T.iloc[0]
print(Teff)
```
```sh
[Output] 0.0005125238763055535
```

Which leads us to the result $T_{\text{eff}}=5125 \text{K}$

For the visualization of variables throughout the star the `visualize()` method is needed.

```sh
model.visualize(x_axis = 'r', which = ["P", "T", "L", "M", "rho"], merge = True)
```

![example](images/example.png) 


License
----

**Free Software!** 
For the benefit of everyone.