<p align="center">
  <img src="images/StellarPyLogo.png" alt="StarPyLogo" width="300px" height="auto">
</p>

---

Python module for solving the stellar-interior equations.

 - Author: **Nuno Cerviño Luridiana**
 - Last update: April, 2025

It provides:

* A Star object to model massive stars.
* Two different functions to optimize modeling.
* A 351 stars data set for Hertzsprung-Russell diagram representation retrived from [Pecaut M. \& Mamajek J](https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt) and Kaggle repository [Chen S.](https://www.kaggle.com/code/salmanhiro/hertzsprung-russell-diagram/input).


# Installation

Using the terminal, clone this repository to your local machine.

```sh
git clone https://github.com/gallati/stellarpy
```

Enter the cloned directory and install the package.

```sh
cd stellarpy/pip install .
```

Now you can access all StellarPy functionalities!


# Star class

Represents a star with a given mass and chemical composition. To preform the stellar-interior numerical calculation, initial values for radius, luminosity and central temperature are required. How ever, StellarPy provides the functions `error_table` and `find_minimum` to optimize which values better depict the star.

## Atributes

The `Star` object takes the following initial parameters:

* `Mtot` (QuantityLike): Total mass of the star.
* `Rtot` (QuantityLike): Total radius of the star.
* `Ltot` (QuantityLike): Total luminosity of the star.
* `Tc` (QuantityLike): Central temperature of the star.
* `X` (QuantityLike): Fraction of star mass in H.
* `Y` (QuantityLike): Fraction of mass in He.

Once the object is initialized, a numerical estimation for radius, pressure, temperature, mass, luminosity, density, energy generation rate and opacity throughout the star is performed.


## Methods

Several built-in methods are provided for the `Star` object.

* `get`

    Returns the requested Star instance data. Arguments:

    * `variable` (string, default = 'all'):
        If default ('all'), a Data Frame object is returned containing the calculated values of the variables. 
        For queries on specific variables you must enter one of the following strings: 'r', 'P', 'T', 'l', 'm', 'rho', 'epsilon' or 'kappa'.

    * `input_units` (bool, default = True):
        If True, requested data will be expressed using the same units as those used to initialize the Star instance. If False, model internal units will be used to express the requested data.

    * `to_csv` (bool, default = False):
        If True, requested data will be stored in a csv file in the current directory.
        If False, no csv file will be created.

    * `name` (string, default = None):
        String containing the name of the csv file, if created.


* `parameters`
    
    Returns Star instance atributes as a list following the order: [Mtot, Rtot, Ltot, Tc, X, Y]


* `redefine`

    Redefines Star instance atributes. Arguments:

    * `Mtot` (QuantityLike, default = None): Total mass of the star.
    * `Rtot` (QuantityLike, default = None): Total radius of the star.
    * `Ltot` (QuantityLike, default = None): Total luminosity of the star.
    * `Tc` (QuantityLike, default = None): Central temperature of the star.
    * `X` (QuantityLike, default = None): Fraction of star mass in H.
    * `Y` (QuantityLike, default = None): Fraction of mass in He.


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

In order to properly estimate the variables of the star, the unit system adopted for internal calculations of the model varies with respect to CGS. However, both input and output values of the model can be expressed in any unit system using `Quantity` objects from astropy.

    radius (r)                         ->   1e10 cm
    pressure (P)                       ->   1e15 dyn cm^-2
    temperature (T)                    ->   1e7 K
    mass (m)                           ->   1e33 g
    luminosity (l)                     ->   1e33 erg s^-1
    density (rho)                      ->   1 g cm^-3
    energy generation rate (epsilon)   ->   1 erg g^-1 s^-1
    opacity (kappa)                    ->   1 cm^2 g^-1


# Optimization functions

StellarPy provides two function to optimize which values better depict the star.

* `error_table`

    Table containing the total relative error for total luminosity and total radius variations. Given a Star object, total relative error for variations of Ltot and Rtot is computed. Arguments:

    * `star` (Star): 
        Star object for which the minimum must be found.

    * `n` (int):
        Size of the maximum variation. The output table length will be (2*n+1).

    * `dR` (float): 
        Total radius variation.

    * `dL` (float): 
        Total luminosity variation.
    
    * `numbering` (bool, default = False):
        Enables table numbering.


* `find_minimum`
    
    Total relative error minimum finder for Star objects. Arguments:

    * `star` (Star):

        Star object for which the minimum must be found.

    * `x0` (list, default = None): 

        List containing specific initial parameters required for the minimum search, listed as [Rtot, Ltot, Tc] in model units. If no list is provided, current parameters of the Star object will be used.

    Returns:

    Optimized parameters and total relative error as a list following the order `[Rtot, Ltot, Tc, error]`.

# Usage example

StellarPy provides a Jupyter Notebook file in which an example is followed through. Here are some of the results achieved in that file.

Let us build a model for a star with the following mass and chemical composition:

$$M = 5.0\times10^{33}\text{g} \quad X=0.80 \quad Y=0.16$$

As initial parameters for the rest of the magnitudes we will consider as follows:

$$R = 11.5\times10^{10}\text{cm} \quad L = 40.0\times10^{33}\text{erg}\,\text{s}^{-1} \quad T_c = 1.5\times10^7\text{K}$$

This way, the `Star` object is initialized as shown.

```sh
star = Star(Mtot=5.0, Rtot=11.5, Ltot=40.0, Tc=1.5, X=0.80, Y=0.16)
```

In order to quantify how great our initial choice of parameters is, let us invoke the `error` method.

```sh
star.error()
```
```sh
[Output]: 103.03429829195771
```

This result can be improved significatibly optimizing our model by using the function `find_minimum`.

```sh
Rmin, Lmin, Tmin, error_min = find_minimum(star=star, x0=None)
star.redefine(Rtot=Rmin, Ltot=Lmin, Tc=Tmin)
```
```sh
[Output]: 
Minimum found at (model units): 
   Rtot = 11.2518
   Ltot = 42.5688
   Tc   = 1.8569
Error: 0.0166 %
```
Now that the star is well defined, let's plot some cool graphs! 

## Variables throughout the star

```sh
star.visualize(merge=True, figsize=(10, 6))
```

<p align="center">
  <img src="images/usage_example1.png" alt="StarPyLogo" width="600px" height="auto">
</p>

## Energy generation rate and opacity

```sh
star.visualize(which=["kappa", "epsilon"], normalize=False)
```

<div align="center">
  <img src="images/usage_example2.png" width="45%">
  <img src="images/usage_example3.png" width="45%">
</div>


## Teperature-Density Diagram

```sh
star.TDD()
```

<p align="center">
  <img src="images/usage_example4.png" alt="StarPyLogo" width="600px" height="auto">
</p>

## Hertzsprung–Russell Diagram

```sh
star.HR()
```
<p align="center">
  <img src="images/usage_example5.png" alt="StarPyLogo" width="450px" height="auto">
</p>

License
----

**Free Software!** 
For the benefit of everyone.