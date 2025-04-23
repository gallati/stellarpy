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


## Data table

Finally, to access star-interior data we invoke the get method.

```sh
star.get()
```


|     | E   | fase   |       r |            P |           T |           l |           m |          rho |     epsilon |    kappa |
|----:|:----|:-------|--------:|-------------:|------------:|------------:|------------:|-------------:|------------:|---------:|
| -11 | --  | ^^^^^^ | 11.4885 |  1.40029e-13 | 0.000475788 | 40          | 5           |  2.03432e-10 | 0           | 270.555  |
| -10 | --  | ^^^^^^ | 11.385  |  2.58774e-09 | 0.00480113  | 40          | 5           |  3.72557e-07 | 0           | 151.801  |
|  -9 | --  | ^^^^^^ | 11.2815 |  4.11604e-08 | 0.00920585  | 40          | 5           |  3.09052e-06 | 0           | 129.001  |
|  -8 | --  | ^^^^^^ | 11.178  |  2.22439e-07 | 0.0136921   | 40          | 5           |  1.12294e-05 | 0           | 116.813  |
|  -7 | --  | ^^^^^^ | 11.0745 |  7.56512e-07 | 0.0182623   | 40          | 5           |  2.86337e-05 | 0           | 108.698  |
|  -6 | --  | ^^^^^^ | 10.971  |  1.98614e-06 | 0.0229186   | 40          | 5           |  5.99015e-05 | 0           | 102.698  |
|  -5 | --  | ^^^^^^ | 10.8675 |  4.41902e-06 | 0.0276637   | 40          | 5           |  0.000110416 | 0           |  97.9787 |
|  -4 | --  | ^^^^^^ | 10.764  |  8.76421e-06 | 0.0325      | 40          | 5           |  0.0001864   | 0           |  94.1106 |
|  -3 | --  | ^^^^^^ | 10.6605 |  1.59736e-05 | 0.0374302   | 40          | 5           |  0.000294983 | 0           |  90.8456 |
|  -2 | --  | ^^^^^^ | 10.557  |  2.72896e-05 | 0.0424571   | 40          | 5           |  0.000444287 | 0           |  88.0282 |
|  -1 | --  | ^^^^^^ | 10.4535 |  4.42996e-05 | 0.0475835   | 40          | 5           |  0.000643516 | 0           |  85.555  |
|   0 | --  | START  | 10.35   |  6.89985e-05 | 0.0528125   | 40          | 5           |  0.000903065 | 0           |  83.3538 |
|   1 | --  | START  | 10.2465 |  0.000103861 | 0.0581471   | 40          | 5           |  0.00123464  | 0           |  81.3724 |
|   2 | --  | START  | 10.143  |  0.000151924 | 0.0635905   | 40          | 5           |  0.00165139  | 0           |  79.5722 |
|   3 | --  | RADIAT | 10.0395 |  0.000217319 | 0.0691545   | 40          | 4.99986     |  0.00217217  | 0           |  78.0381 |
|   4 | --  | RADIAT |  9.936  |  0.000304296 | 0.0748457   | 40          | 4.99953     |  0.00281026  | 0           |  76.5503 |
|   5 | --  | RADIAT |  9.8325 |  0.000418256 | 0.0806609   | 40          | 4.99912     |  0.00358423  | 0           |  75.1377 |
|   6 | --  | RADIAT |  9.729  |  0.000565651 | 0.086601    | 40          | 4.99862     |  0.00451484  | 0           |  73.8068 |
|   7 | --  | RADIAT |  9.6255 |  0.000754129 | 0.0926675   | 40          | 4.99799     |  0.00562515  | 0           |  72.556  |
|   8 | --  | RADIAT |  9.522  |  0.00099274  | 0.098864    | 40          | 4.99724     |  0.00694086  | 0           |  71.3778 |
|   9 | --  | RADIAT |  9.4185 |  0.00129214  | 0.105194    | 40          | 4.99633     |  0.0084905   | 0           |  70.2655 |
|  10 | --  | RADIAT |  9.315  |  0.00166485  | 0.111663    | 40          | 4.99525     |  0.0103058   | 0           |  69.2119 |
|  11 | --  | RADIAT |  9.2115 |  0.00212552  | 0.118275    | 40          | 4.99398     |  0.0124219   | 0           |  68.2105 |
|  12 | --  | RADIAT |  9.108  |  0.00269128  | 0.125034    | 40          | 4.99248     |  0.0148781   | 0           |  67.2564 |
|  13 | --  | RADIAT |  9.0045 |  0.00338203  | 0.131946    | 40          | 4.99073     |  0.0177173   | 0           |  66.3435 |
|  14 | --  | RADIAT |  8.901  |  0.00422097  | 0.139015    | 40          | 4.9887      |  0.0209878   | 0           |  65.4689 |
|  15 | --  | RADIAT |  8.7975 |  0.00523506  | 0.146247    | 40          | 4.98635     |  0.0247428   | 0           |  64.6285 |
|  16 | --  | RADIAT |  8.694  |  0.00645551  | 0.153648    | 40          | 4.98366     |  0.0290415   | 0           |  63.8199 |
|  17 | --  | RADIAT |  8.5905 |  0.00791846  | 0.161224    | 40          | 4.98058     |  0.0339491   | 0           |  63.0392 |
|  18 | --  | RADIAT |  8.487  |  0.00966563  | 0.168979    | 40          | 4.97707     |  0.0395378   | 0           |  62.2845 |
|  19 | --  | RADIAT |  8.3835 |  0.0117452   | 0.176921    | 40          | 4.97309     |  0.0458878   | 0           |  61.5531 |
|  20 | --  | RADIAT |  8.28   |  0.0142128   | 0.185056    | 40          | 4.96859     |  0.0530873   | 0           |  60.843  |
|  21 | --  | RADIAT |  8.1765 |  0.0171322   | 0.193391    | 40          | 4.96353     |  0.0612342   | 0           |  60.1523 |
|  22 | --  | RADIAT |  8.073  |  0.0205772   | 0.201932    | 40          | 4.95784     |  0.0704366   | 0           |  59.4793 |
|  23 | --  | RADIAT |  7.9695 |  0.0246323   | 0.210686    | 40          | 4.95146     |  0.0808137   | 0           |  58.8222 |
|  24 | --  | RADIAT |  7.866  |  0.0293947   | 0.219662    | 40          | 4.94435     |  0.0924975   | 0           |  58.1798 |
|  25 | --  | RADIAT |  7.7625 |  0.0349757   | 0.228866    | 40          | 4.93643     |  0.105633    | 0           |  57.5506 |
|  26 | --  | RADIAT |  7.659  |  0.0415031   | 0.238307    | 40          | 4.92763     |  0.120382    | 0           |  56.9334 |
|  27 | --  | RADIAT |  7.5555 |  0.0491232   | 0.247992    | 40          | 4.91788     |  0.136919    | 0           |  56.3268 |
|  28 | --  | RADIAT |  7.452  |  0.0580036   | 0.257931    | 40          | 4.9071      |  0.155442    | 0           |  55.7299 |
|  29 | --  | RADIAT |  7.3485 |  0.0683359   | 0.268132    | 40          | 4.89521     |  0.176164    | 0           |  55.1416 |
|  30 | --  | RADIAT |  7.245  |  0.0803389   | 0.278604    | 40          | 4.88212     |  0.199322    | 0           |  54.5608 |
|  31 | --  | RADIAT |  7.1415 |  0.0942628   | 0.289357    | 40          | 4.86774     |  0.225176    | 0           |  53.9865 |
|  32 | --  | RADIAT |  7.038  |  0.110393    | 0.300401    | 40          | 4.85197     |  0.254013    | 0           |  53.4179 |
|  33 | --  | RADIAT |  6.9345 |  0.129054    | 0.311745    | 40          | 4.83471     |  0.286147    | 0           |  52.8539 |
|  34 | --  | RADIAT |  6.831  |  0.150617    | 0.323399    | 40          | 4.81585     |  0.321923    | 0           |  52.2938 |
|  35 | --  | RADIAT |  6.7275 |  0.175504    | 0.335375    | 40          | 4.79528     |  0.36172     | 0           |  51.7367 |
|  36 | --  | RADIAT |  6.624  |  0.204194    | 0.347684    | 40          | 4.77287     |  0.405953    | 0           |  51.1818 |
|  37 | --  | RADIAT |  6.5205 |  0.237233    | 0.360336    | 40          | 4.74852     |  0.455076    | 0           |  50.6281 |
|  38 | --  | RADIAT |  6.417  |  0.275238    | 0.373343    | 40          | 4.72209     |  0.509584    | 0           |  50.0751 |
|  39 | --  | RADIAT |  6.3135 |  0.31891     | 0.386717    | 40          | 4.69345     |  0.57002     | 0           |  49.5217 |
|  40 | PP  | RADIAT |  6.21   |  0.369043    | 0.400471    | 40          | 4.66247     |  0.636974    | 0.00024305  |  48.9673 |
|  41 | PP  | RADIAT |  6.1065 |  0.426533    | 0.414616    | 40          | 4.62899     |  0.711087    | 0.000334155 |  48.4111 |
|  42 | PP  | RADIAT |  6.003  |  0.492395    | 0.429165    | 40          | 4.59289     |  0.793059    | 0.000458354 |  47.8524 |
|  43 | PP  | RADIAT |  5.8995 |  0.56777     | 0.444131    | 40          | 4.55401     |  0.883644    | 0.000627328 |  47.2902 |
|  44 | PP  | RADIAT |  5.796  |  0.653944    | 0.459527    | 39.9999     | 4.5122      |  0.983661    | 0.000856766 |  46.724  |
|  45 | PP  | RADIAT |  5.6925 |  0.75236     | 0.475366    | 39.9999     | 4.46732     |  1.09399     | 0.0011677   |  46.1529 |
|  46 | PP  | RADIAT |  5.589  |  0.864637    | 0.491661    | 39.9998     | 4.4192      |  1.21558     | 0.00158829  |  45.5763 |
|  47 | PP  | RADIAT |  5.4855 |  0.992587    | 0.508425    | 39.9997     | 4.36771     |  1.34945     | 0.00215611  |  44.9932 |
|  48 | PP  | RADIAT |  5.382  |  1.13823     | 0.525672    | 39.9996     | 4.31269     |  1.49669     | 0.00292128  |  44.4031 |
|  49 | PP  | RADIAT |  5.2785 |  1.30381     | 0.543414    | 39.9994     | 4.254       |  1.65844     | 0.0039504   |  43.8052 |
|  50 | PP  | RADIAT |  5.175  |  1.49184     | 0.561666    | 39.9991     | 4.1915      |  1.83595     | 0.0053319   |  43.1988 |
|  51 | PP  | RADIAT |  5.0715 |  1.70506     | 0.580438    | 39.9987     | 4.12506     |  2.03048     | 0.00718278  |  42.5831 |
|  52 | PP  | RADIAT |  4.968  |  1.94652     | 0.599745    | 39.9981     | 4.05455     |  2.24341     | 0.00965755  |  41.9576 |
|  53 | PP  | RADIAT |  4.8645 |  2.21957     | 0.619598    | 39.9972     | 3.97987     |  2.47614     | 0.0131974   |  41.3214 |
|  54 | PP  | RADIAT |  4.761  |  2.52786     | 0.640008    | 39.996      | 3.90092     |  2.73013     | 0.017111    |  40.6741 |
|  55 | PP  | RADIAT |  4.6575 |  2.87535     | 0.660987    | 39.9944     | 3.81762     |  3.00686     | 0.0221435   |  40.0149 |
|  56 | PP  | RADIAT |  4.554  |  3.26635     | 0.682545    | 39.9921     | 3.72993     |  3.30786     | 0.0286004   |  39.3433 |
|  57 | PP  | RADIAT |  4.4505 |  3.70548     | 0.704691    | 39.9891     | 3.6378      |  3.63464     | 0.0368659   |  38.6589 |
|  58 | PP  | RADIAT |  4.347  |  4.19767     | 0.727434    | 39.985      | 3.54123     |  3.9887      | 0.0474209   |  37.961  |
|  59 | PP  | RADIAT |  4.2435 |  4.74817     | 0.75078     | 39.9795     | 3.44025     |  4.37149     | 0.0608648   |  37.2494 |
|  60 | PP  | RADIAT |  4.14   |  5.36247     | 0.774736    | 39.9722     | 3.33493     |  4.7844      | 0.0779418   |  36.5237 |
|  61 | PP  | RADIAT |  4.0365 |  6.04627     | 0.799305    | 39.9624     | 3.22537     |  5.22867     | 0.0995704   |  35.7834 |
|  62 | PP  | RADIAT |  3.933  |  6.80545     | 0.82449     | 39.9495     | 3.11172     |  5.70542     | 0.12688     |  35.0286 |
|  63 | PP  | RADIAT |  3.8295 |  7.64602     | 0.850292    | 39.9325     | 2.99416     |  6.21561     | 0.161252    |  34.2594 |
|  64 | PP  | RADIAT |  3.726  |  8.57396     | 0.876709    | 39.9102     | 2.87294     |  6.75993     | 0.204363    |  33.4758 |
|  65 | PP  | RADIAT |  3.6225 |  9.59511     | 0.903739    | 39.8813     | 2.74834     |  7.33877     | 0.258239    |  32.6782 |
|  66 | PP  | RADIAT |  3.519  | 10.7151      | 0.931374    | 39.8439     | 2.6207      |  7.9522      | 0.325307    |  31.8667 |
|  67 | PP  | RADIAT |  3.4155 | 11.939       | 0.959606    | 39.7965     | 2.49043     |  8.59986     | 0.398197    |  31.042  |
|  68 | PP  | RADIAT |  3.312  | 13.2715      | 0.988426    | 39.7374     | 2.35798     |  9.28097     | 0.490943    |  30.2046 |
|  69 | PP  | RADIAT |  3.2085 | 14.7164      | 1.01782     | 39.6636     | 2.22384     |  9.99421     | 0.603192    |  29.3553 |
|  70 | PP  | RADIAT |  3.105  | 16.2765      | 1.04777     | 39.5725     | 2.08858     | 10.7377      | 0.73841     |  28.4951 |
|  71 | PP  | RADIAT |  3.0015 | 17.9534      | 1.07825     | 39.4607     | 1.95279     | 11.5091      | 0.900494    |  27.625  |
|  72 | PP  | RADIAT |  2.898  | 19.7473      | 1.10926     | 39.3247     | 1.81712     | 12.3053      | 1.09378     |  26.7464 |
|  73 | PP  | RADIAT |  2.7945 | 21.6567      | 1.14075     | 39.1609     | 1.68224     | 13.1225      | 1.32303     |  25.8605 |
|  74 | PP  | RADIAT |  2.691  | 23.6783      | 1.17271     | 38.9655     | 1.54886     | 13.9566      | 1.59341     |  24.9689 |
|  75 | PP  | RADIAT |  2.5875 | 25.807       | 1.2051      | 38.7347     | 1.41768     | 14.8023      | 1.90818     |  24.0731 |
|  76 | PP  | RADIAT |  2.484  | 28.0352      | 1.2379      | 38.4666     | 1.28943     | 15.6543      | 2.24686     |  23.1746 |
|  77 | PP  | RADIAT |  2.3805 | 30.3534      | 1.27109     | 38.1609     | 1.16483     | 16.5062      | 2.63357     |  22.2747 |
|  78 | PP  | RADIAT |  2.277  | 32.7496      | 1.30464     | 37.8159     | 1.04458     | 17.3513      | 3.07244     |  21.3748 |
|  79 | PP  | RADIAT |  2.1735 | 35.2099      | 1.33852     | 37.4313     | 0.929332    | 18.1826      | 3.56736     |  20.4762 |
|  80 | PP  | RADIAT |  2.07   | 37.7179      | 1.37273     | 37.0075     | 0.819708    | 18.9923      | 4.122       |  19.5799 |
|  81 | PP  | RADIAT |  1.9665 | 40.2557      | 1.40726     | 36.5466     | 0.716266    | 19.7729      | 4.73967     |  18.6871 |
|  82 | PP  | CONVEC |  1.863  | 28.6199      | 1.22765     |  1.86464    | 0.498088    | 16.1143      | 2.2371      |  24.5604 |
|  83 | PP  | CONVEC |  1.7595 | 30.225       | 1.25474     |  1.69602    | 0.427661    | 16.6506      | 2.52245     |  23.5114 |
|  84 | PP  | CONVEC |  1.656  | 31.8164      | 1.28076     |  1.52218    | 0.363027    | 17.1712      | 2.82389     |  22.5658 |
|  85 | PP  | CONVEC |  1.5525 | 33.3838      | 1.30563     |  1.34589    | 0.304256    | 17.6738      | 3.13902     |  21.7141 |
|  86 | PP  | CONVEC |  1.449  | 34.9166      | 1.32929     |  1.17018    | 0.251362    | 18.1564      | 3.46486     |  20.9481 |
|  87 | PP  | CONVEC |  1.3455 | 36.4041      | 1.35166     |  0.998245   | 0.204297    | 18.6166      | 3.79793     |  20.2605 |
|  88 | PP  | CONVEC |  1.242  | 37.8356      | 1.37267     |  0.833311   | 0.162953    | 19.0524      | 4.13426     |  19.6449 |
|  89 | PP  | CONVEC |  1.1385 | 39.2005      | 1.39227     |  0.678475   | 0.127164    | 19.4619      | 4.46949     |  19.0958 |
|  90 | PP  | CONVEC |  1.035  | 40.4884      | 1.41039     |  0.536551   | 0.0967048   | 19.843       | 4.79891     |  18.6083 |
|  91 | PP  | CONVEC |  0.9315 | 41.6892      | 1.42697     |  0.409911   | 0.0712935   | 20.1941      | 5.11763     |  18.1782 |
|  92 | PP  | CONVEC |  0.828  | 42.7936      | 1.44198     |  0.300347   | 0.0505955   | 20.5134      | 5.42061     |  17.802  |
|  93 | PP  | CONVEC |  0.7245 | 43.7924      | 1.45535     |  0.208943   | 0.0342255   | 20.7993      | 5.70287     |  17.4764 |
|  94 | PP  | CONVEC |  0.621  | 44.6776      | 1.46704     |  0.135996   | 0.0217516   | 21.0506      | 5.95956     |  17.1988 |
|  95 | PP  | CONVEC |  0.5175 | 45.4418      | 1.47703     |  0.0809587  | 0.0127      | 21.2659      | 6.18612     |  16.967  |
|  96 | PP  | CONVEC |  0.414  | 46.0784      | 1.48527     |  0.0424441  | 0.00655958  | 21.4441      | 6.37838     |  16.7793 |
|  97 | PP  | CONVEC |  0.3105 | 46.5812      | 1.49173     |  0.0182651  | 0.00278746  | 21.5842      | 6.53251     |  16.6342 |
|  98 | --  | CENTER |  0.207  | 46.942       | 1.49634     |  0.0055268  | 0.000814872 | 21.6843      | 6.69316     |  16.5319 |
|  99 | --  | CENTER |  0.1035 | 47.1574      | 1.49909     |  0.00069085 | 0.000101859 | 21.744       | 6.72388     |  16.4714 |
| 100 | --  | CENTER |  0      | 47.2293      | 1.5         |  0          | 0           | 21.7639      | 6.73413     |  16.4513 |

License
----

**Free Software!** 
For the benefit of everyone.