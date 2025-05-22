<p align="center">
  <img src="images/StellarPyLogo.png" alt="StarPyLogo" width="300px" height="auto">
</p>

---

Python module for solving the stellar-interior equations.

 - Author: **Nuno Cerviño Luridiana**
 - Last update: May, 2025

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

    * `dR` (QuantityLike): 
        Total radius variation.

    * `dL` (QuantityLike): 
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
```sh
[Output]: Effective temperature: Teff = 8288 K.
```
<p align="center">
  <img src="images/usage_example5.png" alt="StarPyLogo" width="450px" height="auto">
</p>


## Data table

Finally, to access stellar-interior data using model units we invoke the get method.

```sh
star.get()
```


|     | E   | fase   |         r |            P |           T |          l |           m |          rho |       epsilon |    kappa |
|----:|:----|:-------|----------:|-------------:|------------:|-----------:|------------:|-------------:|--------------:|---------:|
| -11 | --  | ^^^^^^ | 11.2406   |  1.48924e-13 | 0.000486281 | 42.5688    | 5           |  2.11687e-10 |   0           | 26.0884  |
| -10 | --  | ^^^^^^ | 11.1393   |  2.75213e-09 | 0.00490702  | 42.5688    | 5           |  3.87675e-07 |   0           | 14.6375  |
|  -9 | --  | ^^^^^^ | 11.0381   |  4.37752e-08 | 0.00940887  | 42.5688    | 5           |  3.21593e-06 |   0           | 12.439   |
|  -8 | --  | ^^^^^^ | 10.9368   |  2.36569e-07 | 0.0139941   | 42.5688    | 5           |  1.1685e-05  |   0           | 11.2638  |
|  -7 | --  | ^^^^^^ | 10.8355   |  8.04571e-07 | 0.018665    | 42.5688    | 5           |  2.97956e-05 |   0           | 10.4812  |
|  -6 | --  | ^^^^^^ | 10.7343   |  2.11232e-06 | 0.0234241   | 42.5688    | 5           |  6.23322e-05 |   0           |  9.90272 |
|  -5 | --  | ^^^^^^ | 10.633    |  4.69975e-06 | 0.0282738   | 42.5688    | 5           |  0.000114897 |   0           |  9.44765 |
|  -4 | --  | ^^^^^^ | 10.5317   |  9.32097e-06 | 0.0332168   | 42.5688    | 5           |  0.000193964 |   0           |  9.07467 |
|  -3 | --  | ^^^^^^ | 10.4305   |  1.69884e-05 | 0.0382557   | 42.5688    | 5           |  0.000306953 |   0           |  8.75983 |
|  -2 | --  | ^^^^^^ | 10.3292   |  2.90233e-05 | 0.0433935   | 42.5688    | 5           |  0.000462315 |   0           |  8.48817 |
|  -1 | --  | ^^^^^^ | 10.2279   |  4.71138e-05 | 0.0486329   | 42.5688    | 5           |  0.000669628 |   0           |  8.24968 |
|   0 | --  | START  | 10.1267   |  7.33818e-05 | 0.0539772   | 42.5688    | 5           |  0.00093971  |   0           |  8.03743 |
|   1 | --  | START  | 10.0254   |  0.000110459 | 0.0594295   | 42.5688    | 5           |  0.00128474  |   0           |  7.84638 |
|   2 | --  | START  |  9.92413  |  0.000161576 | 0.064993    | 42.5688    | 5           |  0.0017184   |   0           |  7.67279 |
|   3 | --  | RADIAT |  9.82286  |  0.000231124 | 0.0706796   | 42.5688    | 4.99986     |  0.00226031  |   0           |  7.52487 |
|   4 | --  | RADIAT |  9.7216   |  0.000323627 | 0.0764963   | 42.5688    | 4.99954     |  0.00292429  |   0           |  7.3814  |
|   5 | --  | RADIAT |  9.62033  |  0.000444828 | 0.0824399   | 42.5688    | 4.99914     |  0.00372967  |   0           |  7.2452  |
|   6 | --  | RADIAT |  9.51906  |  0.000601587 | 0.0885109   | 42.5688    | 4.99865     |  0.00469805  |   0           |  7.11687 |
|   7 | --  | RADIAT |  9.4178   |  0.000802041 | 0.0947112   | 42.5688    | 4.99804     |  0.00585344  |   0           |  6.99627 |
|   8 | --  | RADIAT |  9.31653  |  0.00105582  | 0.101045    | 42.5688    | 4.99731     |  0.00722257  |   0           |  6.88267 |
|   9 | --  | RADIAT |  9.21526  |  0.00137424  | 0.107515    | 42.5688    | 4.99643     |  0.00883512  |   0           |  6.77542 |
|  10 | --  | RADIAT |  9.114    |  0.00177065  | 0.114126    | 42.5688    | 4.99537     |  0.0107241   |   0           |  6.67384 |
|  11 | --  | RADIAT |  9.01273  |  0.0022606   | 0.120884    | 42.5688    | 4.99413     |  0.0129263   |   0           |  6.5773  |
|  12 | --  | RADIAT |  8.91146  |  0.00286233  | 0.127792    | 42.5688    | 4.99267     |  0.0154822   |   0           |  6.48531 |
|  13 | --  | RADIAT |  8.8102   |  0.00359702  | 0.134856    | 42.5688    | 4.99096     |  0.0184369   |   0           |  6.3973  |
|  14 | --  | RADIAT |  8.70893  |  0.00448934  | 0.142082    | 42.5688    | 4.98898     |  0.0218404   |   0           |  6.31299 |
|  15 | --  | RADIAT |  8.60766  |  0.00556798  | 0.149474    | 42.5688    | 4.9867      |  0.0257482   |   0           |  6.23198 |
|  16 | --  | RADIAT |  8.5064   |  0.00686614  | 0.157039    | 42.5688    | 4.98407     |  0.0302219   |   0           |  6.15403 |
|  17 | --  | RADIAT |  8.40513  |  0.00842227  | 0.164782    | 42.5688    | 4.98107     |  0.0353294   |   0           |  6.07878 |
|  18 | --  | RADIAT |  8.30386  |  0.0102808   | 0.172709    | 42.5688    | 4.97765     |  0.0411461   |   0           |  6.00604 |
|  19 | --  | RADIAT |  8.2026   |  0.012493    | 0.180827    | 42.5688    | 4.97377     |  0.0477552   |   0           |  5.93554 |
|  20 | --  | RADIAT |  8.10133  |  0.015118    | 0.189142    | 42.5688    | 4.96939     |  0.0552489   |   0           |  5.86712 |
|  21 | --  | RADIAT |  8.00006  |  0.018224    | 0.197662    | 42.5688    | 4.96445     |  0.063729    |   0           |  5.80056 |
|  22 | --  | RADIAT |  7.8988   |  0.0218892   | 0.206392    | 42.5688    | 4.9589      |  0.0733082   |   0           |  5.73571 |
|  23 | --  | RADIAT |  7.79753  |  0.0262038   | 0.215341    | 42.5688    | 4.95269     |  0.084111    |   0           |  5.67241 |
|  24 | --  | RADIAT |  7.69626  |  0.0312711   | 0.224516    | 42.5688    | 4.94575     |  0.0962746   |   0           |  5.61053 |
|  25 | --  | RADIAT |  7.595    |  0.0372101   | 0.233926    | 42.5688    | 4.93803     |  0.109951    |   0           |  5.54993 |
|  26 | --  | RADIAT |  7.49373  |  0.0441567   | 0.243577    | 42.5688    | 4.92945     |  0.125307    |   0           |  5.49048 |
|  27 | --  | RADIAT |  7.39246  |  0.0522671   | 0.253479    | 42.5688    | 4.91995     |  0.142529    |   0           |  5.43208 |
|  28 | --  | RADIAT |  7.2912   |  0.0617197   | 0.26364     | 42.5688    | 4.90944     |  0.161818    |   0           |  5.37461 |
|  29 | --  | RADIAT |  7.18993  |  0.072719    | 0.27407     | 42.5688    | 4.89784     |  0.183401    |   0           |  5.31798 |
|  30 | --  | RADIAT |  7.08867  |  0.0854985   | 0.284778    | 42.5688    | 4.88508     |  0.207524    |   0           |  5.26208 |
|  31 | --  | RADIAT |  6.9874   |  0.100325    | 0.295773    | 42.5688    | 4.87105     |  0.234459    |   0           |  5.20682 |
|  32 | --  | RADIAT |  6.88613  |  0.117504    | 0.307066    | 42.5688    | 4.85567     |  0.264506    |   0           |  5.15211 |
|  33 | --  | RADIAT |  6.78487  |  0.137382    | 0.318667    | 42.5688    | 4.83883     |  0.297994    |   0           |  5.09787 |
|  34 | --  | RADIAT |  6.6836   |  0.160355    | 0.330586    | 42.5688    | 4.82044     |  0.335284    |   0           |  5.04402 |
|  35 | --  | RADIAT |  6.58233  |  0.186875    | 0.342836    | 42.5688    | 4.80037     |  0.376774    |   0           |  4.99046 |
|  36 | --  | RADIAT |  6.48107  |  0.217454    | 0.355426    | 42.5688    | 4.77851     |  0.422897    |   0           |  4.93712 |
|  37 | --  | RADIAT |  6.3798   |  0.252677    | 0.368369    | 42.5688    | 4.75475     |  0.474132    |   0           |  4.88392 |
|  38 | --  | RADIAT |  6.27853  |  0.293206    | 0.381677    | 42.5688    | 4.72895     |  0.530998    |   0           |  4.83079 |
|  39 | --  | RADIAT |  6.17727  |  0.339792    | 0.395362    | 42.5688    | 4.701       |  0.594066    |   0           |  4.77765 |
|  40 | PP  | RADIAT |  6.076    |  0.393288    | 0.409436    | 42.5688    | 4.67075     |  0.663957    |   0.000289337 |  4.72442 |
|  41 | PP  | RADIAT |  5.97473  |  0.454656    | 0.423913    | 42.5688    | 4.63806     |  0.741349    |   0.000397953 |  4.67103 |
|  42 | PP  | RADIAT |  5.87347  |  0.524988    | 0.438806    | 42.5688    | 4.60281     |  0.826976    |   0.000546105 |  4.61742 |
|  43 | PP  | RADIAT |  5.7722   |  0.605513    | 0.454128    | 42.5687    | 4.56483     |  0.921639    |   0.000747788 |  4.56349 |
|  44 | PP  | RADIAT |  5.67093  |  0.697615    | 0.469893    | 42.5687    | 4.52398     |  1.0262      |   0.00102182  |  4.50919 |
|  45 | PP  | RADIAT |  5.56967  |  0.802855    | 0.486114    | 42.5686    | 4.48011     |  1.1416      |   0.00139345  |  4.45444 |
|  46 | PP  | RADIAT |  5.4684   |  0.922981    | 0.502805    | 42.5686    | 4.43308     |  1.26885     |   0.00189652  |  4.39917 |
|  47 | PP  | RADIAT |  5.36713  |  1.05996     | 0.519981    | 42.5685    | 4.38273     |  1.40902     |   0.00257627  |  4.34331 |
|  48 | PP  | RADIAT |  5.26587  |  1.21597     | 0.537655    | 42.5683    | 4.32891     |  1.56327     |   0.00349309  |  4.28679 |
|  49 | PP  | RADIAT |  5.1646   |  1.39347     | 0.555842    | 42.5681    | 4.27148     |  1.73285     |   0.0047274   |  4.22953 |
|  50 | PP  | RADIAT |  5.06333  |  1.59517     | 0.574555    | 42.5677    | 4.2103      |  1.91907     |   0.00638608  |  4.17147 |
|  51 | PP  | RADIAT |  4.96207  |  1.82408     | 0.593809    | 42.5672    | 4.14523     |  2.12331     |   0.00861088  |  4.11255 |
|  52 | PP  | RADIAT |  4.8608   |  2.08354     | 0.613616    | 42.5665    | 4.07616     |  2.34704     |   0.0119169   |  4.05268 |
|  53 | PP  | RADIAT |  4.75953  |  2.3772      | 0.63399     | 42.5655    | 4.00296     |  2.59179     |   0.0154944   |  3.99183 |
|  54 | PP  | RADIAT |  4.65827  |  2.7091      | 0.654944    | 42.5641    | 3.92554     |  2.85915     |   0.0201104   |  3.92991 |
|  55 | PP  | RADIAT |  4.557    |  3.0836      | 0.67649     | 42.5622    | 3.84381     |  3.15075     |   0.0260546   |  3.86687 |
|  56 | PP  | RADIAT |  4.45573  |  3.50549     | 0.698639    | 42.5596    | 3.75771     |  3.46826     |   0.0336932   |  3.80266 |
|  57 | PP  | RADIAT |  4.35447  |  3.97989     | 0.721402    | 42.5561    | 3.6672      |  3.81338     |   0.0434876   |  3.73722 |
|  58 | PP  | RADIAT |  4.2532   |  4.51233     | 0.744789    | 42.5514    | 3.57226     |  4.18778     |   0.0560171   |  3.67051 |
|  59 | PP  | RADIAT |  4.15193  |  5.1087      | 0.768809    | 42.545     | 3.47293     |  4.59312     |   0.0720064   |  3.60249 |
|  60 | PP  | RADIAT |  4.05067  |  5.7752      | 0.793468    | 42.5365    | 3.36923     |  5.03099     |   0.0923579   |  3.53313 |
|  61 | PP  | RADIAT |  3.9494   |  6.51832     | 0.818773    | 42.5251    | 3.26128     |  5.50285     |   0.11819     |  3.46239 |
|  62 | PP  | RADIAT |  3.84813  |  7.34479     | 0.844727    | 42.5099    | 3.14919     |  6.01007     |   0.150882    |  3.39025 |
|  63 | PP  | RADIAT |  3.74687  |  8.2616      | 0.871333    | 42.4899    | 3.03315     |  6.55384     |   0.19213     |  3.31674 |
|  64 | PP  | RADIAT |  3.6456   |  9.27572     | 0.898591    | 42.4637    | 2.91336     |  7.13513     |   0.244001    |  3.24187 |
|  65 | PP  | RADIAT |  3.54433  | 10.3941      | 0.926498    | 42.4295    | 2.79011     |  7.75459     |   0.309004    |  3.16565 |
|  66 | PP  | RADIAT |  3.44307  | 11.6235      | 0.955051    | 42.3856    | 2.66371     |  8.41253     |   0.381269    |  3.08811 |
|  67 | PP  | RADIAT |  3.3418   | 12.9703      | 0.984243    | 42.3302    | 2.53455     |  9.10888     |   0.472729    |  3.00931 |
|  68 | PP  | RADIAT |  3.24053  | 14.4404      | 1.01406     | 42.2604    | 2.40306     |  9.84308     |   0.58427     |  2.92929 |
|  69 | PP  | RADIAT |  3.13927  | 16.039       | 1.0445      | 42.173     | 2.26972     | 10.6141      |   0.719721    |  2.84815 |
|  70 | PP  | RADIAT |  3.038    | 17.77        | 1.07554     | 42.0646    | 2.13507     | 11.4203      |   0.883462    |  2.76596 |
|  71 | PP  | RADIAT |  2.93673  | 19.6365      | 1.10716     | 41.931     | 1.9997      | 12.2594      |   1.08046     |  2.68285 |
|  72 | PP  | RADIAT |  2.83547  | 21.64        | 1.13934     | 41.7679    | 1.86423     | 13.1286      |   1.3163      |  2.59893 |
|  73 | PP  | RADIAT |  2.7342   | 23.7801      | 1.17205     | 41.5705    | 1.72933     | 14.0243      |   1.59714     |  2.51435 |
|  74 | PP  | RADIAT |  2.63293  | 26.0548      | 1.20528     | 41.3338    | 1.59569     | 14.9423      |   1.92735     |  2.42924 |
|  75 | PP  | RADIAT |  2.53167  | 28.4597      | 1.23898     | 41.0547    | 1.46402     | 15.8775      |   2.28684     |  2.34377 |
|  76 | PP  | RADIAT |  2.4304   | 30.9882      | 1.27314     | 40.7313    | 1.33505     | 16.8243      |   2.70167     |  2.25805 |
|  77 | PP  | RADIAT |  2.32913  | 33.6312      | 1.30772     | 40.3602    | 1.2095      | 17.7765      |   3.17756     |  2.17224 |
|  78 | PP  | RADIAT |  2.22787  | 36.3773      | 1.3427      | 39.9392    | 1.08806     | 18.7271      |   3.72025     |  2.08645 |
|  79 | PP  | RADIAT |  2.1266   | 39.2125      | 1.37805     | 39.4669    | 0.97143     | 19.6687      |   4.33536     |  2.00082 |
|  80 | PP  | RADIAT |  2.02533  | 42.1203      | 1.41376     | 38.9435    | 0.86023     | 20.5935      |   5.02832     |  1.91546 |
|  81 | PP  | RADIAT |  1.92407  | 45.0821      | 1.44982     | 38.3707    | 0.755046    | 21.4934      |   5.80428     |  1.83046 |
|  82 | PP  | RADIAT |  1.8228   | 48.0774      | 1.48624     | 37.7521    | 0.656389    | 22.3599      |   6.66804     |  1.74591 |
|  83 | PP  | CONVEC |  1.72153  | 51.0822      | 1.52299     | 37.1123    | 0.564714    | 23.184       |   7.62354     |  1.66191 |
|  84 | CN  | CONVEC |  1.62027  | 54.0806      | 1.55814     | 36.3193    | 0.480282    | 23.9912      |  11.4869      |  1.58778 |
|  85 | CN  | CONVEC |  1.519    | 57.0492      | 1.5918      | 35.1817    | 0.403254    | 24.7729      |  18.1882      |  1.52133 |
|  86 | CN  | CONVEC |  1.41773  | 59.9664      | 1.62387     | 33.6338    | 0.333713    | 25.5253      |  26.6038      |  1.46183 |
|  87 | CN  | CONVEC |  1.31647  | 62.8101      | 1.65425     | 31.6319    | 0.271657    | 26.2449      |  38.1853      |  1.40864 |
|  88 | CN  | CONVEC |  1.2152   | 65.5583      | 1.68283     | 29.1381    | 0.217       | 26.928       |  53.3285      |  1.3612  |
|  89 | CN  | CONVEC |  1.11393  | 68.1887      | 1.70952     | 26.1599    | 0.169571    | 27.5711      |  72.4802      |  1.31903 |
|  90 | CN  | CONVEC |  1.01267  | 70.6794      | 1.73423     | 22.7573    | 0.129114    | 28.171       |  95.8848      |  1.28171 |
|  91 | CN  | CONVEC |  0.9114   | 73.0093      | 1.75687     | 19.0492    | 0.0952942   | 28.7246      | 123.485       |  1.24888 |
|  92 | CN  | CONVEC |  0.810133 | 75.1581      | 1.77738     | 15.2094    | 0.0676969   | 29.2289      | 154.838       |  1.22024 |
|  93 | CN  | CONVEC |  0.708867 | 77.1067      | 1.79567     | 11.4516    | 0.045835    | 29.6812      | 189.053       |  1.1955  |
|  94 | CN  | CONVEC |  0.6076   | 78.8375      | 1.81168     |  8.00326   | 0.0291527   | 30.0792      | 224.792       |  1.17446 |
|  95 | CN  | CONVEC |  0.506333 | 80.3346      | 1.82537     |  5.0709    | 0.0170326   | 30.4206      | 260.318       |  1.15692 |
|  96 | CN  | CONVEC |  0.405067 | 81.5837      | 1.83667     |  2.8037    | 0.00880235  | 30.7036      | 293.612       |  1.14272 |
|  97 | CN  | CONVEC |  0.3038   | 82.5717      | 1.84553     |  1.26231   | 0.0037423   | 30.9261      | 322.514       |  1.13177 |
|  98 | --  | CENTER |  0.202533 | 83.2812      | 1.85186     |  0.400695  | 0.00109459  | 31.0853      | 361.013       |  1.12405 |
|  99 | --  | CENTER |  0.101267 | 83.7052      | 1.85563     |  0.0500869 | 0.000136823 | 31.1802      | 362.851       |  1.11949 |
| 100 | --  | CENTER |  0        | 83.8469      | 1.85688     |  0         | 0           | 31.2118      | 363.465       |  1.11798 |

License
----

This repository is licensed under the MIT License. See the LICENSE file for further details.