from stellarpy import Star
from stellarpy import error_table
from stellarpy import find_minimum

# Initializing a Star object named after Mizar Aa
Altair = Star(Mtot=5.0, Rtot=11.5, Ltot=40.0, Tc=1.5, X=0.80, Y=0.16, solar_units=False)

# Printing total relative error of the calculation
print(Altair.error())

# Finding the minimum relative error using error_table
error_table(star=Altair, n=5, dR=0.5, dL=5.0, numbering=True)

# Finding the minimum relative error using find_minimum and redefining star atributes
Rmin, Lmin, Tmin, error_min = find_minimum(star=Altair, x0=None)
Altair.redefine(Rtot=Rmin, Ltot=Lmin, Tc=Tmin)

# Now that the star is well defined, let's plot some cool graphs!
# First, variables throughout the star
Altair.visualize(merge=True, figsize=(10, 6))

# Plotting both energy generation rate and opacity
Altair.visualize(which=["kappa", "epsilon"], normalize=False)
Altair.visualize(x_axis="T", which=["kappa", "epsilon"], normalize=False)

# Plotting the Teperature-Density Diagram
Altair.TDD()

# Plotting Altair in the Hertzsprung–Russell Diagram
Altair.HR()

# Redefining Altair using real parameters
Mtot=3.70
Rtot=12.425
Ltot=40.7
Tc = 1.5
X=0.710
Y=0.271

Altair.redefine(Mtot=Mtot, Rtot=Rtot, Ltot=Ltot, Tc=Tc, X=X, Y=Y)
find_minimum(star=Altair, x0=None)



