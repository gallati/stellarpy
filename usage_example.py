from stellarpy.star import Star
from stellarpy.error_table import error_table
from stellarpy.find_minimum import find_minimum

# Initializing a Star object named after Mizar Aa
Altair = Star(Mtot=5.0, Rtot=11.5, Ltot=40.0, Tc=1.5, X=0.80, Y=0.16, solar_units=False)

# Printing total relative error of the calculation
print(Altair.error())

# Finding the minimum relative error using error_table
error_table(star=Altair, n=5, dR=0.5, dL=5.0, numbering=True)

# Finding the minimum relative error using find_minimum
Rmin, Lmin, Tmin, error_min = find_minimum(star=Altair, x0=None)

# Redefining star atributes
Altair.redefine(Rtot=Rmin, Ltot=Lmin, Tc=Tmin)

# Now that the star is well defined, let's plot some cool graphs!
# First, variables throughout the star
Altair.visualize(merge=True)

# Plotting both energy generation rate and opacity
Altair.visualize(which=["kappa", "epsilon"], normalize=False)

# Ploting the Teperature-Density Diagram
Altair.TDD()



