import numpy as np
from model import Model
import matplotlib.pyplot as plt
from scipy import optimize

def error_table(n, Mtot=5.0, Rtot=11.5, Ltot=70.0, Tc=2.0, dR=0.5, dL=5.0, tol=0.01):
    """
    Table containing the total relative error for total luminosity and total radius variations.
    Given the initial values of a model, it computes the smaller total relative error for variations
    of the total luminosity and total radius.

    ## Parameters:
        * n (float): Size of the maximum variation. The output table length will be (2*n+1).
        * Mtot (float, default = 5.0) : total mass of the star.
        * Rtot (float, default = 11.5) : total radius of the star.
        * Ltot (float, default = 70.0) : total luminosity of the star.
        * Tc (float, default = 2.0) : central temperature of the star.
        * dR (float, default = 0.5) : Total radius variation.
        * dL (float, default = 5.0) : Total luminosity variation.
    """

    # Defining radius and luminosity variations
    R = [Rtot + j*dR for j in range(-n, n+1)]
    L = [Ltot + i*dL for i in range(-n, n+1)]

    # Defining the error matrix
    error_matrix = np.array([[optimize.minimize(lambda x: Model(Mtot=Mtot, Rtot=r, Ltot=l, Tc=x[0]).error(), x0=Tc, tol=tol).x[0] for r in R] for l in L])
    # error_matrix = np.array([[Model(Mtot=Mtot, Rtot=r, Ltot=l, Tc=Tc).error() for r in R] for l in L])

    # Customizing the plot
    plt.rcParams["font.family"] = "serif"                   # Changing font
    plt.figure(figsize=(10, 7))                             # Changing figure size
    plt.imshow(error_matrix, cmap="coolwarm", vmax=100.0)   # Plotting error matrix

    # Adding anotations for each element in the matrix
    for (i, j), value in np.ndenumerate(error_matrix):
        plt.text(j, i, f"{value:.2f} %", ha="center", va="center", color="black", fontsize=12)

    plt.title("Summary table", fontsize=20)                                                         # Title
    plt.xticks(ticks=range(2*n+1), labels=[f"{j}·$\\delta$R" for j in range(-n, n+1)], fontsize=16) # Minor ticks (x axis)
    plt.yticks(ticks=range(2*n+1), labels=[f"{i}·$\\delta$L" for i in range(-n, n+1)], fontsize=16) # Minor ticks (y axis)

    plt.show()