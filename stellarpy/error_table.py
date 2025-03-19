from stellarpy.star import Star
import matplotlib.pyplot as plt
from scipy import optimize
import numpy as np

def error_table(star:Star, n:int, dR=0.5, dL=5.0, numbering=False):
    """
    Table containing the total relative error for total luminosity and total radius variations.
    Given a Star object, total relative error for variations of Ltot and Rtot is computed.

    ## Arguments:

        * star (Star): 
            Star object for which the minimum must be found.

        * n (int): 
            Size of the maximum variation. The output table length will be (2*n+1).

        * dR (float, default = 0.5): 
            Total radius variation.

        * dL (float, default = 5.0): 
            Total luminosity variation.
        
        * numbering (bool, default = False):
            Enables table numbering.
    """

    # Selecting values from the star
    Mtot, Rtot, Ltot, Tc, X, Y = star.parameters()
    
    # Defining radius and luminosity variations
    R = [Rtot + j*dR for j in range(-n, n+1)]
    L = [Ltot + i*dL for i in range(-n, n+1)]

    # Defining the error matrix
    error_matrix = np.array([[optimize.minimize(lambda x: Star(Mtot=Mtot, Rtot=r, Ltot=l, Tc=x[0], X=X, Y=Y).error(), x0=Tc).fun for r in R] for l in L])

    # Customizing the plot
    plt.rcParams["font.family"] = "serif"                   # Changing font
    plt.figure(figsize=(7, 7))                              # Changing figure size
    plt.imshow(error_matrix, cmap="coolwarm", vmax=np.max(error_matrix), origin="lower")     # Plotting error matrix

    if numbering:
        # Adding anotations for each element in the matrix
        for (i, j), value in np.ndenumerate(error_matrix):
            plt.text(j, i, f"{value/100:.2f}", ha="center", va="center", color="black", fontsize=10)

    plt.title("Summary table", fontsize=20, weight="bold")                                                         # Title
    plt.xlabel("R$_\\text{total}\\pm\\delta$R", fontsize=18)
    plt.ylabel("L$_\\text{total}\\pm\\delta$L", fontsize=18)
    plt.xticks(ticks=range(2*n+1), labels=R, fontsize=14) # Minor ticks (x axis)
    plt.yticks(ticks=range(2*n+1), labels=L, fontsize=14) # Minor ticks (y axis)

    # plt.xticks(ticks=range(2*n+1), labels=[f"{j}·$\\delta$R" for j in range(-n, n+1)], fontsize=16) # Minor ticks (x axis)
    # plt.yticks(ticks=range(2*n+1), labels=[f"{i}·$\\delta$L" for i in range(-n, n+1)], fontsize=16) # Minor ticks (y axis)

    # Printing the minimum
    min_index = np.unravel_index(np.argmin(error_matrix), error_matrix.shape)
    Lmin, Rmin = L[min_index[0]], R[min_index[1]]
    result = optimize.minimize(lambda x: Star(Mtot=Mtot, Rtot=Rmin, Ltot=Lmin, Tc=x[0], X=X, Y=Y).error(), x0=Tc)
    print(f"Minimum found at: \n   Rtot = {Rmin:.4f}\n   Ltot = {Lmin:.4f}\n   Tc =   {result.x[0]:.4f}")
    print(f"Error: {result.fun:.4f} %")

    plt.show()