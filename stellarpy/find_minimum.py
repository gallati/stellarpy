from stellarpy import Star
from scipy import optimize

def find_minimum(star, x0=None):
    """
    Total relative error minimum finder for Star objects.

    ## Arguments:
        * star (Star): 
            Star object for which the minimum must be found.

        * x0 (list, default = None): 
            List containing specific initial parameters required for the minimum search,
            listed as [Rtot, Ltot, Tc] in model units. If no list is provided, current 
            parameters of the Star object will be used.
    
    ## Returns:
        Optimized parameters and total relative error as a list following the order [Rtot, Ltot, Tc, error]

    """

    # Selecting values from the star
    Mtot, Rtot, Ltot, Tc, X, Y = star.parameters()

    if x0 is not None:
        Rtot, Ltot, Tc = x0

    # Finding the minimum
    def f(x):
        return Star(Mtot=Mtot, Rtot=x[0], Ltot=x[1], Tc=x[2], X=X, Y=Y).error()

    # Printing result
    result = optimize.minimize(f, x0=[Rtot, Ltot, Tc])
    print(f"Minimum found at: \n   Rtot = {result.x[0]:.4f}\n   Ltot = {result.x[1]:.4f}\n   Tc   = {result.x[2]:.4f}")
    print(f"Error: {result.fun:.4f} %")

    # Returning result
    return result.x[0], result.x[1], result.x[2], result.fun