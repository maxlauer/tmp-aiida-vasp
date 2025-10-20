import numpy as np

def Murnaghan(V, E0, V0, K0, K1):
    """
    Input:
        V - pd.Series - Volume - A^3
    ------------------------------------------
    Parameters
    ----------
    E0 - float - Energy at equilibrium Volume - eV
    V0 - float - equilibrium volume - A^3
    K0 - float - Bulk modulus - GPa
    K1 - float - pressure derivative of K0 - no dim

    ------------------------------------------
    Return - energy as function of volume - pd.Series - eV
    """
    e = 1.60217663 * 10 ** (-19)
    e_val = 1.60217663 #C divide Joule to get eV
    conversion_factor = e_val * 10**(-2)## G[Pa] (10**(9)) * [A^3](10**(-30)) * [1/e] 1/1.602 * 10**19

    return E0 + ( (K0 * V) / K1 * (((V0 / V) ** K1) / (K1 - 1) + 1) - (K0 * V0) / (K1 - 1) ) * conversion_factor
