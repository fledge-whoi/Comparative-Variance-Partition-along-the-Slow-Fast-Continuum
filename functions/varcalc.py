import numpy as np


def varcalc(X, PI):
    """
    **********************************************************************
    Function that returns the between group variance of a life history
    outcome, i.e. the variance of the group means weighted by the group
    proportions (law of total variance):
        Vb = PI' * (X .* X) - (PI' * X)^2
    INPUT
    X: withingroup mean of any life history outcome (one value per group).
    PI: relative distribution of juvenile stages (group proportions, sum to 1).
    OUTPUT
    VarX: Between group variance.
    **********************************************************************
    """
    X = np.asarray(X, dtype=float)
    PI = np.asarray(PI, dtype=float)
    return float(PI @ (X * X) - (PI @ X) ** 2)
