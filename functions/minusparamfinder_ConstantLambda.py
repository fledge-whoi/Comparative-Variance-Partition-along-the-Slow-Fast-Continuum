import numpy as np


def minusparamfinder_ConstantLambda(thetax, Y, one_ngroup_nparam, Py, Rcomb):
    """
    **********************************************************************
    Function that computes the difference between the population growth rate
    of any baseline population and the growth rate of a corresponding
    heterogeneous population (for which we know the phenotype distribution at
    birth and the vital rates of high quality individuals).
    In practice the condition lambda = 1 is written with the net reproductive
    rate: for this life cycle the R0 of a group is
        R0 = f * sj * g / ((1 - sa) * (1 - sj * (1 - g)))
    and, because newborns are allocated to groups independently of their
    parent's group, lambda = 1 if and only if sum(Rcomb .* R0) = 1.
    INPUT
    thetax: Vector containing the vital rate values for high-quality individuals.
    Y: Multiplicative factor that defines the vital rates of low quality individuals.
    one_ngroup_nparam: Matrix of dimension ngroupxnparam filled with ones.
    Py: Template that maps low vital rates values in each heterogeneity group.
    Rcomb: Phenotype distribution at birth
    OUTPUT
    diff: R0 of the heterogeneous population minus 1 (zero for a stable population).
    **********************************************************************
    """
    Ymat = Y * one_ngroup_nparam
    thetaxy = thetax * (Ymat ** Py)

    sjvect = thetaxy[:, 0]
    savect = thetaxy[:, 1]
    gvect = thetaxy[:, 2]
    fvect = thetaxy[:, 3]

    rzerovect = (fvect * sjvect * gvect) / ((1 - savect) * (1 - (sjvect * (1 - gvect))))

    diff = (Rcomb @ rzerovect) - 1  # - 1 means we are looking for a stable population
    return float(diff)
