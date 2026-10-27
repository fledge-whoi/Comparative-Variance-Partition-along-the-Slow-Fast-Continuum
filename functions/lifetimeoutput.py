import numpy as np


def lifetimeoutput(N, I_ngroupstage, selectionvect):
    """
    **********************************************************************
    Function that computes and returns the mean and variance of life expectancy.
    (longevity = number of time steps lived), computed from the fundamental matrix
    N = (I - U)^-1 (Caswell 2009, Oikos, doi:10.1111/j.1600-0706.2009.17620.x):
        eta1' = 1' N                     (mean longevity)
        eta2' = eta1' (2N - I)           (second moment)
        var   = eta2 - eta1 .* eta1
    INPUT
    N: The fundamental matrix.
    I_ngroupstage: An identity matrix of dimension ngroupstagexngroupstage
    selectionvect: Vector containing the stage from which we compute expected longevity
                   (boolean; here the juvenile state of each group).
    OUTPUT
    meanLE: mean longevity for each heterogeneity group.
    varLE: within-group variance in longevity for each heterogeneity group
           (individual stochasticity).
    **********************************************************************
    """
    eta1tilde = N.sum(axis=0)
    eta2tilde = eta1tilde @ (2 * N - I_ngroupstage)
    varetatilde = eta2tilde - eta1tilde * eta1tilde
    meanLE = eta1tilde[selectionvect]
    varLE = varetatilde[selectionvect]
    return meanLE, varLE
