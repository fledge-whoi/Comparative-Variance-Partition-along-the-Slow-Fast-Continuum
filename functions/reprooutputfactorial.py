import numpy as np


def reprooutputfactorial(U, F, N, selectionvect, s1, t0, Z):
    """
    **********************************************************************
    Function that computes and returns the mean and variance of lifetime
    reproductive output (LRO) with a Markov chain with rewards. For details
    about the equations, see Van Daalen and Caswell 2017
    (doi: 10.1007/s12080-017-0335-2) and Van Daalen and Caswell 2020.
    INPUT
    U: The transition matrix (t x t, t = number of transient states).
    F: The fertility matrix (t x t).
    N: The fundamental matrix, N = (I - U)^-1.
    selectionvect: Boolean vector selecting the starting states for which
                   results are returned (here the juvenile state of each group).
    s1, t0, Z: utility vectors and matrices to perform calculations
               (s1: (t+1) ones; t0: t zeros; Z = [I 0], t x (t+1)).
    OUTPUT
    meanLRO: mean lifetime reproductive output for each heterogeneity group.
    varLRO: within-group variance in lifetime reproductive output for each
            heterogeneity group (individual stochasticity).
    **********************************************************************
    """
    # theta = sj1 sa1 g1 f1 sj2 sa2 g2 f2 r
    # We consider a situation in which the only absorbing state is death

    f = F.sum(axis=0)       # Horizontal fertility vect
    m = 1 - U.sum(axis=0)   # Horizontal mortality vect
    R1 = np.outer(s1, np.append(f, 0.0))  # Rewards correspond to state-based reproductive output
    R1tilde = Z @ R1 @ Z.T
    P = np.block([[U, t0.reshape(-1, 1)],
                  [m.reshape(1, -1), np.ones((1, 1))]])  # Markov chain with absorbing state

    # We suppose a poisson distribution of offspring
    R2f = R1 ** 2
    # R3f = R1 ** 3

    R2p = R1 + R2f
    # R2ptilde = Z @ R2p @ Z.T
    # R3p = R1 + 3 * R2f + R3f
    # R3ptilde = Z @ R3p @ Z.T

    # calculation of moments vectors rho for the lifetime reproductive output
    # Formulas given in Van Daalen and Caswell 2017 (Theoretical Ecology) and
    # Van Daalen and Caswell 2020 (Ecological modelling)
    rhotilde1 = N.T @ Z @ (P * R1).T @ s1
    rhotilde2 = N.T @ ((Z @ (P * R2p).T @ s1) + 2 * (U * R1tilde).T @ rhotilde1)

    rhotildemean = rhotilde1
    rhotildevar = rhotilde2 - rhotilde1 * rhotilde1

    # We did not gather juveniles and adults together in the model structure,
    # so equations 26 and 27 from Van Daalen and Caswell 2020 cannot be applied
    # directly, but the principle here is simple: we decide that the juvenile
    # population is the starting point of our cohort
    rhomean = rhotildemean[selectionvect]
    rhovar = rhotildevar[selectionvect]

    # partitioning of variance
    meanLRO = rhomean
    varLRO = rhovar
    return meanLRO, varLRO
