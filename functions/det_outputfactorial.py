import numpy as np
from scipy import linalg


def _dominant_eig(Amat):
    """Dominant eigenvalue with its right (w) and left (v) eigenvectors.
    Equivalent of [wmat, dmat, vmat] = eig(Amat) followed by selecting the
    largest eigenvalue (MATLAB's max compares complex numbers by modulus)."""
    lam, vl, vr = linalg.eig(Amat, left=True, right=True)
    imax = int(np.argmax(np.abs(lam)))
    return lam[imax], vr[:, imax], vl[:, imax]


def generation_time(Umat, Fmat, Amat):
    """Generation time G = 1 + (v'Uw)/(v'Fw) (Bienvenu & Legendre 2015).
    Shared by det_outputfactorial and minusparamfinder_ConstantG."""
    lambda1, w, v = _dominant_eig(Amat)
    vc = np.conj(v)
    GT = 1 + (vc @ Umat @ w) / (vc @ Fmat @ w)
    return np.real(GT), lambda1, w


def det_outputfactorial(Umat, Fmat, Amat):
    """
    **********************************************************************
    Function that returns the population's growth rate, it's stable stage
    distribution and its generation time.
    INPUT
    Umat: transition matrix (survival and maturation).
    Fmat: fertility matrix.
    Amat: population matrix (Amat = Umat + Fmat).
    OUTPUT
    lambda1: population's growth rate (dominant eigenvalue of Amat).
    SSD: stable stage distribution (right eigenvector, scaled to sum to 1).
         States are ordered group by group: [J1 A1 J2 A2 ...].
    GT: generation time (in years), G = 1 + (v'Uw)/(v'Fw)
        (Bienvenu & Legendre 2015).
    **********************************************************************
    """
    # Calculate the rate of population change (lambda1 = the dominant eigenvalue)
    # v: left eigenvector (reproductive values), w: right eigenvector (stable structure)
    GT, lambda1, w = generation_time(Umat, Fmat, Amat)

    RawAge = w
    SSD = np.real(RawAge / RawAge.sum())

    return float(np.real(lambda1)), SSD, float(GT)
