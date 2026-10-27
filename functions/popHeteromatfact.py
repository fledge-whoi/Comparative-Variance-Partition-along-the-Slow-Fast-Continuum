import numpy as np

# The vec-permutation matrix only depends on (s, g): it is computed once and
# cached because this function is called many times by fzero/brentq.
_vecperm_cache = {}


def _vecperm(s, g):
    # vec permutation matrix
    key = (s, g)
    if key not in _vecperm_cache:
        p = np.zeros((s * g, s * g))
        for i in range(s):
            for j in range(g):
                e = np.zeros((s, g))
                e[i, j] = 1
                p = p + np.kron(e, e.T)
        _vecperm_cache[key] = p
    return _vecperm_cache[key]


def popHeteromatfact(thetamat, Rcomb, g, s):
    """
    **********************************************************************
    Function that computes the transition matrix (Utilde), fertility matrix
    (Ftilde), and returns Utilde, Ftilde and Atilde (Atilde = Ftilde + Utilde).
    To construct the matrix population models we follow formulas presented in
    Caswell et al. 2018 (doi: 10.1002/ecm.1306). The code we implemented
    reflect that omega (the number of age classes) = 1.
    States are ordered group by group: [J1 A1 J2 A2 ... Jg Ag].
    Heterogeneity is fixed: individuals never move between groups (U is block
    diagonal), and newborns of every group are distributed among the juvenile
    stages of all groups according to Rcomb (offspring group does not depend
    on the parent's group).
    INPUT
    thetamat: Matrix containing vital rates for each heterogeneity group.
               (g x 4, one row per group, columns sj sa g f)
    Rcomb: Phenotype distribution at birth (g x 1 vector, sums to 1).
    g: number of groups.
    s: number of stages (2: juveniles, adults).
    OUTPUT
    Utilde: (s*g) x (s*g) transition matrix (survival and maturation).
    Ftilde: (s*g) x (s*g) fertility matrix.
    Atilde: population projection matrix, Atilde = Utilde + Ftilde.
    **********************************************************************
    """
    thetamat = np.atleast_2d(np.asarray(thetamat, dtype=float))
    Rcomb = np.asarray(Rcomb, dtype=float).reshape(-1, 1)
    Is = np.eye(s)
    Ig = np.eye(g)

    # vec permutation matrix
    K = _vecperm(s, g)

    H = np.repeat(Rcomb, g, axis=1)
    bbH = np.kron(np.outer(Is[:, 0], Is[0, :]), H)
    for i in range(1, s):
        bbH = bbH + np.kron(np.outer(Is[:, i], Is[i, :]), H)

    # 1: sj, 2: sa, 3: g, 4: F  (columns 0, 1, 2, 3 in Python)
    # Block-diagonal matrix: block i is the 2x2 life cycle of group i.
    # (Same result as the MATLAB loop bbU = bbU + kron(Ig(:,i)*Ig(i,:), block),
    #  written with direct indexing because it is much faster.)
    sj, sa, gamma, f = thetamat[:, 0], thetamat[:, 1], thetamat[:, 2], thetamat[:, 3]
    J = np.arange(g) * s       # index of the juvenile state of each group
    A = J + 1                  # index of the adult state of each group
    bbU = np.zeros((s * g, s * g))
    bbU[J, J] = sj * (1 - gamma)
    bbU[A, J] = sj * gamma
    bbU[A, A] = sa

    # create block fertilty matrix
    bbF = np.zeros((s * g, s * g))
    bbF[J, A] = f

    Utilde = bbU  # Fixed heterogeneity, so no transitions between groups and no age classes
    Ftilde = K.T @ bbH @ K @ bbF
    Atilde = Utilde + Ftilde
    return Utilde, Ftilde, Atilde
