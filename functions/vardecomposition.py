import itertools
import numpy as np

from .varcalc import varcalc


def _unique_rows_first(A):
    """Equivalent of MATLAB [C, ia] = unique(A, "rows", "first"):
    unique rows sorted in ascending lexicographic order, and the index of the
    first occurrence of each of them in A."""
    rows = [tuple(r) for r in np.asarray(A)]
    uniq = sorted(set(rows))
    first = [rows.index(u) for u in uniq]
    return np.array(uniq, dtype=float), np.array(first, dtype=int)


def vardecomposition(M, PI, Px, Py, ngroup, excluded):
    """
    **********************************************************************
    This is the code for factorial variance partition. The code is abstract in
    order to apply to all heterogeneity structures in the
    populations (i.e. all possible structures of the Px matrices). Method: Caswell (2025).
    Principle: for each subset S of vital rates, the variance of the group
    means marginalised over the vital rates NOT in S is computed; the
    contribution of the interaction among the factors in S is then obtained
    by subtracting the contributions of all its subsets of lower order.
    INPUT
    M: Matrix with mean values of any life history outcome. (2x2x2x2 array)
    PI: Matrix (2x2x2x2 array) with relative distribution of juvenile stages.
    Px, Py: Template matrices that allow us to differenciate between
            heterogeneity groups (ngroup x nparam).
    ngroup: number of groups
    excluded: Vector stores the vital rates for which there is or is not
              heterogeneity (1 = heterogeneity, 0 = none).
    OUTPUT
    sortedvarfactoriallist: vector of length ngroup with the variance
        components, in the row order of Px. In the quality scenario, element
        i is the component of the set of vital rates marked 1 in Px[i, :]:
        Python indices 7 = sj, 11 = sa, 13 = g, 14 = f (main effects;
        MATLAB rows 8, 12, 14, 15); the others are interactions; index 15
        (no factor) is 0. The components sum to the between-group variance.
    **********************************************************************
    M and PI should be in compact form (not a vectorized form), built with
    reshape(..., order="F") to reproduce MATLAB's column-major reshape.
    """
    Px = np.asarray(Px, dtype=float)
    Py = np.asarray(Py, dtype=float)

    # Here comb (for combination) refers to Py (not Px) because we sum on
    # dimensions that are NOT part of the combination of heterogeneity groups
    # we use.

    # Here we permute because reshape flipped the dimensions for the
    # parameters (e.g. sj corresponds to the 4th dimension instead of the 1st
    # and conversely f corresponds to the 1st dimension instead of the 4th).
    comb = np.flip(Py, axis=1)
    varfactoriallist = np.zeros(ngroup)
    MPI = M * PI

    for i in range(ngroup):
        contractionvect = tuple(np.nonzero(comb[i, :])[0])
        if contractionvect:
            PIpartial = PI.sum(axis=contractionvect, keepdims=True)
            with np.errstate(invalid="ignore", divide="ignore"):
                Mpartial = MPI.sum(axis=contractionvect, keepdims=True) / PIpartial
        else:
            Mpartial = M
            PIpartial = PI
        varfactoriallist[i] = varcalc(Mpartial.ravel(order="F"), PIpartial.ravel(order="F"))

    varfactoriallist[np.isnan(varfactoriallist)] = 0

    # We remove groups that are identical so we can automatize variance
    # partition calculations
    sortedPx, hetindex = _unique_rows_first(Px)
    sortedPx = np.flip(sortedPx, axis=0)
    hetindex = np.flip(hetindex, axis=0)
    sortedvarfactoriallist = varfactoriallist[hetindex].copy()

    newngroup = len(sortedvarfactoriallist)
    degreesvect = sortedPx.sum(axis=1)

    # We calculate variance for interactions (more than one factor)s
    for deg in range(2, int(degreesvect.max()) + 1):
        degreeindex = np.nonzero(degreesvect == deg)[0]
        # We get all subcombinations of factors (lower degrees)
        indexvarsubdegree = degreesvect < deg

        for ind in degreeindex:
            # We all factor combinations with a common factor with heterogeneity
            combi = sortedPx[ind, :]  # Now we need sortedPx
            indexvarcommonfact = (np.sum(np.abs(sortedPx - combi) - combi, axis=1) + degreesvect) == 0

            # We combine the two pieces of information
            indexvartotake = indexvarcommonfact & indexvarsubdegree
            vartotake = sortedvarfactoriallist[indexvartotake]
            sumvartotake = vartotake.sum()

            # Be aware that when variances are low, are very close to one
            # another the substraction sometimes yield results
            # that are very very close to zero, but negative.
            sortedvarfactoriallist[ind] = sortedvarfactoriallist[ind] - sumvartotake

    # We reinsert sortedvarfactoriallist into a 16x1 vector to ease code
    # implementation down the road.

    # This code does not work yet for more than 1 excluded parameter
    excludedpos = np.nonzero(np.asarray(excluded) == 0)[0] + 1  # MATLAB 1-based positions
    initialpowers = (2 ** (4 - excludedpos)).astype(int)
    powers = []
    supplement = np.full(ngroup, np.nan)
    queue = list(range(1, ngroup + 1))  # MATLAB 1-based group numbers

    for i in range(1, len(initialpowers) + 1):
        for dummy in itertools.combinations(initialpowers, i):
            powers.append(int(sum(dummy)))

    for i in range(newngroup):
        k = queue[0]
        supplement[k - 1] = sortedvarfactoriallist[i]
        queue = [q for q in queue if q != k]
        for p in powers:
            # p
            kbis = k + p
            queue = [q for q in queue if q != kbis]
            supplement[kbis - 1] = sortedvarfactoriallist[i]

    # supplement
    supplement[np.isnan(supplement)] = 0
    sortedvarfactoriallist = supplement
    return sortedvarfactoriallist
