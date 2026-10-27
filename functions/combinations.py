import itertools
import numpy as np


def combinations(*vectors):
    """
    **********************************************************************
    Python equivalent of MATLAB's combinations(...).Variables: all
    combinations of the elements of the input vectors, one combination per
    row. As in MATLAB, the FIRST input varies slowest and the LAST input
    varies fastest, e.g. combinations([1, 0], [1, 0]) gives
        [[1, 1], [1, 0], [0, 1], [0, 0]].
    This row ordering defines the order of the heterogeneity groups and of
    the birth distributions throughout the project.
    INPUT
    *vectors: any number of 1-D sequences.
    OUTPUT
    2-D float array of shape (prod(len(v)), len(vectors)).
    **********************************************************************
    """
    return np.array(list(itertools.product(*vectors)), dtype=float)
