from .popHeteromatfact import popHeteromatfact
from .det_outputfactorial import generation_time


def minusparamfinder_ConstantG(thetax, Y, one_ngroup_nparam, Py, Rcomb, G, ngroup, nstage):
    """
    **********************************************************************
    Function that computes the difference between the generation time of any
    baseline population and the generation time of a corresponding
    heterogeneous population (for which we know the phenotype distribution at
    birth and the vital rates of high quality individuals).
    Generation time follows Bienvenu & Legendre (2015): G = 1 + (v'Uw)/(v'Fw).
    INPUT
    thetax: Vector containing the vital rate values for high-quality individuals.
    Y: Multiplicative factor that defines the vital rates of low quality individuals.
    one_ngroup_nparam: Matrix of dimension ngroupxnparam filled with ones.
    Py: Template that maps low vital rates values in each heterogeneity group.
    Rcomb: Phenotype distribution at birth
    G: Generation time of the baseline population (target value, in years)
    ngroup: number of heterogeneity groups
    nstage: number of stages in the matrix model
    OUTPUT
    diff: difference between generation time in the baseline population and in
          the heterogeneous population.
    **********************************************************************
    """
    Ymat = Y * one_ngroup_nparam
    thetaxy = thetax * (Ymat ** Py)

    Uhet, Fhet, Ahet = popHeteromatfact(thetaxy, Rcomb, ngroup, nstage)
    g, _, _ = generation_time(Uhet, Fhet, Ahet)
    diff = G - g
    return float(diff)
