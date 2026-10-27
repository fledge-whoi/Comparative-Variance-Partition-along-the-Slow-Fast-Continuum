"""
main.py -- Simulations for "Differential Contributions of Heterogeneity and
Stochasticity to Individuals' Life Histories along the Slow-fast Continuum"

Python translation of main.m (same structure, same variable names).

What this script does
---------------------
For each of the six synthetic species (rows of thetamatpop) it builds a set of
heterogeneous populations (two-stage matrix population models, 16 fixed
heterogeneity groups), computes the mean and variance of longevity (L) and
lifetime reproductive output (LRO) of newborns, partitions the total variance
into within-group (individual stochasticity) and between-group (fixed
individual heterogeneity) components, and decomposes the between-group
variance factorially (Caswell 2025). Results are saved in ./data and are read
by plotresults_articlescript.py.

How to reproduce the article
----------------------------
    python main.py                       # the three scenarios (default)
    python main.py --scenario quality    # a single scenario
    python main.py --quick               # small grid, ~1 min, to test the installation
Options: --workers N (number of CPU cores, default: all).
Then run:  python plotresults_articlescript.py

Notation (see the Methods of the article)
    sj : juvenile survival          sa : adult survival
    g  : maturation rate (gamma)    f  : fertility (recruits per adult)
    X  : multiplicative coefficient for "+" vital rates (Q+ or T+ in the article)
    Y  : multiplicative coefficient for "-" vital rates (Q- or T- in the article)

Indexing: Python arrays start at 0. "Species 1" of the article is index 0,
"group 1" (MATLAB row 1) is index 0, etc.
"""
import argparse
import csv
import gzip
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import brentq

from functions import (combinations, popHeteromatfact, det_outputfactorial,
                       lifetimeoutput, reprooutputfactorial, varcalc,
                       vardecomposition, minusparamfinder_ConstantG,
                       minusparamfinder_ConstantLambda)
from functions.det_outputfactorial import generation_time

HERE = os.path.dirname(os.path.abspath(__file__))

SCENARIOS = {
    # scenario name: output file name (without extension)
    "quality": "simulation_results_Quality_Gconstant",
    "tradeoff_sasj": "simulation_results_SaSjTradeoff_LAMBDAconstant",
    "tradeoff_saf": "simulation_results_SaFTradeoff_LAMBDAconstant",
}


def setup(scenario, quick=False):
    """Sections "Choosing vital rates ..." to "Basic utilities" of main.m.
    Returns every quantity needed by the core loops in a dictionary."""

    # %% Choosing vital rates, and determining where heterogeneity is introduced
    # Default code is for QUALITY heterogeneity. Use the "scenario" argument
    # ("tradeoff_sasj" or "tradeoff_saf") for TRADEOFF heterogeneity.

    # Number of demo parameters considered: sj sa g f
    # Population demographic rates
    # One row per species, ordered from fast to slow (Table 1a of the article).
    thetamatpop = np.array([  # sj    sa    g      f
        [0.01, 0.05, 0.9, 105.4500],
        [0.2, 0.1, 0.6, 6.9000],
        [0.5, 0.5, 0.4, 1.7500],
        [0.7, 0.75, 0.1, 1.3214],
        [0.8, 0.9, 0.075, 0.4333],
        [0.85, 0.95, 0.05, 0.2265],
    ])

    npop, nparam = thetamatpop.shape

    nstage = 2  # Two stages: Juveniles and Adults; within each class
    nclass = 2  # Number of classes per phenotypecombinations
    heterogeneitymat = np.array([1, 1, 1, 1])  # We select demographic parameters on which to introduce individual heterogeneity
    ngroup = nclass ** nparam  # Total number of heterogeneity groups
    heterogeneityfilter = np.tile(heterogeneitymat, (ngroup, 1))

    # Combinations of frequencies for each phenotype at birth that we want to
    # try. We could choose different frequency vects for different parameters
    Rvect = np.array([0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99])
    if quick:
        Rvect = np.array([0.1, 0.5, 0.9])  # reduced grid for a quick test
    Rallcomb = combinations(Rvect, Rvect, Rvect, Rvect)  # columns: Rsj Rsa Rg Rf
    ncomb = Rallcomb.shape[0]

    # Px: Combinations of parameter values in the groups (1 means greater parameter
    # value than base model, 0 means lower parameter value than base model)
    # Py(i,k) = 1 if vital rate k of group i is multiplied by Y ("-" value).
    # If both are 0, group i keeps the baseline value of vital rate k.
    # Rows follow the ordering of combinations(): sj varies slowest, f fastest.
    if scenario == "quality":
        # QUALITY
        Pxunfiltered = combinations([1, 0], [1, 0], [1, 0], [1, 0])  # Sj Sa g f
        Pyunfiltered = 1 - Pxunfiltered  # Py: Complementary matrix of Px
    # Or the matrices to be used when introducing trade-off heterogeneity:
    elif scenario == "tradeoff_saf":
        # TRADEOFF sa F
        # Groups 1-8: trade-off group, groups 9-16: baseline group.
        # NOTE: X is applied to f (column 4) and Y to sa (column 2), i.e. the
        # trade-off group has HIGHER fertility and LOWER adult survival.
        Pxunfiltered = np.vstack([np.hstack([np.zeros((8, 3)), np.ones((8, 1))]),
                                  np.zeros((8, 4))])  # Tradeoff sa F
        Pyunfiltered = np.vstack([np.hstack([np.zeros((8, 1)), np.ones((8, 1)), np.zeros((8, 2))]),
                                  np.zeros((8, 4))])  # Tradeoff sa F
    elif scenario == "tradeoff_sasj":
        # TRADEOFF sj sa
        # Groups 1-8: trade-off group (sa * X, sj * Y), groups 9-16: baseline.
        # The probability of being born in the trade-off group is Rsj.
        Pxunfiltered = np.vstack([np.hstack([np.zeros((8, 1)), np.ones((8, 1)), np.zeros((8, 2))]),
                                  np.zeros((8, 4))])  # Tradeoff sj sa
        Pyunfiltered = np.vstack([np.hstack([np.ones((8, 1)), np.zeros((8, 3))]),
                                  np.zeros((8, 4))])  # Tradeoff sj sa
    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    # #########################################################################
    # Note that there are other possible configurations for the combinations of parameter values. For instance:
    # Only best and worst possible individuals
    # Pxunfiltered = np.vstack([np.tile([1, 1, 1, 1], (8, 1)), np.tile([0, 0, 0, 0], (8, 1))])
    # #########################################################################

    # We modify Px and Py to only account for parameters on which we want to
    # introduce heterogeneity
    # The final Px matrix should remind you of figure 1 (from the method
    # section of the article)
    Px = Pxunfiltered * heterogeneityfilter
    Py = Pyunfiltered * heterogeneityfilter

    # For all combinations of sj,sa,g,f values, we create a Rmat that shows the
    # frequencies of for each group at birth (we consider independance between
    # phenotypes)
    Rmat = np.zeros((ngroup, nparam, ncomb))
    for comb in range(ncomb):
        Rj, Ra, Rg, Rf = Rallcomb[comb]
        # Rmatcombunfiltered = combinations([Rj, 1-Rj], [Ra, 1-Ra], [Rg, 1-Rg], [Rf, 1-Rf])
        # Rmatcombfiltered = Rmatcombunfiltered*heterogeneityfilter + (1 - heterogeneityfilter)  # We remove groups we do not need
        # Rmat[:, :, comb] = Rmatcombfiltered
        Rmat[:, :, comb] = combinations([Rj, 1 - Rj], [Ra, 1 - Ra], [Rg, 1 - Rg], [Rf, 1 - Rf])

    # %% Amplitude of heterogeneity
    # Here we choose the demographic parameter value in the (+) group first,
    # then, latter in the process, we calculate the value of parameters in the (-) group
    # according to specific criteria (either constante growth rate or
    # generation time)
    # demoratetocheck is a variable that simply allows us to determine the
    # vital rates to account for when calculating the maximum increase in vital
    # rate. Here we constructed it in a way that includes all parameters from
    # heterogeneitymat and excludes parameters you don't want (set their value
    # to 0, for instance here we excluded fertility because it is not a transition rate)
    transitionrate = np.array([1, 1, 1, 0])
    demoratetocheck = (transitionrate * heterogeneitymat).astype(bool)
    xmaxmat = 1.0 / thetamatpop[:, demoratetocheck] - 1  # We try to get the maximum increase in vital rates
    # for each species (since vital rates cannot exceed one).

    # While xmat is the maximum threshold for the vital rates of high-quality individuals
    # nx is the number of vital rates we want to test for between baseline parameters and
    # the threshold. Change this parameter according to the resolution you need for exploring the parameter space.
    xmat = xmaxmat.min(axis=1)
    nx = 10 if quick else 50
    # NOTE: as in main.m, the core loop uses xmaxmat(pop) (MATLAB linear
    # indexing = first column = bound set by juvenile survival), not xmat(pop).
    # Values of x that push another probability to 1 or above are discarded by
    # the validity check, so fewer than nx amplitudes are explored for most
    # species (only x = 0 for species 1). Kept to reproduce the article; use
    # xmat to explore the full range.
    xmaxloop = xmaxmat[:, 0]

    # %% Basic utilities and population statistics

    # Useful vectors and matrices to conduct matrix calculations
    one_ngroup_nparam = np.ones((ngroup, nparam))
    I_ngroupstage = np.eye(ngroup * nstage)

    # We choose a stage for which to compute the expected longevity and
    # lifetime reproductive output (here juveniles).
    # States are ordered group by group: [J1 A1 J2 A2 ... J16 A16].
    selectionvect = np.tile([True, False], ngroup)

    # Useful vectors and matrices to conduct calculations for the rewards
    t = nstage * ngroup   # number of transient states
    s = t + 1             # transient states + death (absorbing state)
    s1 = np.ones(s)
    t0 = np.zeros(t)
    Z = np.hstack([np.eye(t), t0.reshape(-1, 1)])  # removes the absorbing state

    # Getting basic information on the homogeneous populations with baseline vital rates.
    GTbasevect = np.full(npop, np.nan)
    initialstats = np.full((4, npop), np.nan)  # rows: mean L, var L, mean LRO, var LRO
    for pop in range(npop):
        if not np.any(np.isnan(thetamatpop[pop, :])):
            thetamat = thetamatpop[pop:pop + 1, :]
            Uhet, Fhet, Ahet = popHeteromatfact(thetamat, [1.0], 1, nstage)
            gt, _, _ = generation_time(Uhet, Fhet, Ahet)
            # NOTE: generation time is rounded to 2 decimals and this rounded
            # value is the target of the quality scenario.
            GTbasevect[pop] = np.round(gt, 2)

            # Ntilde is the fundamental matrix (see Caswell et al. 2018, doi:10.1002/ecm.1306)
            Ntilde = np.linalg.solve(np.eye(2) - Uhet, np.eye(2))
            # NOTE: [False, True] gives longevity statistics for ADULTS, whereas
            # LRO statistics ([True, False]) are for JUVENILES (as in main.m).
            # These baseline values are stored only.
            meanLE, varLE = lifetimeoutput(Ntilde, np.eye(2), np.array([False, True]))
            meanLRO, varLRO = reprooutputfactorial(Uhet, Fhet, Ntilde, np.array([True, False]),
                                                   np.ones(3), np.zeros(2),
                                                   np.hstack([np.eye(2), np.zeros((2, 1))]))
            initialstats[:, pop] = [meanLE[0], varLE[0], meanLRO[0], varLRO[0]]

    return dict(scenario=scenario, thetamatpop=thetamatpop, npop=npop, nparam=nparam,
                nstage=nstage, nclass=nclass, heterogeneitymat=heterogeneitymat,
                ngroup=ngroup, Rvect=Rvect, Rallcomb=Rallcomb, ncomb=ncomb, Px=Px, Py=Py,
                Rmat=Rmat, xmaxmat=xmaxmat, xmat=xmat, xmaxloop=xmaxloop, nx=nx,
                one_ngroup_nparam=one_ngroup_nparam, I_ngroupstage=I_ngroupstage,
                selectionvect=selectionvect, s1=s1, t0=t0, Z=Z,
                GTbasevect=GTbasevect, initialstats=initialstats)


def find_Y(f):
    """Equivalent of Y = fzero(f, [0 1]) with the error handling of main.m.
    Search Y in [0, 1] ("-" rates can only be lower than baseline).
    Options to pass to fzero function: brentq is used with tolerances at
    machine precision, like fzero's defaults. (In MATLAB,
    options = optimset('Display','iter') would show iterations.)"""
    with np.errstate(all="ignore"):
        fa, fb = f(0.0), f(1.0)
    if not (np.isfinite(fa) and np.isfinite(fb)):
        # At Y = 0 all "-" rates are zero. The dominant eigenvalue can then
        # belong to groups that never reproduce, so the generation time is
        # infinite (happens for species 2). Retry just above zero.
        # (In MATLAB v1.0.0, Y kept its value from the previous iteration here.)
        a = 1e-9
        with np.errstate(all="ignore"):
            fa = f(a)
        if not np.isfinite(fa):
            return np.nan
    else:
        a = 0.0
    if fa == 0:
        return a
    if fb == 0:
        return 1.0
    if np.sign(fa) == np.sign(fb):
        # Checks if a solution does not exist between 0 and 1
        return np.nan
    return brentq(f, a, 1.0, xtol=1e-15, rtol=8.9e-16, maxiter=500)


_P = None  # parameters shared with the worker processes


def _init_worker(scenario, quick):
    global _P
    _P = setup(scenario, quick)


def run_amplitude(task):
    """Core loop of main.m for ONE species (pop) and ONE amplitude (xind),
    over all birth distributions (comb). Runs in a worker process."""
    pop, xind = task
    P = _P
    ngroup, nstage, nclass, nparam = P["ngroup"], P["nstage"], P["nclass"], P["nparam"]
    Px, Py, ncomb = P["Px"], P["Py"], P["ncomb"]
    one_ngroup_nparam, I_ngroupstage = P["one_ngroup_nparam"], P["I_ngroupstage"]
    selectionvect = P["selectionvect"]
    shape4 = [nclass] * nparam  # 2x2x2x2 array of groups

    out = {k: np.full(ncomb, np.nan) for k in (
        "gap", "lambda1", "GT", "meanLE", "varLE", "IHeteroLE",
        "meanLRO", "varLRO", "IHeteroLRO")}
    out.update(SSD=np.full((ngroup * nstage, ncomb), np.nan),
               PI=np.full((ngroup, ncomb), np.nan),
               LElist=np.full((ngroup, ncomb), np.nan),
               LROlist=np.full((ngroup, ncomb), np.nan))

    theta = P["thetamatpop"][pop, :]
    thetamat = np.tile(theta, (ngroup, 1))

    # We set the amount of heterogeneity amplitude we need to apply
    xvect = np.linspace(0, P["xmaxloop"][pop], P["nx"])

    # Define an X value for the phenotype with better quality
    X = 1 + xvect[xind]
    Xmat = X * one_ngroup_nparam
    # thetax contains the vital rate values for high-quality individuals.
    thetax = thetamat * (Xmat ** Px)

    # We check that probabilities / rates are not greater than 1
    # or lower than 0. We ignore the 4th parameter since it is
    # fertility
    if np.any(np.array([1, 1, 1, 1]) - np.append(thetax[0, :3], 0) <= 1e-15) or np.any(thetax[0, :4] < 0):
        return pop, xind, False, out

    for comb in range(ncomb):
        Rmatsmall = P["Rmat"][:, :, comb]
        # Rcomb is the phenotype distribution at birth.
        Rcomb = Rmatsmall[:, 0] * Rmatsmall[:, 1] * Rmatsmall[:, 2] * Rmatsmall[:, 3]

        # Find the corresponding Y value to keep the population
        # stable or with the same generation time
        if P["scenario"] == "quality":
            f = lambda Y: minusparamfinder_ConstantG(thetax, Y, one_ngroup_nparam, Py, Rcomb,
                                                     P["GTbasevect"][pop], ngroup, nstage)
        else:
            f = lambda Y: minusparamfinder_ConstantLambda(thetax, Y, one_ngroup_nparam, Py, Rcomb)
        Y = find_Y(f)  # When X > 1

        if np.isnan(Y) or np.any(thetax[15, :4] < 0):
            continue
        Ymat = Y * one_ngroup_nparam
        thetaxy = thetax * (Ymat ** Py)

        # We construct the matrix population model
        Uhet, Fhet, Ahet = popHeteromatfact(thetaxy, Rcomb, ngroup, nstage)
        lambda1, SSD, GT = det_outputfactorial(Uhet, Fhet, Ahet)

        # Generate the fundamental matrix
        Ntilde = np.linalg.solve(I_ngroupstage - Uhet, I_ngroupstage)

        # Calculate stable distribution of the groups for a
        # particular stages
        PI = SSD[selectionvect]
        PI = PI / PI.sum()
        reshapedPI = PI.reshape(shape4, order="F")

        # Calculate mean and variance in longevity
        meanLE, varLE = lifetimeoutput(Ntilde, I_ngroupstage, selectionvect)

        # Variance decomposition in longevity
        LEwithinvar = PI @ varLE                 # Vw(L): individual stochasticity
        LEbetweenvartot = varcalc(meanLE, PI)    # Vb(L): individual heterogeneity
        IHeteroLE = LEbetweenvartot / (LEwithinvar + LEbetweenvartot)

        reshapedmeanLE = meanLE.reshape(shape4, order="F")
        LEbetweenvarlist = vardecomposition(reshapedmeanLE, reshapedPI, Px, Py, ngroup, P["heterogeneitymat"])

        # Calculate mean and variance in lifetime reproductive
        # success
        meanLRO, varLRO = reprooutputfactorial(Uhet, Fhet, Ntilde, selectionvect, P["s1"], P["t0"], P["Z"])

        # Variance decomposition in lifetime reproductive
        # success
        LROwithinvar = PI @ varLRO
        LRObetweenvartot = varcalc(meanLRO, PI)
        IHeteroLRO = LRObetweenvartot / (LROwithinvar + LRObetweenvartot)

        reshapedmeanLRO = meanLRO.reshape(shape4, order="F")
        LRObetweenvarlist = vardecomposition(reshapedmeanLRO, reshapedPI, Px, Py, ngroup, P["heterogeneitymat"])

        # Saving results
        # array_to_save_results_to  # nxn = dimension of the data that is saved
        out["gap"][comb] = X - Y                          # 1x1
        out["lambda1"][comb] = lambda1                    # 1x1
        out["SSD"][:, comb] = SSD                         # ngroup*nstage x 1
        out["PI"][:, comb] = PI                           # 1 x ngroup
        out["GT"][comb] = GT                              # 1x1
        out["meanLE"][comb] = PI @ meanLE                 # 1x1
        out["varLE"][comb] = LEbetweenvartot + LEwithinvar  # 1x1
        out["IHeteroLE"][comb] = IHeteroLE                # 1x1
        out["LElist"][:, comb] = LEbetweenvarlist         # 1 x variable
        out["meanLRO"][comb] = PI @ meanLRO               # 1x1
        out["varLRO"][comb] = LRObetweenvartot + LROwithinvar  # 1x1
        out["IHeteroLRO"][comb] = IHeteroLRO              # 1x1
        out["LROlist"][:, comb] = LRObetweenvarlist       # 1 x variable

    return pop, xind, True, out


def run_scenario(scenario, outdir, workers=None, quick=False):
    """Results storage, core loops and saving (sections of main.m)."""
    P = setup(scenario, quick)
    npop, nx, ncomb, ngroup, nstage = P["npop"], P["nx"], P["ncomb"], P["ngroup"], P["nstage"]

    # Results storage
    sortedPx = np.unique(P["Px"], axis=0)  # We get all unique combinations of heterogeneity groups
    nhet = sortedPx.shape[0]
    # Dimensions: comb (birth distribution), xind (amplitude), pop (species).
    # NaN = combination not computed (invalid or no solution for Y).
    gapvalresults = np.full((ncomb, nx, npop), np.nan)                  # 1x1
    lambda1valsresults = np.full((ncomb, nx, npop), np.nan)             # 1x1
    SSDresults = np.full((ngroup * nstage, 1, ncomb, nx, npop), np.nan)  # ngroupnstagex1
    PIresults = np.full((1, ngroup, ncomb, nx, npop), np.nan)           # 1xngroup
    GTresults = np.full((ncomb, nx, npop), np.nan)                      # 1x1

    meanLEresults = np.full((ncomb, nx, npop), np.nan)                  # 1x1
    varLEresults = np.full((ncomb, nx, npop), np.nan)                   # 1x1
    IHeteroLEresults = np.full((ncomb, nx, npop), np.nan)               # 1x1
    LEbetweenvarlistresults = np.full((ngroup, 1, ncomb, nx, npop), np.nan)  # 1xvariable

    meanLROresults = np.full((ncomb, nx, npop), np.nan)                 # 1x1
    varLROresults = np.full((ncomb, nx, npop), np.nan)                  # 1x1
    IHeteroLROresults = np.full((ncomb, nx, npop), np.nan)              # 1x1
    LRObetweenvarlistresults = np.full((ngroup, 1, ncomb, nx, npop), np.nan)  # 1xvariable

    # %% Core loops
    # Following code is written to reproduce results for QUALITY heterogeneity
    # If producing results for TRADE-OFFS, "minusparamfinder_ConstantG" is
    # replaced by "minusparamfinder_ConstantLambda" (in run_amplitude).
    # Please notice the arguments of these two functions are different.
    # The name of the save file also changes. Both are handled automatically
    # by the "scenario" argument.
    # One task per (species, amplitude); each task loops over all birth
    # distributions. Tasks are independent, so they are run in parallel.
    tasks = [(pop, xind) for pop in range(npop) for xind in range(nx)]
    tic = time.time()  # setting a timer
    nvalid = np.zeros(npop, dtype=int)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(scenario, quick)) as ex:
        for done, (pop, xind, valid, o) in enumerate(ex.map(run_amplitude, tasks, chunksize=1), 1):
            nvalid[pop] += valid
            gapvalresults[:, xind, pop] = o["gap"]
            lambda1valsresults[:, xind, pop] = o["lambda1"]
            SSDresults[:, 0, :, xind, pop] = o["SSD"]
            PIresults[0, :, :, xind, pop] = o["PI"]
            GTresults[:, xind, pop] = o["GT"]
            meanLEresults[:, xind, pop] = o["meanLE"]
            varLEresults[:, xind, pop] = o["varLE"]
            IHeteroLEresults[:, xind, pop] = o["IHeteroLE"]
            LEbetweenvarlistresults[:, 0, :, xind, pop] = o["LElist"]
            meanLROresults[:, xind, pop] = o["meanLRO"]
            varLROresults[:, xind, pop] = o["varLRO"]
            IHeteroLROresults[:, xind, pop] = o["IHeteroLRO"]
            LRObetweenvarlistresults[:, 0, :, xind, pop] = o["LROlist"]
            if done % nx == 0:
                print(f"  [{scenario}] {done}/{len(tasks)} amplitude steps done "
                      f"({time.time() - tic:.0f} s)", flush=True)
    toc = time.time() - tic

    # Amplitudes x explored for every species (one row per species)
    xvect = np.vstack([np.linspace(0, P["xmaxloop"][pop], nx) for pop in range(npop)])

    # %% EXECUTE THIS SECTION TO SAVE YOUR RESULTS
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, SCENARIOS[scenario])
    np.savez_compressed(
        base + ".npz",
        xvect=xvect, Rvect=P["Rvect"], Rallcomb=P["Rallcomb"], gapvalresults=gapvalresults,
        meanLROresults=meanLROresults, meanLEresults=meanLEresults, varLEresults=varLEresults,
        varLROresults=varLROresults, IHeteroLEresults=IHeteroLEresults,
        IHeteroLROresults=IHeteroLROresults, LEbetweenvarlistresults=LEbetweenvarlistresults,
        LRObetweenvarlistresults=LRObetweenvarlistresults, Px=P["Px"], Py=P["Py"],
        thetamatpop=P["thetamatpop"], lambda1valsresults=lambda1valsresults,
        GTresults=GTresults, SSDresults=SSDresults, PIresults=PIresults,
        GTbasevect=P["GTbasevect"], initialstats=P["initialstats"],
        heterogeneitymat=P["heterogeneitymat"])
    write_csv(base + ".csv.gz", P, xvect, gapvalresults, lambda1valsresults, GTresults,
              meanLEresults, varLEresults, IHeteroLEresults, LEbetweenvarlistresults,
              meanLROresults, varLROresults, IHeteroLROresults, LRObetweenvarlistresults)

    computed = int(np.sum(~np.isnan(IHeteroLEresults)))
    print(f"[{scenario}] {computed} populations computed in {toc:.0f} s; "
          f"amplitudes explored per species: {nvalid.tolist()}")
    print(f"[{scenario}] results saved to {base}.npz and {base}.csv.gz")
    return base


def write_csv(path, P, xvect, gap, lam, GT, mL, vL, KL, listL, mR, vR, KR, listR):
    """Long-format table: one row per computed population (non-proprietary,
    readable by any software). Column definitions are given in the README."""
    names = ["sj", "sa", "g", "f"]
    quality = P["scenario"] == "quality"
    # labels of the 16 factorial components (row order of Px)
    labels = ["x".join(n for n, b in zip(names, row) if b) or "none" for row in P["Px"].astype(int)]
    header = ["species", "G_baseline", "comb", "xind", "R_sj", "R_sa", "R_g", "R_f", "X", "Y",
              "lambda", "G", "mean_L", "var_L", "Vb_L", "Vw_L", "K_L",
              "mean_LRO", "var_LRO", "Vb_LRO", "Vw_LRO", "K_LRO"]
    if quality:
        header += [f"VbL_{l}" for l in labels[:15]] + [f"VbLRO_{l}" for l in labels[:15]]
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for pop in range(P["npop"]):
            for xind in range(P["nx"]):
                X = 1 + xvect[pop, xind]
                for comb in np.nonzero(~np.isnan(KL[:, xind, pop]))[0]:
                    VbL = vL[comb, xind, pop] * KL[comb, xind, pop]
                    VbR = vR[comb, xind, pop] * KR[comb, xind, pop]
                    row = [pop + 1, P["GTbasevect"][pop], comb, xind, *P["Rallcomb"][comb],
                           X, X - gap[comb, xind, pop], lam[comb, xind, pop], GT[comb, xind, pop],
                           mL[comb, xind, pop], vL[comb, xind, pop], VbL, vL[comb, xind, pop] - VbL,
                           KL[comb, xind, pop],
                           mR[comb, xind, pop], vR[comb, xind, pop], VbR, vR[comb, xind, pop] - VbR,
                           KR[comb, xind, pop]]
                    if quality:
                        row += list(listL[:15, 0, comb, xind, pop]) + list(listR[:15, 0, comb, xind, pop])
                    w.writerow([f"{v:.10g}" if isinstance(v, (float, np.floating)) else v for v in row])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenario", choices=list(SCENARIOS) + ["all"], default="all",
                        help="heterogeneity scenario to simulate (default: all three)")
    parser.add_argument("--workers", type=int, default=None,
                        help="number of parallel processes (default: all CPU cores)")
    parser.add_argument("--quick", action="store_true",
                        help="small parameter grid (3 birth probabilities, 10 amplitudes) "
                             "written to ./data/quick, to test the installation")
    parser.add_argument("--outdir", default=None, help="output folder (default: ./data)")
    args = parser.parse_args()

    outdir = args.outdir or os.path.join(HERE, "data", "quick" if args.quick else "")
    scenarios = list(SCENARIOS) if args.scenario == "all" else [args.scenario]
    for sc in scenarios:
        run_scenario(sc, outdir, workers=args.workers, quick=args.quick)


if __name__ == "__main__":
    main()
