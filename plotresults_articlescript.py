"""
plotresults_articlescript.py -- Figures and summary statistics of the article
=============================================================================

Python translation of plotresults_articlescript.m.

What this script does
---------------------
It reads the simulation results written by main.py in ./data and reproduces
the figures and the table of the Results section of the article:

  Section                                     Output (in ./figures)          Article
  ------------------------------------------  -----------------------------  --------
  1. LE PIE                                   LE_Pie.pdf                     Fig. 3a
  2. LRO PIE                                  LRO_Pie.pdf                    Fig. 3b
     (values of both pie figures)             fig3_pie_percentages.csv
  3. Vb(LRO) vs K(LRO) ("moving mean")        movmean_LRO.pdf                Fig. 2b
  4. Vb(LE)  vs K(LE)  ("moving mean")        movmean_LE.pdf                 Fig. 2a
     (values of both curves)                  fig2_movmean.csv.gz
  5. Quantiles of K for each species          table2_Kstats.csv              Table 2
  6. Correlation between K(L) and K(LRO)      printed to the console         Results 3.1

Every number shown in a figure or in Table 2 is also printed to the console
and written to a CSV file, so that the figures can be checked or redrawn with
any software. The columns of these CSV files are described in the README.

Usage
-----
    python plotresults_articlescript.py                    # reads ./data, writes ./figures
    python plotresults_articlescript.py --datadir D --figdir F

Inputs (written by main.py)
---------------------------
    data/simulation_results_Quality_Gconstant.npz            quality heterogeneity
    data/simulation_results_SaSjTradeoff_LAMBDAconstant.npz  trade-off sa / sj
    data/simulation_results_SaFTradeoff_LAMBDAconstant.npz   trade-off sa / F

Reminder on array dimensions in the result files (see README for details)
--------------------------------------------------------------------------
    IHeteroLEresults, IHeteroLROresults, varLEresults, varLROresults, ...
        shape (ncomb, nx, npop) = (birth distribution, amplitude, species)
    LEbetweenvarlistresults, LRObetweenvarlistresults
        shape (16, 1, ncomb, nx, npop): the 16 factorial components of Vb
    IHetero* contain K(H) = Vb(H) / V(H) as a PROPORTION (0-1); it is
    multiplied by 100 below to obtain percentages.
    NaN entries are populations that were not computed (no valid solution).

Notation
--------
    H    : life history outcome, L (longevity) or LRO (lifetime reproductive output)
    V(H) : total variance of H = Vw(H) + Vb(H)
    Vw   : within-group variance  (individual stochasticity)
    Vb   : between-group variance (fixed individual heterogeneity)
    K(H) : Vb(H) / V(H), share of the variance due to fixed heterogeneity
    G    : generation time of the species (years)
"""
import argparse
import csv
import gzip
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")  # no screen needed: figures are written to files
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


# =============================================================================
# %% SETTING UP GENERAL VARIABLES
# =============================================================================

def load(datadir, name):
    """Load one result file of main.py as a dictionary {variable name: array}."""
    path = os.path.join(datadir, f"simulation_results_{name}.npz")
    if not os.path.exists(path):
        raise SystemExit(f"Missing {path}: run 'python main.py' first.")
    return dict(np.load(path))


# Legend labels: generation time of the six species (Table 1a), fast to slow
legendlabels = ["G = 2.05", "G = 2.20", "G = 3.43", "G = 6.70", "G = 13.85", "G = 25.19"]

# Colours. viridis.m of the MATLAB version is replaced by matplotlib's built-in
# viridis colour map (same colours). The unused "flipud(hot)" colour map of the
# MATLAB script (which reverses the colorbar) is not needed.
colorvect = plt.cm.viridis(np.linspace(0, 1, 6))  # one colour per species (Fig. 2)
piecolors = plt.cm.viridis(np.linspace(0, 1, 5))  # colororder(viridis(5)): sj, sa, g, f, interactions (Fig. 3)
piecomponents = ["sj", "sa", "g", "f", "interactions"]

# Font sizes (points)
labelfontsize = 16  # axis labels
tickfontsize = 10   # tick labels

# Position of the factorial components in LE/LRObetweenvarlistresults.
# The 16 components follow the row order of Px (see main.py); in the quality
# scenario row i of Px lists the vital rates involved in component i:
#   Python index 7  -> Px row [1 0 0 0] -> sj alone (main effect)
#   Python index 11 -> Px row [0 1 0 0] -> sa alone
#   Python index 13 -> Px row [0 0 1 0] -> g alone
#   Python index 14 -> Px row [0 0 0 1] -> f alone
# (MATLAB rows 8, 12, 14, 15). All other rows are interactions between two or
# more vital rates; index 15 (Px row [0 0 0 0], "no factor") is always 0.
mainrows = [7, 11, 13, 14]
interactionrows = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15]

# Number of consecutive populations averaged in the moving mean of Fig. 2
movmeanwindow = 600


# =============================================================================
# Computations for Fig. 3 (pie charts)
# =============================================================================

def pie_shares(betweenvarlist, heterogeneitymat):
    """
    Mean contribution of heterogeneity in sj, sa, g, f and of all their
    interactions to the between-group variance Vb, for each species.

    Steps, for each species:
      1. take the 16 factorial components of Vb of every simulated population;
      2. express each component as a proportion of Vb of that population
         (so that every population has the same weight, whatever its Vb);
      3. average these proportions over all populations of the species;
      4. rescale the five averages (sj, sa, g, f, interactions) to sum to 1.

    INPUT
    betweenvarlist: LEbetweenvarlistresults or LRObetweenvarlistresults,
                    shape (16, 1, ncomb, nx, npop).
    heterogeneitymat: vital rates carrying heterogeneity (1) or not (0).
    OUTPUT
    shares: array (npop, 5), columns sj, sa, g, f, interactions (proportions).
    """
    npop = betweenvarlist.shape[-1]
    shares = np.full((npop, 5), np.nan)
    for j in range(npop):
        # 16 x (ncomb*nx) matrix: one column per simulated population
        # (order="F" reproduces MATLAB's column-major reshape)
        newmat = betweenvarlist[:, 0, :, :, j].reshape(16, -1, order="F").copy()

        # Remove tiny negative values produced by rounding errors in the
        # subtraction step of the factorial decomposition
        newmat[newmat < 0] = 0

        # Step 2: express each component as a proportion of Vb of that
        # population (the 16 components sum to Vb)
        with np.errstate(invalid="ignore", divide="ignore"):
            newmat = newmat / newmat.sum(axis=0, keepdims=True)
        newmat[np.isinf(newmat)] = np.nan  # populations with Vb = 0 are ignored

        # Main effects: share of Vb due to heterogeneity in each vital rate alone
        sjalone, saalone, galone, falone = (newmat[r] for r in mainrows)

        # Sum of all interaction terms. The division corrects for duplicated
        # rows when some vital rates carry no heterogeneity (it equals 1 when
        # heterogeneitymat = [1 1 1 1]). sum(..., "omitmissing") in MATLAB.
        interactions = np.nansum(newmat[interactionrows], axis=0) / 2 ** (4 - np.sum(heterogeneitymat))
        # populations that were not computed stay missing (not 0)
        interactions[np.all(np.isnan(newmat[interactionrows]), axis=0)] = np.nan

        # Step 3: mean proportion over all populations of the species
        X = np.array([np.nanmean(sjalone), np.nanmean(saalone), np.nanmean(galone),
                      np.nanmean(falone), np.nanmean(interactions)])

        # Step 4: rescale so that the five slices of the pie sum to 1
        if np.nansum(X) != 0 and not np.isnan(np.sum(X)):
            shares[j] = X / X.sum()
    return shares


def pie_figure(shares, title, path):
    """
    Draw one row of pie charts (one pie per species, fast to slow) and save it.
    INPUT
    shares: output of pie_shares(), array (npop, 5).
    title: figure title.
    path: output file (PDF).
    """
    fig, axes = plt.subplots(1, len(shares), figsize=(16, 3.2))
    for j, ax in enumerate(axes):
        if not np.any(np.isnan(shares[j])):
            # Slices in the order sj, sa, g, f, interactions. As MATLAB's
            # piechart, the first slice starts at 12 o'clock and slices run
            # clockwise. Black edges separate the slices.
            ax.pie(shares[j], colors=piecolors, startangle=90, counterclock=False,
                   wedgeprops=dict(edgecolor="k", linewidth=1.5))
        ax.set_title(legendlabels[j], fontsize=11)  # generation time of the species

    # Common legend below the pies (one coloured square per slice)
    handles = [matplotlib.patches.Patch(facecolor=piecolors[i], edgecolor="k", label=l)
               for i, l in enumerate([r"$s_j$", r"$s_a$", r"$\gamma$", r"$F$", r"$i$ (interactions)"])]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=11)
    fig.suptitle(title, fontsize=12)
    fig.subplots_adjust(bottom=0.18, top=0.80)  # leave room for the legend and the title
    fig.savefig(path)
    plt.close(fig)


# =============================================================================
# Computations for Fig. 2 ("moving mean" of K as a function of Vb)
# =============================================================================

def movmean_matlab(y, k=movmeanwindow):
    """
    Moving mean, identical to MATLAB's movmean(y, k) for an even window k:
      - the window of element i covers elements i-k/2 to i+k/2-1;
      - near both ends the window is truncated (it shrinks);
      - if the window contains a NaN, the result is NaN ('includenan').
    Implemented with cumulative sums, so that it is fast for long vectors.
    """
    n = len(y)
    kb, kf = k // 2, k // 2 - 1                         # elements before / after i
    c = np.concatenate([[0.0], np.cumsum(np.nan_to_num(y))])   # cumulative sum of values
    cn = np.concatenate([[0], np.cumsum(np.isnan(y))])         # cumulative number of NaN
    i = np.arange(n)
    lo, hi = np.maximum(0, i - kb), np.minimum(n, i + kf + 1)  # window [lo, hi) of each element
    m = (c[hi] - c[lo]) / (hi - lo)                     # mean over the window
    m[(cn[hi] - cn[lo]) > 0] = np.nan                   # NaN if the window contains a NaN
    return m


def movmean_curves(varresults, IHetero):
    """
    Data of one panel of Fig. 2, for each species.
    Steps, for each species:
      1. compute Vb = V * K for every simulated population (because K = Vb / V);
      2. sort all populations of the species by increasing Vb
         (populations that were not computed, NaN, are put at the end);
      3. compute the moving mean of K (in %) over 600 consecutive populations.
    The curve therefore shows how the share of variance due to heterogeneity
    changes with the amplitude of heterogeneity (measured by Vb).
    OUTPUT
    list of (sortedx, ymovmean) per species: Vb sorted, moving mean of K in %.
    """
    # Processing data
    testmat = varresults * IHetero  # Step 1: Vb = V * K
    curves = []
    for j in range(IHetero.shape[2]):
        # All populations of species j as vectors (column-major, as MATLAB's (:))
        xvalues = testmat[:, :, j].ravel(order="F")
        yvalues = IHetero[:, :, j].ravel(order="F")

        # Step 2: sort by Vb (NaN values are sorted last, as in MATLAB)
        sortedxindexes = np.argsort(np.where(np.isnan(xvalues), np.inf, xvalues), kind="stable")
        sortedx = xvalues[sortedxindexes]
        sortedy = yvalues[sortedxindexes]

        # Step 3: computing the 600 point average, in %
        ymovmean = 100 * movmean_matlab(sortedy, movmeanwindow)
        curves.append((sortedx, ymovmean))
    return curves


def movmean_figure(curves, ymax, label, path):
    """
    Draw one panel of Fig. 2 and save it.
    x-axis: between-group variance Vb(H) (amplitude of heterogeneity), log scale.
    y-axis: K(H), percentage of the total variance due to heterogeneity.
    One curve per species (colour = generation time).
    INPUT
    curves: output of movmean_curves().
    ymax: upper limit of the y-axis (same for both panels).
    label: "L" or "LRO" (used in the axis labels).
    path: output file (PDF).
    """
    fig, ax = plt.subplots(figsize=(5, 5))  # Creating plot grid (one panel)

    # plotting the moving average: one line per species
    for j, (sortedx, ymovmean) in enumerate(curves):
        ax.plot(sortedx, ymovmean, linestyle="-", color=colorvect[j], marker=".",
                markersize=2, linewidth=2, label=legendlabels[j])

    # Formatting axes
    for spine in ax.spines.values():
        spine.set_linewidth(2)  # Increases linewidth for all axes
    ax.tick_params(width=2, length=6)  # Increases tick length for all axes
    ax.set_xscale("log")  # Vb spans many orders of magnitude
    # x-range: smallest positive to largest Vb over all species
    allx = np.concatenate([c[0] for c in curves])
    xmin, xmax = np.nanmin(allx[allx > 0]), np.nanmax(allx)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, ymax)
    # One tick per power of ten, labelled 10^k
    ax.set_xticks(10.0 ** np.arange(np.floor(np.log10(xmin)), np.ceil(np.log10(xmax)) + 1))

    # Setting fonts
    ax.tick_params(axis="x", labelrotation=45, labelsize=tickfontsize)  # changes fontsize for label and ticks
    ax.tick_params(axis="y", labelsize=tickfontsize)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)  # keep only the left and bottom axes

    # Formatting labels
    ax.set_xlabel(rf"$V_b({label})$", fontsize=labelfontsize)  # changes fontsize for label only
    ax.set_ylabel(rf"$K({label})$", fontsize=labelfontsize)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# =============================================================================
# Computations for Table 2
# =============================================================================

def Kstats(IHetero):
    """
    Distribution of K(H) (in %) over all simulated populations of each species
    (all birth distributions x all amplitudes).
    OUTPUT
    array (npop, 7), columns: number of populations, D1 (10th percentile),
    Q1 (25th), median, Q3 (75th), D9 (90th percentile), max.
    Percentiles use the same definition as MATLAB's quantile()
    (numpy method "hazen"), so that values are identical to the MATLAB version.
    NOTE: the submitted version of Table 2 labelled the 25th and 75th
    percentiles as D1 and D9; both quartiles and deciles are now reported.
    """
    rows = []
    for j in range(IHetero.shape[2]):
        k = 100 * IHetero[:, :, j].ravel()  # K in % for all populations of species j
        k = k[~np.isnan(k)]                 # keep computed populations only
        if len(k) == 0:
            rows.append([0] + [np.nan] * 6)
            continue
        q = np.quantile(k, [0.10, 0.25, 0.50, 0.75, 0.90], method="hazen")
        rows.append([len(k), *q, k.max()])
    return np.array(rows)


# =============================================================================
# Main script
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datadir", default=os.path.join(HERE, "data"), help="folder with the results of main.py")
    parser.add_argument("--figdir", default=os.path.join(HERE, "figures"), help="output folder")
    args = parser.parse_args()
    os.makedirs(args.figdir, exist_ok=True)

    # Quality heterogeneity results (used for Figs. 2 and 3 and Table 2)
    q = load(args.datadir, "Quality_Gconstant")
    G = [f"{g:.2f}" for g in q["GTbasevect"]]  # generation times as text, for tables

    # -------------------------------------------------------------------------
    # %% 1. LE PIE (Fig. 3a): factorial decomposition of Vb(L)
    # -------------------------------------------------------------------------
    sharesL = pie_shares(q["LEbetweenvarlistresults"], q["heterogeneitymat"])
    pie_figure(sharesL, r"Partition of $V_b(L)$", os.path.join(args.figdir, "LE_Pie.pdf"))

    # -------------------------------------------------------------------------
    # %% 2. LRO PIE (Fig. 3b): factorial decomposition of Vb(LRO)
    # -------------------------------------------------------------------------
    sharesLRO = pie_shares(q["LRObetweenvarlistresults"], q["heterogeneitymat"])
    pie_figure(sharesLRO, r"Partition of $V_b(LRO)$", os.path.join(args.figdir, "LRO_Pie.pdf"))

    # Values of the pie slices: console and CSV file (in %)
    print("Mean contribution to the between-group variance (%)  [Fig. 3]")
    print("            G      sj      sa       g       f   inter.")
    with open(os.path.join(args.figdir, "fig3_pie_percentages.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["outcome", "G", *[f"share_{c}_percent" for c in piecomponents]])
        for lab, S in (("L", sharesL), ("LRO", sharesLRO)):
            for j in range(len(S)):
                print(f"  {lab:>3} {G[j]:>6} " + " ".join(f"{100 * v:7.2f}" for v in S[j]))
                w.writerow([lab, G[j], *[f"{100 * v:.4f}" for v in S[j]]])

    # -------------------------------------------------------------------------
    # %% 3. CODE for exploring the relantionship between Vb(LRO) and K(LRO) (Fig. 2b)
    # %% 4. CODE for exploring the relantionship between Vb(LE) and K(LE) (Fig. 2a)
    # Code used to represent the "moving mean" for all six populations, in a
    # single graph. Both panels share the same y-range so that they are
    # directly comparable.
    # -------------------------------------------------------------------------
    ymax = 105 * max(np.nanmax(q["IHeteroLEresults"]), np.nanmax(q["IHeteroLROresults"]))
    curvesLRO = movmean_curves(q["varLROresults"], q["IHeteroLROresults"])
    curvesL = movmean_curves(q["varLEresults"], q["IHeteroLEresults"])
    movmean_figure(curvesLRO, ymax, "LRO", os.path.join(args.figdir, "movmean_LRO.pdf"))
    movmean_figure(curvesL, ymax, "L", os.path.join(args.figdir, "movmean_LE.pdf"))

    # Values of the curves of Fig. 2 (computed populations only), compressed CSV
    with gzip.open(os.path.join(args.figdir, "fig2_movmean.csv.gz"), "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["outcome", "species", "G", "rank", "Vb", "K_movmean_percent"])
        for lab, curves in (("L", curvesL), ("LRO", curvesLRO)):
            for j, (x, y) in enumerate(curves):
                keep = ~np.isnan(x) & ~np.isnan(y)
                for rank in np.nonzero(keep)[0]:
                    w.writerow([lab, j + 1, G[j], rank + 1, f"{x[rank]:.10g}", f"{y[rank]:.10g}"])

    # -------------------------------------------------------------------------
    # %% 5. CODE for getting the quantiles, extrema and mean of K(L/LRO) for
    #       each population (QUALITY)  -> Table 2, top section
    # -------------------------------------------------------------------------
    table = [("Quality", "K(L)", q["IHeteroLEresults"]), ("Quality", "K(LRO)", q["IHeteroLROresults"])]

    # %% CODE for getting the quantiles, extrema and mean of K(L) for each
    #    population (TRADEOFFS)  -> Table 2, middle and bottom sections.
    # Trade-off groups have the same fitness, so K(LRO) is 0 and only K(L) is
    # reported for them.
    for scen, name in (("Trade-off sa/sj", "SaSjTradeoff_LAMBDAconstant"),
                       ("Trade-off sa/F", "SaFTradeoff_LAMBDAconstant")):
        table.append((scen, "K(L)", load(args.datadir, name)["IHeteroLEresults"]))

    # Console (populations are stored "vertically": one row per species) and CSV
    print("\nDistribution of K(H) in %  [Table 2]")
    with open(os.path.join(args.figdir, "table2_Kstats.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["heterogeneity", "outcome", "G", "n_populations",
                    "D1_percent", "Q1_percent", "median_percent", "Q3_percent", "D9_percent", "max_percent"])
        for scen, outcome, K in table:
            print(f"\n  {scen} -- {outcome}\n       G       n       D1       Q1   median       Q3       D9      max")
            for j, r in enumerate(Kstats(K)):
                print(f"  {G[j]:>6} {int(r[0]):>7} " + " ".join(f"{v:8.3f}" for v in r[1:]))
                w.writerow([scen, outcome, G[j], int(r[0]), *[f"{v:.6g}" for v in r[1:]]])

    # -------------------------------------------------------------------------
    # %% 6. CORRELATION BETWEEN K(L) AND K(LRO) (quality heterogeneity)
    # Pearson correlation over all simulated populations of all species.
    # -------------------------------------------------------------------------
    A_vec = q["IHeteroLEresults"].ravel()
    B_vec = q["IHeteroLROresults"].ravel()
    mask = ~np.isnan(A_vec) & ~np.isnan(B_vec)  # keep only valid pairs
    r = np.corrcoef(A_vec[mask], B_vec[mask])[0, 1]
    print(f"\nPearson correlation coefficient between K(L) and K(LRO): {r:.3f}")
    print(f"\nFigures and tables written to {args.figdir}")


if __name__ == "__main__":
    main()
