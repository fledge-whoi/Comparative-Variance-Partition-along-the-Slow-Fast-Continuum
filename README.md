This README file was generated on 2026-07-01 and last updated on 2026-10-16.

# GENERAL INFORMATION

**Project name:** Comparative Variance Partition along the Slow-Fast Continuum
**Version:** v2.0.0 (Python translation of the MATLAB version 1.0.0)
**Associated article:** Elise J.-P., Fay R., Garnier J., Caswell H., Gaillard J.-M., Gamelon M., Aubry L. M., Jenouvrier S. *Differential Contributions of Heterogeneity and Stochasticity to Individuals' Life Histories along the Slow-fast Continuum.* Manuscript submitted to Ecology Letters.
**Archive (code and data):** Zenodo, https://doi.org/10.5281/zenodo.23298484 (previous MATLAB version 1.0.0: https://doi.org/10.5281/zenodo.21210108)
**Code repository:** https://github.com/jp-elise/individual_heterogeneity_python
**Short description:** Python code and simulated data for the article. The code simulates heterogeneous populations of six synthetic species along the slow-fast continuum and computes the mean and variance of longevity and lifetime reproductive output. It partitions that variance into individual stochasticity (within-group variance) and fixed individual heterogeneity (between-group variance), and decomposes the latter factorially among vital rates.

**Date of creation:** 2023-02-15 to 2026-10-10
**Languages and formats:** all code is in Python, and all data are in open formats (CSV and NumPy `.npz`). No proprietary software is needed.

___

# PROJECT OVERVIEW

The life cycle is a two-stage (juvenile, adult) matrix population model with four vital rates: juvenile survival `sj`, adult survival `sa`, maturation rate `g` (gamma in the article) and fertility `f` (F in the article). Each species is split into fixed heterogeneity groups. The vital rates of a group are multiplied by a coefficient `X` ("+" value; Q+ or T+ in the article) or `Y` ("-" value; Q- or T- in the article). Three heterogeneity scenarios are analysed:

| Scenario (`--scenario`) | Groups | Constraint used to solve for `Y` |
|---|---|---|
| `quality` | 16 groups: every "+"/"-" combination on the four vital rates | same generation time as the baseline species |
| `tradeoff_sasj` | trade-off group (sa × X, sj × Y) + baseline group | stable population (lambda = 1) |
| `tradeoff_saf` | trade-off group (f × X, sa × Y) + baseline group | stable population (lambda = 1) |

The code covers every species, every amplitude `X` and every distribution of newborns among groups. For each population it computes:
- longevity and lifetime reproductive output (LRO) statistics with Markov chain methods (Caswell 2009; van Daalen & Caswell 2017, 2020);
- the variance ratio K(H) = Vb / (Vb + Vw);
- the factorial decomposition of Vb (Caswell 2025).

## Files

```
.
├── main.py                          Step 1 - simulations: writes the result files to ./data
├── plotresults_articlescript.py     Step 2 - Figures 2 and 3 and Table 2 from ./data -> ./figures
├── functions/                       One module per function (same names as the MATLAB version)
│   ├── popHeteromatfact.py          Builds U, F and A = U + F for a heterogeneous population
│   ├── det_outputfactorial.py       Growth rate, stable stage distribution, generation time
│   ├── lifetimeoutput.py            Mean and variance of longevity (fundamental matrix)
│   ├── reprooutputfactorial.py      Mean and variance of LRO (Markov chain with rewards)
│   ├── varcalc.py                   Between-group variance of group means
│   ├── vardecomposition.py          Factorial decomposition of the between-group variance
│   ├── minusparamfinder_ConstantG.py       Objective for Y: same generation time (quality)
│   ├── minusparamfinder_ConstantLambda.py  Objective for Y: lambda = 1 (trade-offs)
│   └── combinations.py              Equivalent of MATLAB's combinations()
├── tests/
│   ├── test_functions.py            Unit tests against analytical results
│   └── test_reproduce_article.py    Checks that ./data reproduces Table 2
├── data/                            Simulation results (see DATA)
├── figures/                         Figures, plus the values shown in them (see DATA)
├── requirements.txt                 Python dependencies
├── README.md                        This file
├── CITATION.cff                     How to cite the code and data
├── LICENSE                          Licence of the code (GPL-3.0-or-later)
└── DATA_LICENSE.md                  Licence of the data (CC BY 4.0)
```

All code is commented. Every script and function starts with a header describing what it does, its inputs and its outputs. The plotting script explains, step by step, how each figure and table is computed from the data.

___

# INSTALLATION

**Software:** Python 3.9 or later. All software used is free and open source.
**Dependencies:** numpy (≥ 1.24), scipy (≥ 1.10) and matplotlib (≥ 3.7). No other package is needed; the tests use Python's built-in `unittest`.
**Tested with:** Python 3.13.16, numpy 2.5.3, scipy 1.18.1 and matplotlib 3.11.2 on Ubuntu 24.04 (x86_64). The code uses only standard, platform-independent features of these packages, so it should also run on macOS and Windows, but it has not been tested there.

```bash
python -m pip install -r requirements.txt
```

**Setup:** none. Unzip the archive and run the commands below from the folder that contains `main.py`. On macOS or Linux, use `python3` instead of `python` if needed.

___

# USAGE

## Reproducing the results of the article

```bash
python main.py --quick                     # optional: quick test of the installation (< 1 min)
python main.py                             # step 1: simulations, the three scenarios -> ./data
python plotresults_articlescript.py        # step 2: Figures 2-3 and Table 2 -> ./figures
python -m unittest discover -s tests -v    # optional: tests (incl. comparison with Table 2)
```

- **Run time of step 1:** about 17 minutes in total (quality 8 min, each trade-off scenario 4–5 min) on 2 CPU cores. It scales roughly with the number of cores.
- **Run time of step 2:** about 10 seconds.
- **Expected output of step 2:** the console shows the values of Fig. 3, Table 2 and the correlation between K(L) and K(LRO) (Pearson r = 0.938). These numbers are also written to CSV files in `./figures` (see DATA).
- **Without re-running the simulations:** download the files of `./data` from Zenodo, put them in a folder `data` next to `main.py`, and run step 2 only.

## Archive layout on Zenodo

Zenodo does not support folders, so the project is split into three archives:

| File on Zenodo | Content |
|---|---|
| `code.zip` | The Python project: scripts, `functions/`, `tests/`, documentation and licences |
| `data.zip` | A folder `data/` with the simulation results (3 `.csv.gz` and 3 `.npz` files) |
| `figures.zip` | A folder `figures/` with Figures 2 and 3 (PDF) and the values shown in them (CSV) |

`README.md`, `CITATION.cff`, `LICENSE` and `DATA_LICENSE.md` are also available as separate files on the record.

**Unzip the three archives into the same folder.** This recreates the structure shown in "Files" above, with `./data` and `./figures` next to `main.py`:

```bash
unzip code.zip -d project && unzip data.zip -d project && unzip figures.zip -d project
cd project
```

To reproduce the figures and Table 2 without re-running the simulations, `code.zip` and `data.zip` are enough: run `python plotresults_articlescript.py`.

**Options of `main.py`:**
- `--scenario quality|tradeoff_sasj|tradeoff_saf|all`: scenario(s) to simulate (default `all`).
- `--workers N`: number of CPU cores to use (default: all).
- `--quick`: small grid (3 birth probabilities × 10 amplitudes) written to `./data/quick`. The figures from this grid are only meant as a test.
- `--outdir PATH`: output folder.

**Options of `plotresults_articlescript.py`:** `--datadir PATH` (input folder) and `--figdir PATH` (output folder).

## Correspondence between outputs and the article

| Article | File in `./figures` | Values in |
|---|---|---|
| Fig. 2a (K(L) as a function of Vb(L)) | `movmean_LE.pdf` | `fig2_movmean.csv.gz` (outcome = L) |
| Fig. 2b (K(LRO) as a function of Vb(LRO)) | `movmean_LRO.pdf` | `fig2_movmean.csv.gz` (outcome = LRO) |
| Fig. 3a (partition of Vb(L)) | `LE_Pie.pdf` | `fig3_pie_percentages.csv` (outcome = L) |
| Fig. 3b (partition of Vb(LRO)) | `LRO_Pie.pdf` | `fig3_pie_percentages.csv` (outcome = LRO) |
| Table 2 | – | `table2_Kstats.csv` |
| Table 1a (vital rates of the species) | – | `thetamatpop` and `GTbasevect` in any `.npz` file |

## Reusing the code

- Species: edit `thetamatpop` in `setup()` of `main.py` (one row per species: `sj sa g f`).
- Heterogeneity structure: edit `Pxunfiltered`, `Pyunfiltered` and `heterogeneitymat` (see the comments in `main.py`).
- Distribution of newborns among groups: edit `Rvect`.
- Amplitudes of heterogeneity: edit `nx` (number of values) and `xmaxloop` (upper bound).
- Other life cycles: change `nstage`, the vital rate checks in `run_amplitude()` and the matrix construction in `functions/popHeteromatfact.py` (see Caswell et al. 2018, doi:10.1002/ecm.1306). The code currently supports two values ("+" and "-") per vital rate.

___

# DATA

**Origin:** all data are **simulated** with `main.py` (simulation results) and summarised with `plotresults_articlescript.py` (figure and table values). No empirical data are used. Running the two scripts regenerates every file exactly.
**Where:** the files are archived on Zenodo (https://doi.org/10.5281/zenodo.23298484) and are written to `./data` and `./figures` when the code is run.
**Licence:** CC BY 4.0 (see LICENSES).

## Overview of the data files

| File | Format | Size | Content |
|---|---|---|---|
| `data/simulation_results_Quality_Gconstant.csv.gz` | compressed CSV | ~38 MB | One row per simulated population, quality heterogeneity (158,901 rows) |
| `data/simulation_results_SaSjTradeoff_LAMBDAconstant.csv.gz` | compressed CSV | ~3 MB | Same, trade-off between sa and sj (439,377 rows) |
| `data/simulation_results_SaFTradeoff_LAMBDAconstant.csv.gz` | compressed CSV | ~3 MB | Same, trade-off between sa and F (480,140 rows) |
| `data/simulation_results_<scenario>.npz` (3 files) | NumPy archive | 100–150 MB each | Same results as full arrays; input of `plotresults_articlescript.py` |
| `figures/fig2_movmean.csv.gz` | compressed CSV | ~4 MB | Values of the curves of Fig. 2 |
| `figures/fig3_pie_percentages.csv` | CSV | < 1 kB | Values of the pie charts of Fig. 3 |
| `figures/table2_Kstats.csv` | CSV | < 3 kB | Table 2 |

**Reading the files:**
- **CSV files** are comma-separated, with column names in the first row, a dot as the decimal separator, and one row per record. `.csv.gz` files are gzip-compressed CSV files. They open directly in Python (`pandas.read_csv("file.csv.gz")`) or R (`read.csv("file.csv.gz")`), or in a spreadsheet after decompression with any archive tool.
- **`.npz` files** are open NumPy archives. Read them with `numpy.load("file.npz")`, which returns a dictionary of arrays.
- **Missing values:** in the `.npz` arrays, `NaN` marks a population that was not computed. This happens when the amplitude makes a probability ≥ 1, or when no `Y` in [0, 1] satisfies the constraint. The CSV files contain computed populations only, so they have no missing values.

**Units used in all files:**
- One time step is one year. Longevity and generation time are in years.
- Lifetime reproductive output (LRO) and fertility are numbers of recruits, i.e. offspring alive at the next pre-breeding census.
- Variances are in squared units (years² or recruits²).
- Probabilities, proportions, coefficients and K are dimensionless. Columns whose name ends in `_percent` are percentages.

**Indices:**
- Species are numbered 1 to 6 in the CSV files, as in the article (1 = fastest, G = 2.05; 6 = slowest, G = 25.19). In the `.npz` arrays they are indexed 0 to 5.
- `comb` and `xind` start at 0 in all files.

## 1. Simulation results: `data/simulation_results_<scenario>.csv.gz`

One row per simulated population, i.e. one combination of species × amplitude of heterogeneity × distribution of newborns among groups.

| Column | Description | Units |
|---|---|---|
| `species` | Species number (1 = fastest to 6 = slowest; Table 1a) | – |
| `G_baseline` | Generation time of the baseline (homogeneous) species | years |
| `comb` | Index of the distribution of newborns among groups (0–2400; row of `Rallcomb` in the `.npz` file) | – |
| `xind` | Index of the amplitude of heterogeneity (0–49; column of `xvect` in the `.npz` file) | – |
| `R_sj`, `R_sa`, `R_g`, `R_f` | Probability of being born with the "+" value of sj, sa, g and f. In the trade-off scenarios only `R_sj` matters: it is the probability of being born in the trade-off group | probability |
| `X` | Coefficient multiplying the "+" vital rates (Q+ or T+ in the article) | dimensionless |
| `Y` | Coefficient multiplying the "-" vital rates (Q- or T- in the article) | dimensionless |
| `lambda` | Asymptotic growth rate of the heterogeneous population | per year |
| `G` | Generation time of the heterogeneous population | years |
| `mean_L` | Mean longevity of newborns | years |
| `var_L` | Total variance in longevity, V(L) = Vb(L) + Vw(L) | years² |
| `Vb_L` | Between-group variance in longevity (fixed individual heterogeneity) | years² |
| `Vw_L` | Within-group variance in longevity (individual stochasticity) | years² |
| `K_L` | Share of the variance in longevity due to heterogeneity, Vb_L / var_L | proportion (0–1) |
| `mean_LRO`, `var_LRO`, `Vb_LRO`, `Vw_LRO`, `K_LRO` | Same as above for lifetime reproductive output | recruits, recruits², recruits², recruits², proportion |
| `VbL_<rates>` | *Quality scenario only.* Factorial component of `Vb_L` due to heterogeneity in the vital rates listed in `<rates>`. One rate gives a main effect (e.g. `VbL_sa`); several rates joined by `x` give their interaction (e.g. `VbL_sjxsa`, `VbL_sjxsaxgxf`). The 15 columns sum to `Vb_L` | years² |
| `VbLRO_<rates>` | *Quality scenario only.* Same for `Vb_LRO` | recruits² |

## 2. Simulation results as arrays: `data/simulation_results_<scenario>.npz`

These contain the same results as the CSV files, stored as arrays with the same names and dimensions as in the MATLAB version. Most arrays have three dimensions, `(comb, xind, pop)` = (2401 distributions of newborns, 50 amplitudes, 6 species).

| Variable | Shape | Description | Units |
|---|---|---|---|
| `thetamatpop` | 6 × 4 | Baseline vital rates of the six species; columns `sj sa g f` (Table 1a) | probabilities; `f` in recruits per adult per year |
| `GTbasevect` | 6 | Generation time of each baseline species, rounded to 2 decimals | years |
| `initialstats` | 4 × 6 | Baseline species statistics. Rows: mean and variance of longevity (both from the ADULT stage), mean and variance of LRO (both from the juvenile stage) | years, years², recruits, recruits² |
| `Rvect` | 7 | Candidate probabilities of being born with the "+" value of a vital rate | probability |
| `Rallcomb` | 2401 × 4 | All combinations of `Rvect` for `[Rsj Rsa Rg Rf]` (index `comb`) | probability |
| `xvect` | 6 × 50 | Amplitudes `x` (`X = 1 + x`) for each species (index `xind`) | dimensionless |
| `Px` | 16 × 4 | 1 if the vital rate (column `sj sa g f`) of the group (row) is multiplied by `X` | 0/1 |
| `Py` | 16 × 4 | 1 if the vital rate of the group is multiplied by `Y` | 0/1 |
| `heterogeneitymat` | 4 | 1 if heterogeneity is introduced on the vital rate (`sj sa g f`) | 0/1 |
| `gapvalresults` | 2401 × 50 × 6 | `X - Y`, the gap between the "+" and "-" coefficients | dimensionless |
| `lambda1valsresults` | 2401 × 50 × 6 | Population growth rate lambda | per year |
| `GTresults` | 2401 × 50 × 6 | Generation time of the heterogeneous population | years |
| `SSDresults` | 32 × 1 × 2401 × 50 × 6 | Stable stage × group distribution (sums to 1); states ordered `[J1 A1 J2 A2 ... J16 A16]` | proportion |
| `PIresults` | 1 × 16 × 2401 × 50 × 6 | Proportion of each group among juveniles at the stable distribution; the weights used in the variance partition | proportion |
| `meanLEresults` | 2401 × 50 × 6 | Mean longevity of newborns | years |
| `varLEresults` | 2401 × 50 × 6 | Total variance in longevity, V(L) | years² |
| `IHeteroLEresults` | 2401 × 50 × 6 | K(L) = Vb(L) / V(L) | proportion |
| `LEbetweenvarlistresults` | 16 × 1 × 2401 × 50 × 6 | Factorial components of Vb(L), in the row order of `Px` (see below) | years² |
| `meanLROresults` | 2401 × 50 × 6 | Mean LRO of newborns | recruits |
| `varLROresults` | 2401 × 50 × 6 | Total variance in LRO, V(LRO) | recruits² |
| `IHeteroLROresults` | 2401 × 50 × 6 | K(LRO) = Vb(LRO) / V(LRO) | proportion |
| `LRObetweenvarlistresults` | 16 × 1 × 2401 × 50 × 6 | Factorial components of Vb(LRO) | recruits² |

- **Groups** (16) are ordered as the rows of `Px`. In the quality scenario, index 0 is `[+ + + +]` and index 15 is `[- - - -]` for `[sj sa g f]`.
- **Factorial components:** indices 7, 11, 13 and 14 are the main effects of `sj`, `sa`, `g` and `f`. Every other index is the interaction of the rates marked 1 in the same row of `Px`. Index 15 is always 0. The components are not meaningful in the trade-off scenarios.
- **Derived quantities:** Vb = `varLEresults * IHeteroLEresults` and Vw = `varLEresults - Vb` (likewise for LRO).

## 3. Values of Fig. 2: `figures/fig2_movmean.csv.gz`

For each species, all simulated populations of the quality scenario are sorted by increasing between-group variance. K is then averaged over a moving window of 600 consecutive populations. Each row is one point of a curve.

| Column | Description | Units |
|---|---|---|
| `outcome` | `L` (longevity, Fig. 2a) or `LRO` (lifetime reproductive output, Fig. 2b) | – |
| `species` | Species number (1–6) | – |
| `G` | Generation time of the species (legend of the figure) | years |
| `rank` | Rank of the population after sorting by `Vb` (1 = smallest) | – |
| `Vb` | Between-group variance of the population (x-axis) | years² (L) or recruits² (LRO) |
| `K_movmean_percent` | Moving mean of K over 600 populations centred on this one (y-axis) | % |

## 4. Values of Fig. 3: `figures/fig3_pie_percentages.csv`

Mean contribution of each source of fixed heterogeneity to the between-group variance (quality scenario). For each population, each factorial component is divided by Vb. These proportions are averaged over all populations of the species, then rescaled to sum to 100%.

| Column | Description | Units |
|---|---|---|
| `outcome` | `L` (Fig. 3a) or `LRO` (Fig. 3b) | – |
| `G` | Generation time of the species | years |
| `share_sj_percent`, `share_sa_percent`, `share_g_percent`, `share_f_percent` | Main effect of heterogeneity in juvenile survival, adult survival, maturation rate and fertility | % |
| `share_interactions_percent` | Sum of all interaction terms | % |

## 5. Table 2: `figures/table2_Kstats.csv`

Distribution of K across all simulated populations of each species. Percentiles follow MATLAB's definition of `quantile` (numpy method `"hazen"`).

| Column | Description | Units |
|---|---|---|
| `heterogeneity` | `Quality`, `Trade-off sa/sj` or `Trade-off sa/F` | – |
| `outcome` | `K(L)` or `K(LRO)` | – |
| `G` | Generation time of the species | years |
| `n_populations` | Number of simulated populations for this species | – |
| `D1_percent`, `D9_percent` | 10th and 90th percentiles (first and ninth deciles) of K | % |
| `Q1_percent`, `Q3_percent` | 25th and 75th percentiles (quartiles) of K | % |
| `median_percent`, `max_percent` | Median and maximum of K | % |

___

# CHANGES FROM THE MATLAB VERSION 1.0.0

- **Language:** the code is translated from MATLAB to Python, so no proprietary software is needed. All computations are kept, and results were checked against the MATLAB version (see `tests/`). File and function names are unchanged, and the comments of the MATLAB code are kept at the corresponding places.
- **Root finding:** MATLAB's `fzero` is replaced by `scipy.optimize.brentq` on the same interval [0, 1], with tolerances at machine precision.
- **Fix:** when the objective function is infinite at Y = 0, the search is now repeated on [1e-9, 1]. In version 1.0.0, `Y` instead kept its value from the previous iteration. This affects species 2 (G = 2.20) of the quality scenario only; all other results are identical.
- **Data formats:** results are saved as CSV tables and NumPy `.npz` archives (open formats) instead of MATLAB `.mat` files. The values shown in the figures and in Table 2 are also saved as CSV files.
- **Table 2:** both quartiles and deciles of K are reported. In the submitted manuscript, the columns labelled D1 and D9 were the 25th and 75th percentiles.

___

# LICENSES

- **Code** (all `.py` files): GNU General Public License v3.0 or later (GPL-3.0-or-later), see `LICENSE`.
- **Data** (all files in `./data` and `./figures`, and the data files on Zenodo): Creative Commons Attribution 4.0 International (CC BY 4.0), see `DATA_LICENSE.md` and https://creativecommons.org/licenses/by/4.0/. You may share and adapt the data for any purpose, provided you give appropriate credit (cite the article and the Zenodo archive).

___

# HOW TO CITE

Please cite the article and the archived code and data (see `CITATION.cff`):
Elise J.-P. & Jenouvrier S. (2026). *Individual Heterogeneity along the Slow-Fast Continuum* (v2.0.0) [code and data]. Zenodo. https://doi.org/10.5281/zenodo.23298484

___

# CONTACT INFORMATION

**Contact person (questions about the code and data):**
Jean-Pierre Elise, École Normale Supérieure de Lyon, jean-pierre.elise@ens-lyon.fr, ORCID https://orcid.org/0009-0000-3899-5021

**Principal investigator:** Stéphanie Jenouvrier, Woods Hole Oceanographic Institution (WHOI), ORCID https://orcid.org/0000-0003-3324-2383

Bug reports and questions can also be posted as issues on the code repository: https://github.com/jp-elise/individual_heterogeneity_python/issues
