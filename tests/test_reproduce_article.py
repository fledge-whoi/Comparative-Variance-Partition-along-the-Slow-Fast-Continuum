"""Checks that the simulation results in ./data reproduce Table 2 of the
article (skipped if ./data has not been generated yet with main.py).
Species 2 (G = 2.20) of the quality scenario uses the corrected values
(see README, "Changes from the MATLAB version 1.0.0").
Trade-off sa/F, species 3: the manuscript reports Q1 = 0.0268, most likely a
typo for 0.0263 (the other three values of that row match to 4 decimals)."""
import os
import unittest

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# Q1, median, Q3, max of K (%) per species, Table 2 (Q1/Q3 were labelled D1/D9)
EXPECTED = {
    ("Quality_Gconstant", "IHeteroLEresults"): [
        (0.001, 0.005, 0.018, 0.745), (0.060, 0.496, 2.141, 25.397), (1.354, 5.559, 14.089, 48.359),
        (0.737, 2.621, 6.828, 32.025), (0.229, 0.852, 2.088, 17.489), (0.219, 0.523, 0.996, 12.764)],
    ("Quality_Gconstant", "IHeteroLROresults"): [
        (0.005, 0.014, 0.049, 0.795), (0.165, 0.962, 2.771, 27.726), (1.429, 5.440, 14.044, 48.927),
        (0.469, 1.582, 5.096, 32.080), (0.179, 0.527, 1.541, 17.918), (0.163, 0.357, 0.847, 12.841)],
    ("SaSjTradeoff_LAMBDAconstant", "IHeteroLEresults"): [
        (0.000, 0.000, 0.000, 0.000), (0.001, 0.007, 0.029, 0.105), (0.015, 0.079, 0.253, 0.534),
        (0.069, 0.337, 1.031, 2.533), (0.034, 0.185, 0.455, 1.075), (0.017, 0.096, 0.266, 0.626)],
    ("SaFTradeoff_LAMBDAconstant", "IHeteroLEresults"): [
        (0.000, 0.000, 0.000, 0.000), (0.000, 0.000, 0.005, 0.013), (0.0263, 0.1342, 0.4795, 1.4286),
        (0.0026, 0.0124, 0.0467, 0.1599), (0.0018, 0.0083, 0.0341, 0.1254), (0.0010, 0.0042, 0.0186, 0.0718)],
}


@unittest.skipUnless(os.path.exists(os.path.join(DATA, "simulation_results_Quality_Gconstant.npz")),
                     "run 'python main.py' first")
class TestTable2(unittest.TestCase):
    def test_table2(self):
        for (name, var), rows in EXPECTED.items():
            K = np.load(os.path.join(DATA, f"simulation_results_{name}.npz"))[var]
            for j, exp in enumerate(rows):
                k = 100 * K[:, :, j].ravel()
                k = k[~np.isnan(k)]
                got = [*np.quantile(k, [0.25, 0.5, 0.75], method="hazen"), k.max()]
                tol = 6e-5 if exp[0] != round(exp[0], 3) else 6e-4  # rounding of the published table
                for g, e in zip(got, exp):
                    with self.subTest(file=name, var=var, species=j + 1):
                        self.assertAlmostEqual(g, e, delta=tol)


if __name__ == "__main__":
    unittest.main()
