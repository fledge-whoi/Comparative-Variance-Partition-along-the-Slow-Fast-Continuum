"""Unit tests of the functions/ package against analytical results.
Run from the project folder:  python -m unittest discover -s tests -v"""
import itertools
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from functions import (combinations, popHeteromatfact, det_outputfactorial, lifetimeoutput,  # noqa: E402
                       reprooutputfactorial, varcalc, vardecomposition,
                       minusparamfinder_ConstantLambda)

THETA = np.array([[0.01, 0.05, 0.9, 105.4500], [0.2, 0.1, 0.6, 6.9000], [0.5, 0.5, 0.4, 1.7500],
                  [0.7, 0.75, 0.1, 1.3214], [0.8, 0.9, 0.075, 0.4333], [0.85, 0.95, 0.05, 0.2265]])
G_ARTICLE = [2.05, 2.20, 3.43, 6.70, 13.85, 25.19]  # Table 1a (rounded)


def baseline(th):
    U, F, A = popHeteromatfact(th[None, :], [1.0], 1, 2)
    N = np.linalg.inv(np.eye(2) - U)
    return U, F, A, N


class TestBaselineSpecies(unittest.TestCase):
    def test_generation_times_match_table1(self):
        for th, g in zip(THETA, G_ARTICLE):
            lam, ssd, GT = det_outputfactorial(*baseline(th)[:3])
            self.assertAlmostEqual(round(GT, 2), g, places=2)

    def test_baseline_species_are_stable(self):
        for th in THETA:
            lam, _, _ = det_outputfactorial(*baseline(th)[:3])
            self.assertAlmostEqual(lam, 1.0, places=3)

    def test_mean_longevity_analytical(self):
        # expected number of years spent as juvenile + as adult, from birth
        for sj, sa, g, f in THETA:
            stayJ = sj * (1 - g)
            expected = 1 / (1 - stayJ) + (sj * g / (1 - stayJ)) / (1 - sa)
            U, F, A, N = baseline(np.array([sj, sa, g, f]))
            m, v = lifetimeoutput(N, np.eye(2), np.array([True, False]))
            self.assertAlmostEqual(m[0], expected, places=10)
            self.assertGreaterEqual(v[0], 0)

    def test_mean_LRO_equals_R0(self):
        for sj, sa, g, f in THETA:
            R0 = f * sj * g / ((1 - sa) * (1 - sj * (1 - g)))
            U, F, A, N = baseline(np.array([sj, sa, g, f]))
            m, v = reprooutputfactorial(U, F, N, np.array([True, False]), np.ones(3), np.zeros(2),
                                        np.hstack([np.eye(2), np.zeros((2, 1))]))
            self.assertAlmostEqual(m[0], R0, places=10)
            self.assertGreater(v[0], 0)


class TestHeterogeneousPopulation(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        self.Px = combinations([1, 0], [1, 0], [1, 0], [1, 0])
        self.Py = 1 - self.Px
        th = THETA[2]
        self.thetaxy = np.tile(th, (16, 1)) * 1.3 ** self.Px * 0.6 ** self.Py
        R = rng.uniform(0.05, 0.95, 4)
        self.Rcomb = combinations(*[[r, 1 - r] for r in R]).prod(axis=1)

    def test_combinations_ordering(self):
        np.testing.assert_array_equal(combinations([1, 0], [1, 0]), [[1, 1], [1, 0], [0, 1], [0, 0]])

    def test_offspring_distribution_and_block_structure(self):
        U, F, A = popHeteromatfact(self.thetaxy, self.Rcomb, 16, 2)
        # newborns of every group are spread over the juvenile states according to Rcomb
        col = F[:, 1]  # offspring of adults of group 1
        np.testing.assert_allclose(col[0::2] / col.sum(), self.Rcomb, atol=1e-14)
        np.testing.assert_allclose(col[1::2], 0)
        # no transitions between groups
        for i in range(16):
            for j in range(16):
                if i != j:
                    self.assertTrue(np.all(U[2*i:2*i+2, 2*j:2*j+2] == 0))

    def test_factorial_decomposition_against_brute_force(self):
        rng = np.random.default_rng(1)
        m = rng.uniform(1, 10, 16)
        PI = rng.dirichlet(np.ones(16))
        comps = vardecomposition(m.reshape([2] * 4, order="F"), PI.reshape([2] * 4, order="F"),
                                 self.Px, self.Py, 16, np.array([1, 1, 1, 1]))

        def VS(S):  # variance of the conditional means of the factors in S
            if not S:
                return 0.0
            keys = [tuple(self.Px[g, list(S)]) for g in range(16)]
            mu, tot = PI @ m, 0.0
            for k in set(keys):
                idx = [g for g in range(16) if keys[g] == k]
                w = PI[idx].sum()
                tot += w * (PI[idx] @ m[idx] / w - mu) ** 2
            return tot
        brute = {}
        for r in range(1, 5):
            for S in itertools.combinations(range(4), r):
                brute[S] = VS(S) - sum(brute[T] for T in brute if set(T) < set(S))
        for g in range(15):
            S = tuple(np.nonzero(self.Px[g])[0])
            self.assertAlmostEqual(comps[g], brute[S], places=12)
        self.assertAlmostEqual(comps[15], 0.0)
        self.assertAlmostEqual(comps.sum(), varcalc(m, PI), places=12)

    def test_constant_lambda_criterion(self):
        # R0 = 1 for the mixture  <=>  lambda = 1
        from scipy.optimize import brentq
        Px = np.vstack([np.hstack([np.zeros((8, 1)), np.ones((8, 1)), np.zeros((8, 2))]), np.zeros((8, 4))])
        Py = np.vstack([np.hstack([np.ones((8, 1)), np.zeros((8, 3))]), np.zeros((8, 4))])
        thetax = np.tile(THETA[3], (16, 1)) * 1.2 ** Px
        f = lambda Y: minusparamfinder_ConstantLambda(thetax, Y, np.ones((16, 4)), Py, self.Rcomb)
        Y = brentq(f, 0, 1)
        U, F, A = popHeteromatfact(thetax * Y ** Py, self.Rcomb, 16, 2)
        lam, _, _ = det_outputfactorial(U, F, A)
        self.assertAlmostEqual(lam, 1.0, places=10)


if __name__ == "__main__":
    unittest.main()
