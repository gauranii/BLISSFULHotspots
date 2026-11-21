# -----------------------------------------------------------------------------
# File:    python/tests/inputs.py
# Paper:   Gauran, Park, Rattsev, Peterson, Kann, Park (2022). Ann. Appl. Stat.
#          16(3):1459-1475.
# Purpose: Small inputs shared by the tests: a simulated data set, the state of
#          a chain and short sampler settings.
# Usage:   python -m pytest python/tests/inputs.py
# -----------------------------------------------------------------------------
import numpy as np

SETTING = {"run": "run0219", "setting": 17, "N": 100, "eta": 0.4, "lambda": 2.0, "theta": 0.3,
           "pi0": 0.8, "f1": "binomial", "C": 5}

# Short chains for the tests (the paper: 5000 and 1000, 10000 and 1000).
SHORT = {"alpha": 0.05, "replications": 1000, "sim_draws": 150, "sim_burn": 60,
         "real_draws": 120, "real_burn": 40, "real_seed": 1234,
         "real_gamma": {"NP": 1.5, "SP": 0.02, "P": 0.5, "GP": 0.5},
         "figS2_seed": 1439, "figure_size": 7}

# Counts: a null part near zero and a well-separated part.
COUNTS = np.array([0] * 30 + [1] * 14 + [2] * 8 + [3] * 4 + [4] * 2 + [5]
                  + [14, 15, 15, 16, 16, 16, 17, 17, 18, 20], dtype=float)


def state(U, **changes):
    """State of a chain at its start for U distinct counts."""
    s = {"eta": 0.5, "lambda": 2.0, "theta": 0.25, "pi0": 0.5, "xi": 1.0, "delta": 1.0, "nu": 0.5,
         "tau": 1.0, "psi_sp": np.full(U, 1 / 50), "psi_np": np.full(U, 1 / 50)}
    s.update(changes)
    return s
