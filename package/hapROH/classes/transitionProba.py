"""
Compute the transition probabilities from one locus to an other,
depending on their genetic distance:
    P( state_{i+1} | state_{i} )
"""

from typing import NewType

import numpy as np

TransitionProba = NewType("TransitionProba", np.ndarray)
"""
(5, nb_snp-1) array of transition probabilities P(state_{i+1} | state_i)
The 5 probabilities are as follows:
    | stay_out_ROH | enter_ROH | leave_ROH | stay_in_ROH | jump_ROH |
aka |     A_00     |   A_01    |   A_10    |    A_11     |   A_12   |
"""


def get_transi_proba(
    r_in: float, r_out: float, r_jump: float, n_ref: int, r_map: np.ndarray
) -> TransitionProba:
    """Analytical exponentiation of the transition rate matrix."""
    # analytic formula for the exponentiaton of the transition rate
    S = r_in + r_out
    P1 = np.array([r_out, r_in / n_ref, r_out, r_in / n_ref, r_in / n_ref]) / S
    P2 = np.array([r_in, -r_in / n_ref, -r_out, r_out / n_ref, r_out / n_ref]) / S
    P3 = np.array([0, 0, 0, (n_ref - 1), -1]) / n_ref

    return TransitionProba(
        P1[:, None]
        + P2[:, None] * np.exp(-S * r_map)
        + P3[:, None] * np.exp(-(r_out + r_jump) * r_map)
    )


def _get_transi_proba_slow(
    r_in: float, r_out: float, r_jump: float, n_ref: int, r_map: np.ndarray
) -> TransitionProba:
    """Numerical exponentiation of the transition rate matrix (slower, for reference)."""
    from scipy.linalg import expm

    # compute transition rate matrix (3 states: not ROH, given ROH, any of the n-1 other ROH)
    T = np.array(
        [
            [0, r_in / n_ref, (n_ref - 1) / n_ref * r_in],
            [r_out, 0, (n_ref - 1) / n_ref * r_jump],
            [r_out, r_jump / n_ref, 0],
        ]
    )
    for i in range(3):
        T[i, i] = np.sum(-T[i])

    T = T[:, :, None] * r_map  # shape(3, 3, len(r_comb))
    T_exp = expm(T.T).T  # expm works along the last two axis -> need to transpose
    return TransitionProba(
        np.array(
            [
                T_exp[0, 0],
                T_exp[0, 1],
                T_exp[1, 0],
                T_exp[1, 1],
                T_exp[1, 2] / (n_ref - 1),
            ]
        )
    )


if __name__ == "__main__":
    from time import time

    r_in, r_out, r_jump, n_ref = 0.3, 1.2, 0.05, 7
    r_map = np.linspace(1e-4, 50, 1_000_000)

    t1 = time()
    old = _get_transi_proba_slow(r_in, r_out, r_jump, n_ref, r_map)
    t1 = time() - t1

    t2 = time()
    new = get_transi_proba(r_in, r_out, r_jump, n_ref, r_map)
    t2 = time() - t2

    print(f"Numeric method: {t1:.3}, Analytic method: {t2:.3}")

    assert np.abs(new - old).max() < 1e-10, (
        "Both computation methods should return the same result"
    )
    assert np.all(new <= 1), "Found probability > 1"
