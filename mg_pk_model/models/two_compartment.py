import math
import numpy as np
from core.derived import (
    D_IV, D_ORAL, V1, T_INF, C_BASE,
    ALPHA, BETA, A_2C, B_2C,
    K10, K12, K21,
    F_CIT, F_OX, KA_CIT, KA_OX,
)

def cp_2c_iv_bolus(t: np.ndarray) -> np.ndarray:
    scale = D_IV / V1
    return C_BASE + scale * (
        A_2C * np.exp(-ALPHA * t) +
        B_2C * np.exp(-BETA  * t)
    )

def cp_2c_iv_infusion(t: np.ndarray) -> np.ndarray:
    R0    = D_IV / T_INF
    scale = R0 / V1

    # concentration accumulated in each exponential component at end of infusion
    qa = (A_2C / ALPHA) * (1.0 - math.exp(-ALPHA * T_INF))
    qb = (B_2C / BETA)  * (1.0 - math.exp(-BETA  * T_INF))

    c_during = C_BASE + scale * (
        (A_2C / ALPHA) * (1.0 - np.exp(-ALPHA * t)) +
        (B_2C / BETA)  * (1.0 - np.exp(-BETA  * t))
    )
    c_after = C_BASE + scale * (
        qa * np.exp(-ALPHA * (t - T_INF)) +
        qb * np.exp(-BETA  * (t - T_INF))
    )
    return np.where(t <= T_INF, c_during, c_after)


def cp_2c_oral(t: np.ndarray, F: float, ka: float) -> np.ndarray:
    _ka = ka
    if abs(_ka - ALPHA) < 1e-6:
        _ka = ALPHA + 1e-5
    if abs(_ka - BETA) < 1e-6:
        _ka = BETA + 1e-5

    P =  (K21 - ALPHA) / ((_ka - ALPHA) * (BETA  - ALPHA))
    Qc = (K21 - BETA)  / ((_ka - BETA)  * (ALPHA - BETA))
    R  = -(P + Qc) # ensures P + Q + R = 0

    scale = F * D_ORAL * _ka / V1
    return C_BASE + scale * (
        P  * np.exp(-ALPHA * t) +
        Qc * np.exp(-BETA  * t) +
        R  * np.exp(-_ka   * t)
    )

def cp_2c_oral_citrate(t: np.ndarray) -> np.ndarray:
    return cp_2c_oral(t, F_CIT, KA_CIT)

def cp_2c_oral_oxide(t: np.ndarray) -> np.ndarray:
    return cp_2c_oral(t, F_OX, KA_OX)