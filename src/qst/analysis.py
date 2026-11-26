from __future__ import annotations
import numpy as np
from scipy.special import jv

def probabilities(psi_t: np.ndarray) -> np.ndarray:
    return np.abs(psi_t) ** 2

def bessel_reference(n_sites: int, coupling: float, t: float) -> np.ndarray:
    n = np.arange(n_sites)
    return (jv(n, 2 * coupling * t) + jv(n+2, 2 * coupling * t)) ** 2

def peak_transfer(P: np.ndarray, times: np.ndarray, target: int = -1):
    k = int(np.argmax(P[:, target]))
    return float(P[k, target]), float(times[k])
