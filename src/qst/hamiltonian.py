"""Single excitation Hamiltonians for an N-site chain (h_bar = 1)"""
from __future__ import annotations
import numpy as np

def _chain(couplings: np.ndarray, onsite: np.ndarray | None = None) -> np.ndarray:
    # Tridiagonal real-symmetric H from N-1 couplings and optional on-site energies.
    H = np.diag(couplings, 1) + np.diag(couplings, -1)
    if onsite is not None:
        H = H + np.diag(onsite)
    return H.astype(float)

def uniform_chain(n: int = 64, coupling: float = 10.0, onsite: float = -13.6) -> np.ndarray:
    # The 'normal' chain that we had before: constant coupling and constant on-site energy
    return _chain(np.full(n-1, coupling), np.full(n, onsite))

def perfect_chain(n: int = 64, lam: float = 1.0) -> np.ndarray:
    # Christandl et al. couplings J_j = (lam/2) sqrt(j (n-j)), j = 1..n-1.
    j = np.arange(1, n)
    return _chain(0.5 * lam * np.sqrt(j * (n - j)))

def add_barrier(H: np.ndarray, site: int, height: float):
    H = H.copy()
    H[site, site] += height
    return H

def add_coupling_noise(H: np.ndarray, sigma: float, seed: int | None = None):
    # Multiply each neighbor coupling by (1 + sigma * N(0,1)).
    rng = np.random.default_rng(seed)
    H = H.copy()
    for j in range(H.shape[0] - 1):
        H[j, j+1] = H[j+1, j] = H[j, j+1] * (1 + sigma * rng.standard_normal())
    return H


