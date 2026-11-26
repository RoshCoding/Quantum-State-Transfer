from __future__ import annotations
import numpy as np

def evolve(H: np.ndarray, psi0: np.ndarray, times: np.ndarray) -> np.ndarray:
    # Exact evolution psi(t) = V^exp(-iEt) + V^H*psi0. Returns shape (len(times), N)
    w, V = np.linalg.eigh(H)
    c = V.conj().T @ psi0
    phases = np.exp(-1j * np.outer(times, w))
    return (phases * c) @ V.T

def site_state(n: int, site: int = 0) -> np.ndarray:
    psi = np.zeros(n, dtype=complex)
    psi[site] = 1.0
    return psi



