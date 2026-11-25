import numpy as np
from pathlib import Path
from qst.hamiltonian import uniform_chain, perfect_chain
from qst.evolve import evolve, site_state
from qst.analysis import probabilities

DATA = Path(__file__).parent / "data"
T = np.arange(500) * 0.1


def _check(H, fname):
    ref = np.loadtxt(DATA / fname)
    P = probabilities(evolve(H, site_state(64), T))
    assert np.abs(P - ref).max() < 1e-6

def test_matches_c_perfect():
    _check(perfect_chain(64), "c_perfect.txt")
    
def test_matches_c_uniform():
    _check(uniform_chain(64), "c_uniform.txt")


