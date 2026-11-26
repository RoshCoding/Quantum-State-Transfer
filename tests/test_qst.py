"""This file is for testing to ensure that the physics is working properly."""
import numpy as np
from qst.hamiltonian import uniform_chain, perfect_chain
from qst.evolve import evolve, site_state
from qst.analysis import probabilities, bessel_reference, peak_transfer

def test_norm_conserved():
    H = uniform_chain(64)
    t = np.linspace(0, 50, 500, endpoint=False)
    P = probabilities(evolve(H, site_state(64), t))
    assert np.allclose(P.sum(axis=1), 1.0, atol=1e-10)

def test_perfect_transfer_at_pi():
    P = probabilities(evolve(perfect_chain(64), site_state(64), np.array([np.pi])))
    assert P[0, -1] > 1 - 1e-9 and P[0, 0] < 1e-9

def test_perfect_revival_at_2pi():
    P = probabilities(evolve(perfect_chain(64), site_state(64), np.array([2*np.pi])))
    assert P[0, 0] > 1 - 1e-9

def test_uniform_matches_bessel_before_reflection():
    J = 10.0
    H = uniform_chain(64, coupling=J)
    for t in (0.1, 0.5, 1.0):
        P = probabilities(evolve(H, site_state(64), np.array([t])))[0]
        assert np.abs(P - bessel_reference(64, J, t)).max() < 1e-8

def test_uniform_never_transfers_well():
    t = np.arange(0, 50, 0.1)
    P = probabilities(evolve(uniform_chain(64), site_state(64), t))
    assert peak_transfer(P, t)[0] < 0.5
