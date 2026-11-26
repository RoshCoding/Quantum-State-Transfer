import matplotlib

matplotlib.use("Agg") # Gives us headless with no window

import numpy as np

from qst.analysis import probabilities
from qst.evolve import evolve, site_state
from qst.hamiltonian import perfect_chain
from qst.plotting import animate, heatmap

def _data():
    times = np.arange(0, 4, 0.1)
    return probabilities(evolve(perfect_chain(16), site_state(16), times)), times

def test_heatmap_runs():
    P, t = _data()
    ax = heatmap(P, t)
    assert ax.get_xlabel() == "site"

def test_animate_writes_gif(tmp_path):
    P, t = _data()
    animate(P, t, tmp_path / "a.gif", step=4)
    assert (tmp_path / "a.gif").stat().st_size > 0
