"""Figures and animations for QST: ``uv python -m qst.plotting`` for regenerating media"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import PowerNorm

from .analysis import probabilities
from .evolve import evolve, site_state
from .hamiltonian import perfect_chain, uniform_chain

def fidelity_plot(P, times, target=-1, ax=None, label=None, mark_pi=False):
    """Line plot of the probability at one site versus time."""
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(times, P[:, target], label=label)
    if mark_pi:
        ax.axvline(np.pi, color="gray", linestyle="--", linewidth=1, label=r"$t=\pi$")
    ax.set_xlabel("time t (1/energy)")
    ax.set_ylabel("probability at target site")
    ax.set_ylim(0, 1.05)
    if label or mark_pi:
        ax.legend()
    return ax

def heatmap(P, times, ax=None, title=None):
    """Image of P: x = site, y = time, color = probability."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    n = P.shape[1]
    im = ax.imshow(
        P,
        aspect="auto",
        origin="lower",
        extent=[-0.5, n - 0.5, times[0], times[-1]],
        norm=PowerNorm(gamma=0.5, vmin=0, vmax=1),
        cmap="viridis",
    )
    ax.figure.colorbar(im, ax=ax, label="probability")
    ax.set_xlabel("site")
    ax.set_ylabel("time t (1/energy)")
    if title:
        ax.set_title(title)
    return ax


def animate(P, times, out, target=-1, fps=20, step=1):
    """Save a GIF of the probability along the chain; needs no ffmpeg (uses Pillow)."""
    n = P.shape[1]
    sites = np.arange(n)
    fig, ax = plt.subplots(figsize=(6, 3), dpi=80)
    (line,) = ax.plot(sites, P[0], "o-", markersize=3)
    ax.axvline(sites[target], color="red", alpha=0.3, label="target site")
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_xlabel("site")
    ax.set_ylabel("probability")
    ax.legend(loc="upper center")

    def update(k):
        line.set_ydata(P[k])
        ax.set_title(f"t = {times[k]:.2f}")
        return (line,)

    anim = FuncAnimation(fig, update, frames=range(0, len(times), step), blit=False)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)


def make_figures(outdir="media"):
    """Regenerate every figure used in the README."""
    plt.switch_backend("Agg")  # no window; works on servers and CI
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    times = np.arange(0, 20, 0.02)
    psi0 = site_state(64)
    P = {
        "uniform": probabilities(evolve(uniform_chain(64), psi0, times)),
        "perfect": probabilities(evolve(perfect_chain(64), psi0, times)),
    }
    for name, data in P.items():
        fig, ax = plt.subplots(figsize=(6, 5))
        heatmap(data, times, ax=ax, title=f"{name} chain")
        fig.savefig(out / f"heatmap_{name}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 3.5))
    for name, data in P.items():
        fidelity_plot(data, times, ax=ax, label=name)
    ax.axvline(np.pi, color="gray", linestyle="--", linewidth=1)
    fig.savefig(out / "fidelity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    anim_times = np.arange(0, 6.4, 0.1)  # 64 frames, covers t = pi
    for name, H in (("uniform", uniform_chain(64)), ("perfect", perfect_chain(64))):
        data = probabilities(evolve(H, psi0, anim_times))
        animate(data, anim_times, out / f"{name}.gif", fps=15)


if __name__ == "__main__":
    make_figures()
