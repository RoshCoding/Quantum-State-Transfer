"""CLI: uv run python -m qst"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .analysis import peak_transfer, probabilities
from .evolve import evolve, site_state
from .hamiltonian import add_barrier, add_coupling_noise, perfect_chain, uniform_chain

def build_parser() -> argparse.ArgumentParser:
     p = argparse.ArgumentParser(prog="qst", description="Quantum state transfer on a chain")
     p.add_argument("--chain", choices=["uniform", "perfect"], default="perfect")
     p.add_argument("-n", type=int, default=64, help="number of sites")
     p.add_argument("--dt", type=float, default=0.1, help="time step")
     p.add_argument("--steps", type=int, default=500, help="number of time points")
     p.add_argument("--start", type=int, default=0, help="site where the particle starts")
     p.add_argument("--target", type=int, default=-1, help="site to monitor (-1 for last)")
     p.add_argument("--noise", type=float, default=0.0, help="relative coupling noise sigma")
     p.add_argument("--seed", type=int, default=None, help="random seed for noise")
     p.add_argument("--barrier-site", type=int, default=None)
     p.add_argument("--barrier-height", type=float, default=0.0)
     p.add_argument("-o", type=Path, default=None, help="Write probabilities to this file")
     return p

def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    H = perfect_chain(args.n) if args.chain == "perfect" else uniform_chain(args.n)
    if args.barrier_site is not None:
        H = add_barrier(H, args.barrier_site, args.barrier_height)
    if args.noise > 0:
        H = add_coupling_noise(H, args.noise, args.seed)

    times = np.arange(args.steps) * args.dt
    P = probabilities(evolve(H, site_state(args.n, args.start), times))

    if args.o is not None:
        args.o.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.o, P, fmt="%.12e")

    peak, t_peak = peak_transfer(P, times, args.target)
    print(f"chain={args.chain} n={args.n}")
    print(f"peak P(site {args.target}) = {peak:.6f} at t = {t_peak:.3f}")
