"""Backend benchmarks for CNCHASH.

Measures events/second for the Numba and Fortran backends across event
sizes, MC trial counts, grid spacings, and thread counts. Results are
printed as a table; the numbers are measurements, not claims.

Usage:
    python benchmarks/benchmark_backends.py [--events N] [--seconds S]
"""

import argparse
import time

import numpy as np

from cnchash import backend as backend_mod
from cnchash import run_hash, run_hash_batch


def make_event(nsta, nmc, seed=0):
    rng = np.random.default_rng(seed)
    az = rng.uniform(0, 360, nsta)
    the = rng.uniform(30, 95, nsta)
    pol = np.where((np.sin(np.deg2rad(the)) * np.cos(np.deg2rad(az - 60))) > 0, 1, -1).astype(np.int32)
    qual = (rng.random(nsta) < 0.2).astype(np.int32)
    azi = np.repeat(az.reshape(-1, 1), nmc, axis=1)
    the_mc = np.repeat(the.reshape(-1, 1), nmc, axis=1)
    return azi, the_mc, pol, qual


def time_run_hash(backend, events, num_threads, min_seconds=1.0):
    """Time run_hash over many events; returns events/second."""
    backend_mod.set_num_threads(num_threads)
    t0 = time.perf_counter()
    n = 0
    while time.perf_counter() - t0 < min_seconds:
        for azi, the, pol, qual in events:
            run_hash(azi, the, pol, qual, nmc=azi.shape[1], backend=backend, num_threads=num_threads)
            n += 1
    elapsed = time.perf_counter() - t0
    return n / elapsed


def time_run_batch(backend, events, num_threads, min_seconds=1.0):
    backend_mod.set_num_threads(num_threads)
    t0 = time.perf_counter()
    n = 0
    while time.perf_counter() - t0 < min_seconds:
        run_hash_batch(events, nmc=events[0][0].shape[1], backend=backend, num_threads=num_threads)
        n += len(events)
    elapsed = time.perf_counter() - t0
    return n / elapsed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", type=int, default=1000)
    parser.add_argument("--seconds", type=float, default=1.5)
    args = parser.parse_args()

    configs = [
        ("small", 12, 30, 5.0),
        ("medium", 30, 30, 5.0),
        ("large", 80, 30, 5.0),
        ("dense", 30, 30, 2.0),
        ("high_mc", 30, 100, 5.0),
    ]

    backends = backend_mod.available_backends()
    threads = [1, 4]

    print(f"{'config':<10} {'backend':<8} {'threads':<8} {'events/s':>12} {'ms/event':>10}")
    print("-" * 55)
    for name, nsta, nmc, _dang in configs:
        events = [make_event(nsta, nmc, seed=i) for i in range(10)]
        for backend in backends:
            for nt in threads:
                eps = time_run_hash(backend, events, nt, args.seconds)
                ms = 1000.0 / eps
                print(f"{name:<10} {backend:<8} {nt:<8} {eps:>12.1f} {ms:>9.3f}")

    print("\nBatch mode (event-level vs scalar; medium config):")
    events = [make_event(30, 30, seed=i) for i in range(args.events)]
    for backend in backends:
        for nt in threads:
            eps_b = time_run_batch(backend, events, nt, args.seconds)
            eps_s = time_run_hash(backend, events, nt, args.seconds)
            print(f"  {backend:<8} threads={nt}: batch {eps_b:>10.1f} ev/s, scalar {eps_s:>10.1f} ev/s, batch speedup {eps_b / max(eps_s, 1e-9):.2f}x")


if __name__ == "__main__":
    main()
