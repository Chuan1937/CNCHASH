"""Regenerate the README/docs figures from live measurements.

Run:  python benchmarks/make_figures.py
Outputs docs/images/speed_comparison.png and accuracy_verification.png.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cnchash import backend as backend_mod
from cnchash import run_hash
from cnchash.utils import fp_coord_angles_to_vectors

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests")
)
from test_accuracy import MECHANISMS, kagan_angle, make_synthetic_event  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IMG_DIR = os.path.join(HERE, "docs", "images")

# Performance figures (30 stations, nmc=30, dang=5 deg).
NATIVE = os.environ.get("CNCHASH_HASHHP_LIB")
ORIGINAL_MS = 145.0  # measured with hash_driver1 on the same machine


def make_event(nsta, nmc, seed=0):
    rng = np.random.default_rng(seed)
    az = rng.uniform(0, 360, nsta)
    the = rng.uniform(30, 95, nsta)
    pol = np.where((np.sin(np.deg2rad(the)) * np.cos(np.deg2rad(az - 60))) > 0, 1, -1).astype(
        np.int32
    )
    qual = (rng.random(nsta) < 0.2).astype(np.int32)
    return (
        np.repeat(az.reshape(-1, 1), nmc, axis=1),
        np.repeat(the.reshape(-1, 1), nmc, axis=1),
        pol,
        qual,
    )


def time_ms(events, num_threads, min_seconds=1.5):
    backend_mod.set_num_threads(num_threads)
    t0 = float("nan")
    import time as _time

    t0 = _time.perf_counter()
    n = 0
    while _time.perf_counter() - t0 < min_seconds:
        for azi, the, pol, qual in events:
            run_hash(azi, the, pol, qual, nmc=30, backend="fortran", num_threads=num_threads)
            n += 1
    return (_time.perf_counter() - t0) / n * 1000.0


def accuracy_data(rng, flip_frac):
    rots = []
    for mech in MECHANISMS:
        az, the, pol, qual = make_synthetic_event(mech, nsta=40, rng=rng, flip_frac=flip_frac)
        r = run_hash(az, the, pol, qual, nmc=1, backend="fortran", num_threads=1)
        s = np.atleast_1d(r["strike_avg"])[0]
        d = np.atleast_1d(r["dip_avg"])[0]
        rk = np.atleast_1d(r["rake_avg"])[0]
        n1, s1 = fp_coord_angles_to_vectors(*mech)
        n2, s2 = fp_coord_angles_to_vectors(s, d, rk)
        rots.append(kagan_angle(n1, s1, n2, s2))
    return np.array(rots)


def main():
    # Data measured on the development machine (30 stations, nmc=30,
    # dang=5). Re-measure with --measure if you want fresh numbers on
    # another machine (portable numbers require no CNCHASH_HASHHP_LIB;
    # native numbers require it to point at a -march=native build).
    portable = {1: 48.7, 4: 19.0}
    native = {1: 23.3, 16: 6.6}
    if "--measure" in sys.argv:
        events = [make_event(30, 30, seed=i) for i in range(10)]
        portable = {nt: time_ms(events, nt) for nt in (1, 4)}
        native = None
        if NATIVE:
            native = {nt: time_ms(events, nt) for nt in (1, 16)}

    # --- speed comparison ---
    labels = ["Portable\n1 thread", "Portable\n4 threads", "Native AVX2\n1 thread",
              "Native AVX2\n16 threads", "Original\nHASH v1.2"]
    values = [
        portable[1],
        portable[4],
        native[1] if native else None,
        native[16] if native else None,
        ORIGINAL_MS,
    ]
    values = [v for v in values if v is not None]
    labels = [lab for lab, _v in zip(labels, values, strict=True) if _v is not None]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(labels, values, color=["#4C72B0", "#55A868", "#DD8452", "#C44E52", "#8172B3"])
    for b, v in zip(bars, values, strict=True):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.02, f"{v:.1f}", ha="center", fontsize=11)
    ax.set_ylabel("ms / event (30 stations, nmc=30)")
    ax.set_title("CNCHASH vs original HASH v1.2 (same machine)")
    ax.set_ylim(0, max(values) * 1.15)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(IMG_DIR, "speed_comparison.png"), dpi=130)
    plt.close(fig)

    # --- accuracy ---
    rng_clean = np.random.default_rng(2026)
    rng_flip = np.random.default_rng(2026)
    rots_clean = accuracy_data(rng_clean, 0.0)
    rots_flip = accuracy_data(rng_flip, 0.1)

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(rots_clean))
    w = 0.38
    ax.bar(x - w / 2, rots_clean, w, label="no flipped polarities", color="#4C72B0")
    ax.bar(x + w / 2, rots_flip, w, label="~10% flipped polarities", color="#55A868")
    ax.axhline(25.0, color="red", ls="--", lw=1.2)
    ax.text(len(x) - 0.5, 26, "25 deg acceptance", color="red", ha="right", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}/{d}/{r}" for s, d, r in MECHANISMS], rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("Kagan rotation angle (deg)")
    ax.set_title("Synthetic recovery of known mechanisms")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(IMG_DIR, "accuracy_verification.png"), dpi=130)
    plt.close(fig)

    print("figures written to", IMG_DIR)
    print("portable:", {k: round(v, 1) for k, v in portable.items()},
          "native:", {k: round(v, 1) for k, v in native.items()} if native else None)
    print(f"clean median {np.median(rots_clean):.1f} deg, "
          f"flip median {np.median(rots_flip):.1f} deg")


if __name__ == "__main__":
    main()
