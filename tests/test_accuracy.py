"""Accuracy tests.

1. Synthetic recovery: events generated from known mechanisms (uniform
   full-focal-sphere station coverage, ~10% flipped polarities) must be
   recovered within a Kagan rotation angle of 25 degrees.
2. Original-HASH parity: when the HASH_complete Fortran sources can be
   compiled, the acceptable-mechanism sets produced by CNCHASH must
   match the original FOCALMC exactly on identical inputs.
"""

import os
import shutil
import subprocess

import numpy as np
import pytest

from cnchash import backend as backend_mod
from cnchash import run_hash
from cnchash.utils import fp_coord_angles_to_vectors

SEED = 2026
RNG = np.random.default_rng(SEED)

MECHANISMS = [
    (45, 60, 90),
    (120, 40, -60),
    (200, 75, 150),
    (300, 25, -120),
    (0, 90, 0),
    (75, 50, 30),
    (10, 15, -40),
    (250, 80, 20),
    (160, 65, 10),
    (350, 45, -90),
]

KAGAN_TOL = 25.0  # degrees


def predicted_polarity(az, the, strike, dip, rake):
    """Polarity predicted from the moment tensor of a mechanism."""
    sr, dr, rr = np.deg2rad([strike, dip, rake])
    M11 = -np.sin(dr) * np.cos(rr) * np.sin(2 * sr) - np.sin(2 * dr) * np.sin(rr) * np.sin(sr) ** 2
    M22 = np.sin(dr) * np.cos(rr) * np.sin(2 * sr) - np.sin(2 * dr) * np.sin(rr) * np.cos(sr) ** 2
    M33 = np.sin(2 * dr) * np.sin(rr)
    M12 = np.sin(dr) * np.cos(rr) * np.cos(2 * sr) + 0.5 * np.sin(2 * dr) * np.sin(rr) * np.sin(2 * sr)
    M13 = -np.cos(dr) * np.cos(rr) * np.cos(sr) - np.cos(2 * dr) * np.sin(rr) * np.sin(sr)
    M23 = -np.cos(dr) * np.cos(rr) * np.sin(sr) + np.cos(2 * dr) * np.sin(rr) * np.cos(sr)
    M = np.array([[M11, M12, M13], [M12, M22, M23], [M13, M23, M33]])
    out = np.zeros(len(az), dtype=np.int32)
    for k in range(len(az)):
        tr, pr = np.deg2rad(the[k]), np.deg2rad(az[k])
        a = np.array([np.sin(tr) * np.cos(pr), np.sin(tr) * np.sin(pr), -np.cos(tr)])
        out[k] = -1 if (a @ M @ a) < 0 else 1
    return out


def uniform_sphere_stations(nsta, rng):
    """Uniformly distributed takeoff/azimuth directions (10..170 deg takeoff)."""
    az = np.zeros(nsta)
    the = np.zeros(nsta)
    got = 0
    while got < nsta:
        v = rng.normal(size=(nsta * 3, 3))
        v /= np.linalg.norm(v, axis=1, keepdims=True)
        z = v[:, 2]
        sel = np.where(np.abs(z) < 0.985)[0][: nsta - got]
        if len(sel) == 0:
            continue
        the[got : got + len(sel)] = np.degrees(np.arccos(np.clip(z[sel], -1, 1)))
        az[got : got + len(sel)] = np.degrees(np.arctan2(v[sel, 1], v[sel, 0])) % 360.0
        got += len(sel)
    return az, the


def make_synthetic_event(mechanism, nsta=40, rng=None, flip_frac=0.1):
    """Full-focal-sphere synthetic event with ~10% flipped polarities."""
    if rng is None:
        rng = RNG
    az, the = uniform_sphere_stations(nsta, rng)
    pol = predicted_polarity(az, the, *mechanism)
    if flip_frac > 0:
        nbad = max(1, int(len(pol) * flip_frac))
        bad = rng.choice(len(pol), nbad, replace=False)
        pol[bad] = -pol[bad]
    qual = np.zeros(nsta, dtype=np.int32)
    return az, the, pol, qual


def kagan_angle(n1, s1, n2, s2):
    """Minimum rotation angle (degrees) between two mechanisms."""
    best = 180.0
    B1 = np.cross(n1, s1)
    for swap in (False, True):
        for neg in (False, True):
            a, b = (n2, s2) if not swap else (s2, n2)
            if neg:
                a, b = -a, -b
            B2 = np.cross(a, b)
            phi = [
                np.arccos(np.clip(np.dot(n1, a), -1, 1)),
                np.arccos(np.clip(np.dot(s1, b), -1, 1)),
                np.arccos(np.clip(np.dot(B1, B2), -1, 1)),
            ]
            if phi[0] < 1e-4 and phi[1] < 1e-4 and phi[2] < 1e-4:
                rot = 0.0
            elif phi[0] < 1e-4:
                rot = np.degrees(phi[1])
            elif phi[1] < 1e-4:
                rot = np.degrees(phi[2])
            elif phi[2] < 1e-4:
                rot = np.degrees(phi[0])
            else:
                dn = np.stack([n1 - a, s1 - b, B1 - B2])
                dn /= np.linalg.norm(dn, axis=1, keepdims=True)
                q = [np.dot(dn[0], dn[1]), np.dot(dn[0], dn[2]), np.dot(dn[1], dn[2])]
                iout = next((i for i, v in enumerate(q) if v > 0.9999), None)
                if iout is None:
                    iout = int(np.argmin(np.linalg.norm(dn, axis=1)))
                keep = [i for i in range(3) if i != iout]
                R = np.cross(dn[keep[0]], dn[keep[1]])
                R /= np.linalg.norm(R)
                theta = [
                    np.arccos(np.clip(np.dot(n1, R), -1, 1)),
                    np.arccos(np.clip(np.dot(s1, R), -1, 1)),
                    np.arccos(np.clip(np.dot(B1, R), -1, 1)),
                ]
                iuse = int(np.argmin([abs(t - np.pi / 2) for t in theta]))
                c = np.cos(phi[iuse]) - np.cos(theta[iuse]) ** 2
                s = np.sin(theta[iuse]) ** 2
                rot = np.degrees(np.arccos(np.clip(c / s, -1, 1))) if abs(s) > 1e-10 else 0.0
            best = min(best, rot)
    return best


@pytest.fixture(scope="module")
def fortran_available():
    if not backend_mod.get_backend("fortran").available:
        pytest.skip("fortran backend not available")
    return True


def test_synthetic_recovery_success_rate(fortran_available):
    """10 mechanisms x 40 stations, ~10% flips: >= 8/10 recovered within 25 deg."""
    rng = np.random.default_rng(SEED)
    rots = []
    for mech in MECHANISMS:
        az, the, pol, qual = make_synthetic_event(mech, nsta=40, rng=rng)
        r = run_hash(az, the, pol, qual, nmc=1, backend="fortran", num_threads=1)
        assert r["success"], f"no solution for mechanism {mech}"
        s = np.atleast_1d(r["strike_avg"])[0]
        d = np.atleast_1d(r["dip_avg"])[0]
        rk = np.atleast_1d(r["rake_avg"])[0]
        n1, s1 = fp_coord_angles_to_vectors(*mech)
        n2, s2 = fp_coord_angles_to_vectors(s, d, rk)
        rots.append(kagan_angle(n1, s1, n2, s2))
    rots = np.array(rots)
    assert (rots <= KAGAN_TOL).sum() >= 8, (
        f"only {(rots <= KAGAN_TOL).sum()}/10 recovered; rotations={np.round(rots, 1)}"
    )
    assert np.median(rots) < 20.0


def test_synthetic_recovery_no_flips(fortran_available):
    """Without flipped polarities every mechanism is recovered."""
    rng = np.random.default_rng(SEED + 1)
    for mech in MECHANISMS[:6]:
        az, the, pol, qual = make_synthetic_event(mech, nsta=40, rng=rng, flip_frac=0.0)
        r = run_hash(az, the, pol, qual, nmc=1, backend="fortran", num_threads=1)
        assert r["success"]
        s = np.atleast_1d(r["strike_avg"])[0]
        d = np.atleast_1d(r["dip_avg"])[0]
        rk = np.atleast_1d(r["rake_avg"])[0]
        n1, s1 = fp_coord_angles_to_vectors(*mech)
        n2, s2 = fp_coord_angles_to_vectors(s, d, rk)
        assert kagan_angle(n1, s1, n2, s2) < 15.0


# ---------------------------------------------------------------------------
# Original HASH parity (requires a Fortran compiler and the HASH_complete
# sources; skipped when they are not available).
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def original_hash_bin(tmp_path_factory):
    """Compile the original HASH_complete core harness once per session."""
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    hash_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "HASH_complete")
    src_dir = os.path.join(hash_dir, "src")
    obj_dir = os.path.join(hash_dir, "obj")
    if not os.path.isdir(src_dir):
        pytest.skip("HASH_complete sources not present")
    os.makedirs(obj_dir, exist_ok=True)

    # Compile the original subroutines if not already built.
    subs = ["fmech_subs", "pol_subs", "station_subs", "uncert_subs", "util_subs"]
    for sub in subs:
        obj = os.path.join(obj_dir, f"{sub}.o")
        if not os.path.exists(obj):
            r = subprocess.run(
                ["gfortran", "-O", "-I", os.path.join(src_dir, "include"), "-c",
                 os.path.join(src_dir, "subs", f"{sub}.f"), "-o", obj],
                capture_output=True, text=True,
            )
            if r.returncode != 0:
                pytest.skip(f"could not compile HASH_complete/{sub}.f")

    harness = str(tmp_path_factory.mktemp("hashcmp") / "fcompare")
    if not os.path.exists(harness):
        src = """
      program fcompare
      parameter (npick0=500,nmc0=500,nmax0=500,ncoor=31032)
      integer npol, nmc, maxout, nextra, nmismax, nf, nsltn, nout2
      real p_azi_mc(npick0,nmc0), p_the_mc(npick0,nmc0)
      integer p_pol(npick0), p_qual(npick0)
      real strike(nmax0), dip(nmax0), rake(nmax0)
      real faults(3,nmax0), slips(3,nmax0)
      real str_avg(5), dip_avg(5), rak_avg(5)
      real prob(5), var_est(2,5)
      real mfrac, stdr, magap, mpgap
      real badfrac, cangle, prob_max, dang
      real azi, the
      integer pol, qual, i
      read(*,*) npol
      read(*,*) dang
      read(*,*) badfrac
      read(*,*) cangle
      read(*,*) prob_max
      read(*,*) nmc
      read(*,*) maxout
      do i=1,npol
        read(*,*) azi, the, pol, qual
        p_azi_mc(i,1)=azi
        p_the_mc(i,1)=the
        p_pol(i)=pol
        p_qual(i)=qual
      end do
      nmismax=max(nint(npol*badfrac),2)
      nextra=max(nint(npol*badfrac*0.5),2)
      call GET_GAP(npol,p_azi_mc,p_the_mc,magap,mpgap)
      call FOCALMC(p_azi_mc,p_the_mc,p_pol,p_qual,npol,nmc,
     &       dang,maxout,nextra,nmismax,nf,strike,dip,
     &       rake,faults,slips)
      nout2=min(maxout,nf)
      call MECH_PROB(nout2,faults,slips,cangle,prob_max,nsltn,
     &          str_avg,dip_avg,rak_avg,prob,var_est)
      call GET_MISF(npol,p_azi_mc,p_the_mc,p_pol,p_qual,str_avg(1),
     &       dip_avg(1),rak_avg(1),mfrac,stdr)
      write(*,'(A,I6)') 'nf=', nf
      write(*,'(A,I6)') 'nsltn=', nsltn
      write(*,'(A,6F12.6)') 'soln=', str_avg(1),dip_avg(1),rak_avg(1),
     &       prob(1),var_est(1,1),var_est(2,1)
      write(*,'(A,2F12.6)') 'misf=', mfrac, stdr
      do i=1,nf
        write(*,'(A,I5,3F12.6)') 'm ', i, strike(i),dip(i),rake(i)
      end do
      end
"""
        hsrc = os.path.join(os.path.dirname(harness), "fcompare.f")
        with open(hsrc, "w") as f:
            f.write(src)
        objs = " ".join(os.path.join(obj_dir, f"{s}.o") for s in subs)
        r = subprocess.run(
            ["gfortran", "-O", "-I", os.path.join(src_dir, "include"), hsrc,
             *objs.split(), "-o", harness],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            pytest.skip("could not build original HASH harness")
    return harness


def _run_original(bin_path, az, the, pol, qual):
    n = len(pol)
    header = "\n".join([str(n), "5.0", "0.1", "45.0", "0.1", "1", "500"])
    data = header + "\n" + "\n".join(
        [f"{a:.8f} {t:.8f} {p} {q}" for a, t, p, q in zip(az, the, pol, qual, strict=True)]
    )
    r = subprocess.run([bin_path], input=data, capture_output=True, text=True, timeout=60)
    out = {}
    for line in r.stdout.splitlines():
        line = line.strip()
        if line.startswith("nf="):
            out["nf"] = int(line.split("=")[1])
        elif line.startswith("nsltn="):
            out["nsltn"] = int(line.split("=")[1])
        elif line.startswith("soln="):
            out["soln"] = [float(x) for x in line.split("=")[1].split()]
        elif line.startswith("misf="):
            out["misf"] = tuple(float(x) for x in line.split("=")[1].split())
        elif line.startswith("m "):
            parts = line.split()
            out.setdefault("mech", []).append((float(parts[2]), float(parts[3]), float(parts[4])))
    return out


def test_original_hash_mechanism_sets_match(original_hash_bin, fortran_available):
    """Acceptable mechanism sets must match the original HASH exactly."""
    rng = np.random.default_rng(SEED + 2)
    for mech in MECHANISMS[:5]:
        az, the, pol, qual = make_synthetic_event(mech, nsta=40, rng=rng)
        fo = _run_original(original_hash_bin, az, the, pol, qual)
        r = run_hash(az, the, pol, qual, nmc=1, backend="fortran", num_threads=1)
        assert r["nf"] == fo["nf"], f"{mech}: nf {r['nf']} vs original {fo['nf']}"
        so = np.array(fo["mech"])
        sc = np.column_stack([r["strike"], r["dip"], r["rake"]])
        # rake is periodic in 360 degrees (180 == -180)
        so[:, 2] = np.mod(so[:, 2], 360.0)
        sc[:, 2] = np.mod(sc[:, 2], 360.0)
        so = np.sort(so, axis=0)
        sc = np.sort(sc, axis=0)
        assert np.allclose(so[:, 0], sc[:, 0], atol=1e-6)
        assert np.allclose(so[:, 1], sc[:, 1], atol=1e-6)
        assert np.allclose(so[:, 2], sc[:, 2], atol=1e-6)


def test_original_hash_preferred_mechanism_close(original_hash_bin, fortran_available):
    """Preferred mechanisms must agree within nodal-plane ambiguity."""
    rng = np.random.default_rng(SEED + 3)
    for mech in MECHANISMS[:5]:
        az, the, pol, qual = make_synthetic_event(mech, nsta=40, rng=rng)
        fo = _run_original(original_hash_bin, az, the, pol, qual)
        r = run_hash(az, the, pol, qual, nmc=1, backend="fortran", num_threads=1)
        sc = np.atleast_1d(r["strike_avg"])[0]
        dc = np.atleast_1d(r["dip_avg"])[0]
        rc = np.atleast_1d(r["rake_avg"])[0]
        so, do_, ro = fo["soln"][:3]
        d_strike = min(abs(sc - so) % 360.0, abs(sc - (so + 180.0) % 360.0) % 360.0)
        assert d_strike < 5.0, f"{mech}: strike {sc:.1f} vs original {so:.1f}"
        assert abs(dc - do_) < 5.0, f"{mech}: dip {dc:.1f} vs original {do_:.1f}"
        assert min(abs(rc - ro) % 360.0, 360.0 - abs(rc - ro) % 360.0) < 10.0
