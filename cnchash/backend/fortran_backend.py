"""Fortran native backend.

Loads ``libhashhp`` (the Modern Fortran/OpenMP core) through a stable
ISO_C_BINDING ABI and dispatches the full HASH pipeline to it. Falls
back gracefully (``available=False``) when the library is not found.

Library discovery order:
1. ``CNCHASH_HASHHP_LIB`` environment variable (explicit path)
2. ``<package>/lib/libhashhp.so`` (installed layout)
3. ``<repo>/src/hash_hp/libhashhp.so`` (in-tree build)
4. ``ctypes.util.find_library("hashhp")``
"""

import ctypes
import ctypes.util
import os
import threading

import numpy as np

from .base import HashBackend

_QUALITY_LETTERS = "ABCDEF"


class _CEventResult(ctypes.Structure):
    """Mirror of hash_c_api.c_event_result_t (bind(C) layout)."""

    _fields_ = [
        ("success", ctypes.c_int32),
        ("nf", ctypes.c_int32),
        ("nout1", ctypes.c_int32),
        ("nout2", ctypes.c_int32),
        ("nsltn", ctypes.c_int32),
        ("quality", ctypes.c_int32),
        ("strike_avg", ctypes.c_double * 5),
        ("dip_avg", ctypes.c_double * 5),
        ("rake_avg", ctypes.c_double * 5),
        ("prob", ctypes.c_double * 5),
        ("var_est", (ctypes.c_double * 2) * 5),
        ("mfrac", ctypes.c_double),
        ("mavg", ctypes.c_double),
        ("stdr", ctypes.c_double),
    ]


def _synchronized(method):
    """Serialize access to the process-global native state."""

    def wrapper(self, *args, **kwargs):
        with self._lock:
            return method(self, *args, **kwargs)

    return wrapper


def _as_c_double_array(arr):
    return np.ascontiguousarray(arr, dtype=np.float64)


def _as_c_int32_array(arr):
    return np.ascontiguousarray(arr, dtype=np.int32)


def _find_library():
    """Locate libhashhp using the documented search order."""
    env = os.environ.get("CNCHASH_HASHHP_LIB")
    if env and os.path.exists(env):
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    for rel in (
        "../lib/libhashhp.so",
        "../../src/hash_hp/libhashhp.so",
        "../../build/libhashhp.so",
    ):
        candidate = os.path.normpath(os.path.join(here, rel))
        if os.path.exists(candidate):
            return candidate
    found = ctypes.util.find_library("hashhp")
    if found:
        return found
    return None


class FortranBackend(HashBackend):
    """ctypes binding to the Modern Fortran/OpenMP core."""

    name = "fortran"

    def __init__(self):
        self._lib_path = _find_library()
        self._lib = None
        self.available = False
        self._tables = {}
        self._table_ip = 0
        # The native core keeps process-global state (rotation grid cache,
        # thread-private workspaces), so independent Python-thread calls
        # must be serialized at this boundary. OpenMP threading inside a
        # single call is unaffected.
        self._lock = threading.RLock()
        if self._lib_path is not None:
            try:
                self._lib = ctypes.CDLL(self._lib_path)
                self._bind_functions()
                self.available = True
            except (OSError, AttributeError):
                self._lib = None
                self._lib_path = None

    # ---- ctypes plumbing ----------------------------------------------------

    def _bind_functions(self):
        lib = self._lib
        lib.cnchash_version.argtypes = [
            ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32),
        ]
        lib.cnchash_version.restype = None

        lib.cnchash_set_num_threads.argtypes = [ctypes.c_int32]
        lib.cnchash_set_num_threads.restype = None
        lib.cnchash_get_max_threads.argtypes = []
        lib.cnchash_get_max_threads.restype = ctypes.c_int32

        lib.cnchash_get_rotation_grid.argtypes = [
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int32),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.c_int32,
        ]
        lib.cnchash_get_rotation_grid.restype = None

        lib.cnchash_get_gap.argtypes = [
            ctypes.c_int32,
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.cnchash_get_gap.restype = None

        lib.cnchash_get_misfit.argtypes = [
            ctypes.c_int32,
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.cnchash_get_misfit.restype = None

        lib.cnchash_get_misfit_amp.argtypes = [
            ctypes.c_int32,
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.cnchash_get_misfit_amp.restype = None

        lib.cnchash_mech_prob.argtypes = [
            ctypes.c_int32,
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int32),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
        ]
        lib.cnchash_mech_prob.restype = None

        lib.cnchash_mech_rot.argtypes = [
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.cnchash_mech_rot.restype = None

        lib.cnchash_run_event.argtypes = [
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.POINTER(_CEventResult),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
        ]
        lib.cnchash_run_event.restype = None

        lib.cnchash_run_event_amp.argtypes = [
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.POINTER(_CEventResult),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
        ]
        lib.cnchash_run_event_amp.restype = None

        lib.cnchash_run_batch.argtypes = [
            ctypes.c_int32,
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.POINTER(_CEventResult),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
        ]
        lib.cnchash_run_batch.restype = None

        lib.cnchash_run_batch_amp.argtypes = [
            ctypes.c_int32,
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.int32, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.POINTER(_CEventResult),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
        ]
        lib.cnchash_run_batch_amp.restype = None

        lib.cnchash_build_velocity_table.argtypes = [
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int32,
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.POINTER(ctypes.c_int32),
            ctypes.POINTER(ctypes.c_int32),
        ]
        lib.cnchash_build_velocity_table.restype = None

        lib.cnchash_get_tts.argtypes = [
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ctypes.c_int32,
            ctypes.c_int32,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int32),
        ]
        lib.cnchash_get_tts.restype = None

    # ---- interface ----------------------------------------------------------

    def info(self):
        return {
            "backend": "fortran",
            "openmp": True,
            "threads": self.get_num_threads(),
            "compiler": "Modern Fortran (gfortran)",
            "precision": "float64",
            "library": self._lib_path,
            "native_version": self.native_version(),
        }

    @_synchronized
    def native_version(self):
        """Version of the native libhashhp core as 'major.minor.patch'."""
        major, minor, patch = ctypes.c_int32(), ctypes.c_int32(), ctypes.c_int32()
        self._lib.cnchash_version(ctypes.byref(major), ctypes.byref(minor), ctypes.byref(patch))
        return f"{major.value}.{minor.value}.{patch.value}"

    @_synchronized
    def set_num_threads(self, num_threads):
        self._lib.cnchash_set_num_threads(int(num_threads))
        return self.get_num_threads()

    @_synchronized
    def get_num_threads(self):
        return int(self._lib.cnchash_get_max_threads())

    # ---- per-event pipeline -------------------------------------------------

    def _unpack_result(self, cres, npol, nspr):
        """Convert a C result to the public result dict."""
        nsltn = int(cres.nsltn)
        quality = _QUALITY_LETTERS[int(cres.quality)]
        if not cres.success:
            return {
                "success": False,
                "strike_avg": 999.0,
                "dip_avg": 99.0,
                "rake_avg": 999.0,
                "var_est": [99.0, 99.0],
                "mfrac": float(cres.mfrac),
                "mavg": float(cres.mavg),
                "quality": quality,
                "prob": 0.0,
                "stdr": float(cres.stdr),
                "nmult": 0,
                "nout2": int(cres.nout2),
                "nout1": int(cres.nout1),
                "nf": int(cres.nf),
                "npol": npol,
                "nspr": nspr,
            }
        var_est = [float(cres.var_est[0][0]), float(cres.var_est[0][1])]
        output = {
            "success": True,
            "strike_avg": float(cres.strike_avg[0]) if nsltn > 0 else 999.0,
            "dip_avg": float(cres.dip_avg[0]) if nsltn > 0 else 99.0,
            "rake_avg": float(cres.rake_avg[0]) if nsltn > 0 else 999.0,
            "var_est": var_est,
            "mfrac": float(cres.mfrac),
            "mavg": float(cres.mavg),
            "quality": quality,
            "prob": float(cres.prob[0]) if nsltn > 0 else 0.0,
            "stdr": float(cres.stdr),
            "nmult": nsltn,
            "nout2": int(cres.nout2),
            "nout1": int(cres.nout1),
            "nf": int(cres.nf),
            "npol": npol,
            "nspr": nspr,
        }
        if nsltn > 1:
            output["strike_avg"] = [float(cres.strike_avg[i]) for i in range(nsltn)]
            output["dip_avg"] = [float(cres.dip_avg[i]) for i in range(nsltn)]
            output["rake_avg"] = [float(cres.rake_avg[i]) for i in range(nsltn)]
            output["prob"] = [float(cres.prob[i]) for i in range(nsltn)]
            output["rms_diff"] = np.array(
                [[cres.var_est[i][0], cres.var_est[i][1]] for i in range(nsltn)]
            ).T
            output["quality"] = [quality] * nsltn
        return output

    @_synchronized
    def run_event(
        self,
        p_azi_mc,
        p_the_mc,
        p_pol,
        p_qual,
        npol,
        nmc,
        dang,
        maxout,
        badfrac,
        cangle,
        prob_max,
        npolmin,
        max_agap,
        max_pgap,
        selection=0,
    ):
        azi = _as_c_double_array(p_azi_mc).ravel(order="F")
        the = _as_c_double_array(p_the_mc).ravel(order="F")
        pol = _as_c_int32_array(p_pol)
        qual = _as_c_int32_array(p_qual)
        res = _CEventResult()
        maxout = int(maxout)
        s_all = np.zeros(maxout, np.float64)
        d_all = np.zeros(maxout, np.float64)
        r_all = np.zeros(maxout, np.float64)
        f_all = np.zeros(3 * maxout, np.float64)
        sl_all = np.zeros(3 * maxout, np.float64)
        self._lib.cnchash_run_event(
            azi,
            the,
            pol,
            qual,
            int(npol),
            int(nmc),
            float(dang),
            maxout,
            float(badfrac),
            float(cangle),
            float(prob_max),
            int(npolmin),
            float(max_agap),
            float(max_pgap),
            int(selection),
            ctypes.byref(res),
            s_all,
            d_all,
            r_all,
            f_all,
            sl_all,
        )
        output = self._unpack_result(res, npol, 0)
        nf = output["nf"]
        output["strike"] = s_all[:nf].copy()
        output["dip"] = d_all[:nf].copy()
        output["rake"] = r_all[:nf].copy()
        output["faults"] = f_all.reshape(maxout, 3)[:nf, :].T.copy()
        output["slips"] = sl_all.reshape(maxout, 3)[:nf, :].T.copy()
        return output

    @_synchronized
    def run_event_amp(
        self,
        p_azi_mc,
        p_the_mc,
        p_pol,
        sp_amp,
        npsta,
        nmc,
        dang,
        maxout,
        badfrac,
        qbadfac,
        cangle,
        prob_max,
        npolmin,
        max_agap,
        max_pgap,
        selection=0,
        nmismax=None,
    ):
        azi = _as_c_double_array(p_azi_mc).ravel(order="F")
        the = _as_c_double_array(p_the_mc).ravel(order="F")
        pol = _as_c_int32_array(p_pol)
        amp = _as_c_double_array(sp_amp)
        res = _CEventResult()
        maxout = int(maxout)
        s_all = np.zeros(maxout, np.float64)
        d_all = np.zeros(maxout, np.float64)
        r_all = np.zeros(maxout, np.float64)
        f_all = np.zeros(3 * maxout, np.float64)
        sl_all = np.zeros(3 * maxout, np.float64)
        self._lib.cnchash_run_event_amp(
            azi,
            the,
            pol,
            amp,
            int(npsta),
            int(nmc),
            float(dang),
            maxout,
            float(badfrac),
            float(qbadfac),
            float(cangle),
            float(prob_max),
            int(npolmin),
            float(max_agap),
            float(max_pgap),
            int(selection),
            -1 if nmismax is None else int(nmismax),
            ctypes.byref(res),
            s_all,
            d_all,
            r_all,
            f_all,
            sl_all,
        )
        npol = int(np.sum(pol != 0))
        nspr = int(np.sum(amp != 0.0))
        output = self._unpack_result(res, npol, nspr)
        nf = output["nf"]
        output["strike"] = s_all[:nf].copy()
        output["dip"] = d_all[:nf].copy()
        output["rake"] = r_all[:nf].copy()
        output["faults"] = f_all.reshape(maxout, 3)[:nf, :].T.copy()
        output["slips"] = sl_all.reshape(maxout, 3)[:nf, :].T.copy()
        return output

    # ---- batch --------------------------------------------------------------

    @_synchronized
    def run_batch(
        self,
        events,
        dang,
        maxout,
        badfrac,
        cangle,
        prob_max,
        npolmin,
        max_agap,
        max_pgap,
        selection=0,
    ):
        events = list(events)
        nevent = len(events)
        maxout = int(maxout)
        if nevent == 0:
            return []
        npol_arr = np.zeros(nevent, np.int32)
        nmc_arr = np.zeros(nevent, np.int32)
        offsets = np.zeros(nevent + 1, np.int32)
        pick_offsets = np.zeros(nevent + 1, np.int32)
        azi_chunks, the_chunks, pol_chunks, qual_chunks = [], [], [], []
        for i, (azi, the, pol, qual) in enumerate(events):
            azi = _as_c_double_array(azi)
            the = _as_c_double_array(the)
            pol = _as_c_int32_array(pol)
            qual = _as_c_int32_array(qual)
            npol = pol.size
            nmc = azi.size // npol
            npol_arr[i] = npol
            nmc_arr[i] = nmc
            offsets[i + 1] = offsets[i] + npol * nmc
            pick_offsets[i + 1] = pick_offsets[i] + npol
            azi_chunks.append(azi.ravel(order="F"))
            the_chunks.append(the.ravel(order="F"))
            pol_chunks.append(pol.ravel())
            qual_chunks.append(qual.ravel())
        azi_flat = np.concatenate(azi_chunks)
        the_flat = np.concatenate(the_chunks)
        pol_flat = np.concatenate(pol_chunks)
        qual_flat = np.concatenate(qual_chunks)
        cres = (_CEventResult * nevent)()
        s_all = np.zeros(maxout * nevent, np.float64)
        d_all = np.zeros(maxout * nevent, np.float64)
        r_all = np.zeros(maxout * nevent, np.float64)
        f_all = np.zeros(3 * maxout * nevent, np.float64)
        sl_all = np.zeros(3 * maxout * nevent, np.float64)
        self._lib.cnchash_run_batch(
            nevent,
            npol_arr,
            nmc_arr,
            offsets,
            pick_offsets,
            azi_flat,
            the_flat,
            pol_flat,
            qual_flat,
            float(dang),
            maxout,
            float(badfrac),
            float(cangle),
            float(prob_max),
            int(npolmin),
            float(max_agap),
            float(max_pgap),
            int(selection),
            cres,
            s_all,
            d_all,
            r_all,
            f_all,
            sl_all,
        )
        results = []
        for i in range(nevent):
            npol = int(npol_arr[i])
            output = self._unpack_result(cres[i], npol, 0)
            nf = output["nf"]
            sl = slice(i * maxout, i * maxout + nf)
            output["strike"] = s_all[sl].copy()
            output["dip"] = d_all[sl].copy()
            output["rake"] = r_all[sl].copy()
            output["faults"] = f_all.reshape(nevent, maxout, 3)[i, :nf, :].T.copy()
            output["slips"] = sl_all.reshape(nevent, maxout, 3)[i, :nf, :].T.copy()
            results.append(output)
        return results

    @_synchronized
    def run_batch_amp(
        self,
        events,
        dang,
        maxout,
        badfrac,
        qbadfac,
        cangle,
        prob_max,
        npolmin,
        max_agap,
        max_pgap,
        selection=0,
        nmismax=None,
    ):
        events = list(events)
        nevent = len(events)
        maxout = int(maxout)
        if nevent == 0:
            return []
        npol_arr = np.zeros(nevent, np.int32)
        nmc_arr = np.zeros(nevent, np.int32)
        offsets = np.zeros(nevent + 1, np.int32)
        pick_offsets = np.zeros(nevent + 1, np.int32)
        azi_chunks, the_chunks, pol_chunks, amp_chunks = [], [], [], []
        for i, (azi, the, pol, amp) in enumerate(events):
            azi = _as_c_double_array(azi)
            the = _as_c_double_array(the)
            pol = _as_c_int32_array(pol)
            amp = _as_c_double_array(amp)
            npsta = pol.size
            nmc = azi.size // npsta
            npol_arr[i] = npsta
            nmc_arr[i] = nmc
            offsets[i + 1] = offsets[i] + npsta * nmc
            pick_offsets[i + 1] = pick_offsets[i] + npsta
            azi_chunks.append(azi.ravel(order="F"))
            the_chunks.append(the.ravel(order="F"))
            pol_chunks.append(pol.ravel())
            amp_chunks.append(amp.ravel())
        azi_flat = np.concatenate(azi_chunks)
        the_flat = np.concatenate(the_chunks)
        pol_flat = np.concatenate(pol_chunks)
        amp_flat = np.concatenate(amp_chunks)
        cres = (_CEventResult * nevent)()
        s_all = np.zeros(maxout * nevent, np.float64)
        d_all = np.zeros(maxout * nevent, np.float64)
        r_all = np.zeros(maxout * nevent, np.float64)
        f_all = np.zeros(3 * maxout * nevent, np.float64)
        sl_all = np.zeros(3 * maxout * nevent, np.float64)
        self._lib.cnchash_run_batch_amp(
            nevent,
            npol_arr,
            nmc_arr,
            offsets,
            pick_offsets,
            azi_flat,
            the_flat,
            pol_flat,
            amp_flat,
            float(dang),
            maxout,
            float(badfrac),
            float(qbadfac),
            float(cangle),
            float(prob_max),
            int(npolmin),
            float(max_agap),
            float(max_pgap),
            int(selection),
            -1 if nmismax is None else int(nmismax),
            cres,
            s_all,
            d_all,
            r_all,
            f_all,
            sl_all,
        )
        results = []
        for i in range(nevent):
            lo = int(pick_offsets[i])
            hi = int(pick_offsets[i + 1])
            npol = int(np.count_nonzero(pol_flat[lo:hi]))
            nspr = int(np.count_nonzero(amp_flat[lo:hi]))
            output = self._unpack_result(cres[i], npol, nspr)
            nf = output["nf"]
            sl = slice(i * maxout, i * maxout + nf)
            output["strike"] = s_all[sl].copy()
            output["dip"] = d_all[sl].copy()
            output["rake"] = r_all[sl].copy()
            output["faults"] = f_all.reshape(nevent, maxout, 3)[i, :nf, :].T.copy()
            output["slips"] = sl_all.reshape(nevent, maxout, 3)[i, :nf, :].T.copy()
            results.append(output)
        return results

    # ---- low-level pieces ---------------------------------------------------

    @_synchronized
    def get_rotation_grid(self, dang):
        """Rotation grid (b1/b2/b3 column-major) for a grid spacing.

        Returns a dict with keys nrot/b1/b2/b3 (arrays of shape (3, nrot)).
        """
        import math

        # Size generously: nrot grows roughly as dang^-2 (31032 at dang=5).
        max_rot = max(31032, int(math.ceil(31032 * (5.0 / dang) ** 2 * 2)))
        nrot_out = ctypes.c_int32()
        b1 = np.zeros(3 * max_rot, np.float64)
        b2 = np.zeros(3 * max_rot, np.float64)
        b3 = np.zeros(3 * max_rot, np.float64)
        self._lib.cnchash_get_rotation_grid(
            float(dang), ctypes.byref(nrot_out), b1, b2, b3, max_rot
        )
        n = int(nrot_out.value)
        return {
            "nrot": n,
            "b1": b1.reshape(max_rot, 3)[:n, :].T.copy(),
            "b2": b2.reshape(max_rot, 3)[:n, :].T.copy(),
            "b3": b3.reshape(max_rot, 3)[:n, :].T.copy(),
        }

    @_synchronized
    def get_gap(self, npol, p_azi, p_the):
        magap, mpgap = ctypes.c_double(), ctypes.c_double()
        self._lib.cnchash_get_gap(
            int(npol),
            _as_c_double_array(p_azi),
            _as_c_double_array(p_the),
            ctypes.byref(magap),
            ctypes.byref(mpgap),
        )
        return float(magap.value), float(mpgap.value)

    @_synchronized
    def get_misfit(self, npol, p_azi, p_the, p_pol, p_qual, strike, dip, rake):
        mfrac, stdr = ctypes.c_double(), ctypes.c_double()
        self._lib.cnchash_get_misfit(
            int(npol),
            _as_c_double_array(p_azi),
            _as_c_double_array(p_the),
            _as_c_int32_array(p_pol),
            _as_c_int32_array(p_qual),
            float(strike),
            float(dip),
            float(rake),
            ctypes.byref(mfrac),
            ctypes.byref(stdr),
        )
        return float(mfrac.value), float(stdr.value)

    @_synchronized
    def get_misfit_amp(self, npol, p_azi, p_the, sp_amp, p_pol, strike, dip, rake):
        mfrac, mavg, stdr = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
        self._lib.cnchash_get_misfit_amp(
            int(npol),
            _as_c_double_array(p_azi),
            _as_c_double_array(p_the),
            _as_c_double_array(sp_amp),
            _as_c_int32_array(p_pol),
            float(strike),
            float(dip),
            float(rake),
            ctypes.byref(mfrac),
            ctypes.byref(mavg),
            ctypes.byref(stdr),
        )
        return float(mfrac.value), float(mavg.value), float(stdr.value)

    @_synchronized
    def mech_prob(self, nf, faults, slips, cangle, prob_max):
        f = _as_c_double_array(faults).ravel(order="F")
        s = _as_c_double_array(slips).ravel(order="F")
        nsltn = ctypes.c_int32()
        sa = np.zeros(5, np.float64)
        da = np.zeros(5, np.float64)
        ra = np.zeros(5, np.float64)
        pb = np.zeros(5, np.float64)
        rd = np.zeros((5, 2), np.float64)  # C [solution][plane]; Fortran (2,5) col-major
        self._lib.cnchash_mech_prob(
            int(nf),
            f,
            s,
            float(cangle),
            float(prob_max),
            ctypes.byref(nsltn),
            sa,
            da,
            ra,
            pb,
            rd,
        )
        n = int(nsltn.value)
        return {
            "nsltn": n,
            "strike_avg": sa[:n].copy(),
            "dip_avg": da[:n].copy(),
            "rake_avg": ra[:n].copy(),
            "prob": pb[:n].copy(),
            "rms_diff": rd[:n, :].T.copy(),
        }

    @_synchronized
    def mech_rot(self, n1, s1, n2, s2):
        """Minimum rotation angle (degrees) between two mechanisms.

        Each mechanism is given by its fault-normal (n) and slip (s)
        3-vectors; the result is the Kagan angle between the two
        mechanisms (original MECH_ROT).
        """
        a1 = _as_c_double_array(np.atleast_1d(n1).ravel())
        b1 = _as_c_double_array(np.atleast_1d(s1).ravel())
        a2 = _as_c_double_array(np.atleast_1d(n2).ravel())
        b2 = _as_c_double_array(np.atleast_1d(s2).ravel())
        if not (a1.size == 3 and b1.size == 3 and a2.size == 3 and b2.size == 3):
            raise ValueError("n1, s1, n2, s2 must be 3-vectors")
        rota = ctypes.c_double()
        self._lib.cnchash_mech_rot(a1, b1, a2, b2, ctypes.byref(rota))
        return float(rota.value)

    @_synchronized
    def build_velocity_table(self, depth, velocity, params=None):
        if params is None:
            params = _default_table_params()
        depth = _as_c_double_array(depth)
        velocity = _as_c_double_array(velocity)
        ndel_max = int((params["del2"] - params["del1"]) / params["del3"]) + 1
        ndep_max = int((params["dep2"] - params["dep1"]) / params["dep3"]) + 1
        table = np.zeros(ndel_max * ndep_max, np.float64)
        delttab = np.zeros(ndel_max, np.float64)
        deptab = np.zeros(ndep_max, np.float64)
        ndel, ndep = ctypes.c_int32(), ctypes.c_int32()
        self._lib.cnchash_build_velocity_table(
            depth,
            velocity,
            int(depth.size),
            float(params["del1"]),
            float(params["del2"]),
            float(params["del3"]),
            float(params["dep1"]),
            float(params["dep2"]),
            float(params["dep3"]),
            float(params["pmin"]),
            int(params["nump"]),
            table,
            delttab,
            deptab,
            int(ndel_max),
            int(ndep_max),
            ctypes.byref(ndel),
            ctypes.byref(ndep),
        )
        ndel, ndep = int(ndel.value), int(ndep.value)
        table_data = {
            "table": table.reshape(ndel_max, ndep_max, order="F")[:ndel, :ndep].copy(),
            "delttab": delttab[:ndel].copy(),
            "deptab": deptab[:ndep].copy(),
            "ndel": ndel,
            "ndep": ndep,
        }
        self._table_ip += 1
        table_data["ip"] = self._table_ip
        self._tables[self._table_ip] = table_data
        return table_data

    @_synchronized
    def get_tts(self, ip, del_dist, qdep):
        table = self._tables.get(ip)
        if table is None:
            return 999.0, -1
        ndel = table["ndel"]
        ndep = table["ndep"]
        tt, iflag = ctypes.c_double(), ctypes.c_int32()
        self._lib.cnchash_get_tts(
            table["table"].ravel(order="F"),
            table["delttab"],
            table["deptab"],
            int(ndel),
            int(ndep),
            float(del_dist),
            float(qdep),
            ctypes.byref(tt),
            ctypes.byref(iflag),
        )
        return float(tt.value), int(iflag.value)


def _default_table_params():
    """Default takeoff-table generation parameters (Python front-end)."""
    return {
        "del1": 0.0,  # Minimum distance (km)
        "del2": 120.0,  # Maximum distance (km)
        "del3": 1.0,  # Distance step (km)
        "dep1": 0.0,  # Minimum depth (km)
        "dep2": 35.0,  # Maximum depth (km)
        "dep3": 1.0,  # Depth step (km)
        "pmin": 0.0,  # Minimum ray parameter
        "nump": 1000,  # Number of ray parameters
    }
