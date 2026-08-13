"""Numba reference backend.

Pure Python/Numba implementation of the HASH computation contract.
Always available; serves as the fallback when the native backend is
not installed and as the regression oracle for the Fortran backend.
"""

import numpy as np

from .. import core, uncertainty
from ..velocity import _default_table_params, make_table_from_model
from .base import HashBackend

_QUALITY_LETTERS = "ABCDEF"


def _fail_result(quality="F", npol=0, nspr=0):
    """Result dict for events that fail a pre-check (mirrors driver.py)."""
    return {
        "success": False,
        "strike_avg": 999.0,
        "dip_avg": 99.0,
        "rake_avg": 999.0,
        "var_est": [99.0, 99.0],
        "mfrac": 0.99,
        "mavg": 0.99,
        "quality": quality,
        "prob": 0.0,
        "stdr": 0.0,
        "nmult": 0,
        "nout2": 0,
        "nout1": 0,
        "nf": 0,
        "npol": npol,
        "nspr": nspr,
    }


def _rate_quality(prob, var_avg, mfrac, stdr, mavg=None, has_amp=False):
    """HASH quality rating; returns a letter A..D (or the input for failures)."""
    if has_amp:
        if prob > 0.8 and var_avg <= 25.0 and mfrac <= 0.15 and mavg <= 0.2 and stdr >= 0.5:
            return "A"
        if prob > 0.6 and var_avg <= 35.0 and mfrac <= 0.2 and mavg <= 0.3 and stdr >= 0.4:
            return "B"
        if prob > 0.5 and var_avg <= 45.0 and mfrac <= 0.3 and mavg <= 0.4 and stdr >= 0.3:
            return "C"
        return "D"
    if prob > 0.8 and var_avg <= 25.0 and mfrac <= 0.15 and stdr >= 0.5:
        return "A"
    if prob > 0.6 and var_avg <= 35.0 and mfrac <= 0.2 and stdr >= 0.4:
        return "B"
    if prob > 0.5 and var_avg <= 45.0 and mfrac <= 0.3 and stdr >= 0.3:
        return "C"
    return "D"


class NumbaBackend(HashBackend):
    """Reference backend built on the existing Numba kernels."""

    name = "numba"
    available = True

    def __init__(self):
        self._tables = {}
        self._table_ip = 0

    def info(self):
        return {
            "backend": "numba",
            "openmp": False,
            "threads": core.get_numba_threads(),
            "compiler": "Numba JIT",
            "precision": "float64",
        }

    def set_num_threads(self, num_threads):
        return core.set_numba_threads(num_threads)

    def get_num_threads(self):
        return core.get_numba_threads()

    # ---- per-event pipeline -------------------------------------------------

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
        p_azi_mc = np.ascontiguousarray(p_azi_mc, dtype=np.float64)
        p_the_mc = np.ascontiguousarray(p_the_mc, dtype=np.float64)
        p_pol = np.ascontiguousarray(p_pol, dtype=np.int32)
        p_qual = np.ascontiguousarray(p_qual, dtype=np.int32)

        if npol < npolmin:
            return _fail_result("F", npol=npol)
        magap, mpgap = core.get_gap(npol, p_azi_mc[:, 0], p_the_mc[:, 0])
        if magap > max_agap or mpgap > max_pgap:
            return _fail_result("E", npol=npol)

        nmismax = max(int(npol * badfrac), 2)
        nextra = max(int(npol * badfrac * 0.5), 2)
        res = core.focalmc(
            p_azi_mc, p_the_mc, p_pol, p_qual, npol, nmc, dang, maxout, nextra, nmismax
        )
        nf = res["nf"]
        if nf == 0:
            return _fail_result("F", npol=npol)

        pr = uncertainty.mech_prob(nf, res["faults"], res["slips"], cangle, prob_max)
        nsltn = pr["nsltn"]
        if nsltn > 0:
            mfrac, stdr = core.get_misfit(
                npol,
                p_azi_mc[:, 0],
                p_the_mc[:, 0],
                p_pol,
                p_qual,
                pr["strike_avg"][0],
                pr["dip_avg"][0],
                pr["rake_avg"][0],
            )
        else:
            mfrac, stdr = 0.99, 0.0

        var_avg = (pr["rms_diff"][0, 0] + pr["rms_diff"][1, 0]) / 2.0 if nsltn > 0 else 99.0
        quality = _rate_quality(pr["prob"][0] if nsltn > 0 else 0.0, var_avg, mfrac, stdr)

        output = {
            "success": True,
            "strike_avg": pr["strike_avg"][0] if nsltn > 0 else res["strike"][0],
            "dip_avg": pr["dip_avg"][0] if nsltn > 0 else res["dip"][0],
            "rake_avg": pr["rake_avg"][0] if nsltn > 0 else res["rake"][0],
            "var_est": [pr["rms_diff"][0, 0], pr["rms_diff"][1, 0]] if nsltn > 0 else [99.0, 99.0],
            "mfrac": mfrac,
            "mavg": 0.99,
            "quality": quality,
            "prob": pr["prob"][0] if nsltn > 0 else 0.0,
            "stdr": stdr,
            "nmult": nsltn,
            "nout2": nf,
            "nout1": min(nf, maxout),
            "nf": nf,
            "npol": npol,
            "nspr": 0,
        }
        if nsltn > 1:
            output["strike_avg"] = pr["strike_avg"][:nsltn]
            output["dip_avg"] = pr["dip_avg"][:nsltn]
            output["rake_avg"] = pr["rake_avg"][:nsltn]
            output["prob"] = pr["prob"][:nsltn]
            output["rms_diff"] = pr["rms_diff"][:, :nsltn]
            output["quality"] = [quality] * nsltn
        output["faults"] = res["faults"]
        output["slips"] = res["slips"]
        output["strike"] = res["strike"]
        output["dip"] = res["dip"]
        output["rake"] = res["rake"]
        return output

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
    ):
        from .. import amp_subs

        p_azi_mc = np.ascontiguousarray(p_azi_mc, dtype=np.float64)
        p_the_mc = np.ascontiguousarray(p_the_mc, dtype=np.float64)
        p_pol = np.ascontiguousarray(p_pol, dtype=np.int32)
        sp_amp = np.ascontiguousarray(sp_amp, dtype=np.float64)

        npol = int(np.sum(p_pol != 0))
        nspr = int(np.sum(sp_amp != 0.0))
        if npol < npolmin:
            return _fail_result("F", npol=npol, nspr=nspr)
        magap, mpgap = core.get_gap(npsta, p_azi_mc[:, 0], p_the_mc[:, 0])
        if magap > max_agap or mpgap > max_pgap:
            return _fail_result("E", npol=npol, nspr=nspr)

        nmismax = max(int(npol * 0.1), 2)
        nextra = max(int(npol * 0.05), 0)
        qmismax = nspr * qbadfac
        qextra = max(nspr * qbadfac * 0.5, 2.0)
        res = amp_subs.focalamp_mc(
            p_azi_mc, p_the_mc, sp_amp, p_pol, npsta, nmc, dang, maxout, nextra, nmismax, qextra, qmismax
        )
        nf = res["nf"]
        if nf == 0:
            return _fail_result("F", npol=npol, nspr=nspr)

        pr = uncertainty.mech_prob(nf, res["faults"], res["slips"], cangle, prob_max)
        nsltn = pr["nsltn"]
        if nsltn > 0:
            mfrac, mavg, stdr = amp_subs.get_misf_amp(
                npsta,
                p_azi_mc[:, 0],
                p_the_mc[:, 0],
                sp_amp,
                p_pol,
                pr["strike_avg"][0],
                pr["dip_avg"][0],
                pr["rake_avg"][0],
            )
        else:
            mfrac, mavg, stdr = 0.99, 0.99, 0.0

        var_avg = (pr["rms_diff"][0, 0] + pr["rms_diff"][1, 0]) / 2.0 if nsltn > 0 else 99.0
        quality = _rate_quality(
            pr["prob"][0] if nsltn > 0 else 0.0, var_avg, mfrac, stdr, mavg, has_amp=True
        )

        output = {
            "success": True,
            "strike_avg": pr["strike_avg"][0] if nsltn > 0 else res["strike"][0],
            "dip_avg": pr["dip_avg"][0] if nsltn > 0 else res["dip"][0],
            "rake_avg": pr["rake_avg"][0] if nsltn > 0 else res["rake"][0],
            "var_est": [pr["rms_diff"][0, 0], pr["rms_diff"][1, 0]] if nsltn > 0 else [99.0, 99.0],
            "mfrac": mfrac,
            "mavg": mavg,
            "quality": quality,
            "prob": pr["prob"][0] if nsltn > 0 else 0.0,
            "stdr": stdr,
            "nmult": nsltn,
            "nout2": nf,
            "nout1": min(nf, maxout),
            "nf": nf,
            "npol": npol,
            "nspr": nspr,
        }
        if nsltn > 1:
            output["strike_avg"] = pr["strike_avg"][:nsltn]
            output["dip_avg"] = pr["dip_avg"][:nsltn]
            output["rake_avg"] = pr["rake_avg"][:nsltn]
            output["prob"] = pr["prob"][:nsltn]
            output["rms_diff"] = pr["rms_diff"][:, :nsltn]
            output["quality"] = [quality] * nsltn
        output["faults"] = res["faults"]
        output["slips"] = res["slips"]
        output["strike"] = res["strike"]
        output["dip"] = res["dip"]
        output["rake"] = res["rake"]
        return output

    # ---- low-level pieces ---------------------------------------------------

    def get_gap(self, npol, p_azi, p_the):
        return core.get_gap(npol, np.ascontiguousarray(p_azi, np.float64), np.ascontiguousarray(p_the, np.float64))

    def get_misfit(self, npol, p_azi, p_the, p_pol, p_qual, strike, dip, rake):
        return core.get_misfit(
            npol,
            np.ascontiguousarray(p_azi, np.float64),
            np.ascontiguousarray(p_the, np.float64),
            np.ascontiguousarray(p_pol, np.int32),
            np.ascontiguousarray(p_qual, np.int32),
            strike,
            dip,
            rake,
        )

    def get_misfit_amp(self, npol, p_azi, p_the, sp_amp, p_pol, strike, dip, rake):
        from .. import amp_subs

        return amp_subs.get_misf_amp(
            npol,
            np.ascontiguousarray(p_azi, np.float64),
            np.ascontiguousarray(p_the, np.float64),
            np.ascontiguousarray(sp_amp, np.float64),
            np.ascontiguousarray(p_pol, np.int32),
            strike,
            dip,
            rake,
        )

    def mech_prob(self, nf, faults, slips, cangle, prob_max):
        return uncertainty.mech_prob(nf, faults, slips, cangle, prob_max)

    def build_velocity_table(self, depth, velocity, params=None):
        if params is None:
            params = _default_table_params()
        table_data = make_table_from_model(depth, velocity, params)
        self._table_ip += 1
        self._tables[self._table_ip] = table_data
        return table_data

    def get_tts(self, ip, del_dist, qdep):
        table = self._tables.get(ip)
        if table is None:
            return 999.0, -1
        delttab = table["delttab"]
        deptab = table["deptab"]
        tt = table["table"]
        ndel = table["ndel"]
        ndep = table["ndep"]

        if qdep < deptab[0] or qdep > deptab[-1]:
            return 999.0, -1
        id1, id2 = 0, ndep - 1
        for idx in range(1, ndep):
            if deptab[idx] >= qdep:
                id2 = idx
                id1 = idx - 1
                break
        ix1, ix2 = 0, ndel - 1
        for ix in range(1, ndel):
            if delttab[ix] >= del_dist:
                ix1 = ix - 1
                ix2 = ix
                break
        if (
            tt[ix1, id1] != 0.0
            and tt[ix1, id2] != 0.0
            and tt[ix2, id1] != 0.0
            and tt[ix2, id2] != 0.0
            and delttab[ix2] >= del_dist
        ):
            xfrac = (del_dist - delttab[ix1]) / (delttab[ix2] - delttab[ix1])
            t1 = tt[ix1, id1] + xfrac * (tt[ix2, id1] - tt[ix1, id1])
            t2 = tt[ix1, id2] + xfrac * (tt[ix2, id2] - tt[ix1, id2])
            dfrac = (qdep - deptab[id1]) / (deptab[id2] - deptab[id1])
            return t1 + dfrac * (t2 - t1), 0
        return self._get_tts_extrapolate(table, del_dist, qdep, id1, id2, ix1, ix2)

    @staticmethod
    def _get_tts_extrapolate(table, del_dist, qdep, id1, id2, ix1, ix2):
        delttab = table["delttab"]
        deptab = table["deptab"]
        tt = table["table"]
        ndel = table["ndel"]
        xoffmin1, xoffmin2 = 999.0, 999.0
        ixbest1, ixbest2 = -1, -1
        for ix in range(1, ndel):
            if tt[ix - 1, id1] != 0.0 and tt[ix, id1] != 0.0:
                xoff = abs(0.5 * (delttab[ix - 1] + delttab[ix]) - del_dist)
                if xoff < xoffmin1:
                    xoffmin1 = xoff
                    ixbest1 = ix
            if tt[ix - 1, id2] != 0.0 and tt[ix, id2] != 0.0:
                xoff = abs(0.5 * (delttab[ix - 1] + delttab[ix]) - del_dist)
                if xoff < xoffmin2:
                    xoffmin2 = xoff
                    ixbest2 = ix
        if ixbest1 == -1 or ixbest2 == -1:
            return 999.0, -1
        xfrac = (del_dist - delttab[ixbest1 - 1]) / (delttab[ixbest1] - delttab[ixbest1 - 1])
        tt1 = tt[ixbest1 - 1, id1] + xfrac * (tt[ixbest1, id1] - tt[ixbest1 - 1, id1])
        xfrac = (del_dist - delttab[ixbest2 - 1]) / (delttab[ixbest2] - delttab[ixbest2 - 1])
        tt2 = tt[ixbest2 - 1, id2] + xfrac * (tt[ixbest2, id2] - tt[ixbest2 - 1, id2])
        dfrac = (qdep - deptab[id1]) / (deptab[id2] - deptab[id1])
        return tt1 + dfrac * (tt2 - tt1), 1
