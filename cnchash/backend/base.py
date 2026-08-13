"""Abstract backend interface for CNCHASH.

A backend implements the full HASH computation contract: gap checks,
focal-mechanism grid search (with or without S/P amplitudes), preferred
mechanism extraction, misfit statistics, velocity-model takeoff tables,
and thread control.

All inputs are NumPy arrays with float64/int32 dtypes:

- ``p_azi_mc``, ``p_the_mc``: shape ``(n, nmc)``, azimuths (deg E of N)
  and takeoff angles (deg from vertical) per observation and MC trial
- ``p_pol``: shape ``(n,)`` int32, +1/-1 (0 = no polarity for amp runs)
- ``p_qual``: shape ``(n,)`` int32, 0=impulsive, 1=emergent
- ``sp_amp``: shape ``(n,)`` float64, log10(S/P) ratios, 0 = no data
"""

from abc import ABC, abstractmethod

_QUALITY_LETTERS = "ABCDEF"

MAX_SOLUTIONS = 5


def quality_int_to_letter(code):
    """Map a backend quality code (0=A..5=F) to its letter."""
    return _QUALITY_LETTERS[code]


class HashBackend(ABC):
    """Common interface implemented by both backends."""

    name = "base"
    available = False

    @abstractmethod
    def info(self):
        """Return a dict describing this backend (name, threads, flags)."""

    @abstractmethod
    def set_num_threads(self, num_threads):
        """Set the thread count; return the effective count."""

    @abstractmethod
    def get_num_threads(self):
        """Return the current thread count."""

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
        """Run the full polarity-only pipeline for one event.

        Returns a result dict (see ``run_hash`` for the schema).
        """
        raise NotImplementedError

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
        """Run the full polarity + S/P amplitude pipeline for one event."""
        raise NotImplementedError

    # ---- batch --------------------------------------------------------------

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
        """Run the polarity-only pipeline for many events at once.

        ``events`` is an iterable of ``(p_azi_mc, p_the_mc, p_pol, p_qual)``
        tuples. Returns a list of result dicts in the same order.
        """
        return [
            self.run_event(
                azi,
                the,
                pol,
                qual,
                len(pol),
                azi.shape[1],
                dang,
                maxout,
                badfrac,
                cangle,
                prob_max,
                npolmin,
                max_agap,
                max_pgap,
                selection,
            )
            for azi, the, pol, qual in events
        ]

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
    ):
        """Run the polarity + amplitude pipeline for many events at once."""
        return [
            self.run_event_amp(
                azi,
                the,
                pol,
                amp,
                len(pol),
                azi.shape[1],
                dang,
                maxout,
                badfrac,
                qbadfac,
                cangle,
                prob_max,
                npolmin,
                max_agap,
                max_pgap,
                selection,
            )
            for azi, the, pol, amp in events
        ]

    # ---- low-level pieces (used by tests and the numba reference) -----------

    def get_gap(self, npol, p_azi, p_the):
        """Maximum azimuthal and takeoff gaps (degrees)."""
        raise NotImplementedError

    def get_misfit(self, npol, p_azi, p_the, p_pol, p_qual, strike, dip, rake):
        """(mfrac, stdr) for a given mechanism."""
        raise NotImplementedError

    def get_misfit_amp(self, npol, p_azi, p_the, sp_amp, p_pol, strike, dip, rake):
        """(mfrac, mavg, stdr) for a given mechanism with S/P ratios."""
        raise NotImplementedError

    def mech_prob(self, nf, faults, slips, cangle, prob_max):
        """Preferred mechanism(s) from a set of acceptable mechanisms.

        Returns a dict with keys nsltn/strike_avg/dip_avg/rake_avg/prob/rms_diff.
        """
        raise NotImplementedError

    def build_velocity_table(self, depth, velocity, params=None):
        """Takeoff-angle table from a 1D velocity model.

        Returns a dict with keys table/delttab/deptab/ndel/ndep.
        """
        raise NotImplementedError

    def get_tts(self, ip, del_dist, qdep):
        """Interpolate a takeoff angle from a built table: (tt, iflag)."""
        raise NotImplementedError
