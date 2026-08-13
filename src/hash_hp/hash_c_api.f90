!> C ABI for the HASH-HP native backend.
!>
!> Exposes scalar and array interfaces only (double*, int32*, sizes) so
!> the ABI is stable and callable from ctypes without f2py. Fortran
!> derived types, allocatables and character data never cross this
!> boundary.
module hash_c_api
    use iso_c_binding, only: c_int32_t, c_double
    use hash_kinds, only: rk, ik
    use hash_runtime, only: hash_set_num_threads, hash_get_max_threads, &
                            ensure_rotation_grid_cached, grid_cache
    use hash_geometry, only: pi
    use hash_focalmc, only: focalmc, focalmc_workspace_t
    use hash_amplitude, only: focalamp_mc, get_misf_amp, amplitude_workspace_t
    use hash_misfit, only: get_gap, get_misfit
    use hash_uncertainty, only: mech_prob, mech_rot
    use hash_velocity, only: build_takeoff_table, get_tts
    use hash_batch, only: run_event, run_event_amp, run_batch, run_batch_amp, &
                          event_result_t, max_solutions
    implicit none

    !> Mirror of event_result_t with the BIND(C) attribute (no default
    !> initialization allowed in C-interoperable derived types).
    type, bind(C) :: c_event_result_t
        integer(c_int32_t) :: success
        integer(c_int32_t) :: nf
        integer(c_int32_t) :: nout1
        integer(c_int32_t) :: nout2
        integer(c_int32_t) :: nsltn
        integer(c_int32_t) :: quality
        real(c_double) :: strike_avg(5)
        real(c_double) :: dip_avg(5)
        real(c_double) :: rake_avg(5)
        real(c_double) :: prob(5)
        real(c_double) :: var_est(2, 5)
        real(c_double) :: mfrac
        real(c_double) :: mavg
        real(c_double) :: stdr
    end type c_event_result_t

    private
    public :: c_event_result_t

contains

    !> Convert an internal event result to its C-visible mirror.
    subroutine copy_result(fr, cr)
        type(event_result_t), intent(in) :: fr
        type(c_event_result_t), intent(out) :: cr
        cr%success = int(fr%success, c_int32_t)
        cr%nf = int(fr%nf, c_int32_t)
        cr%nout1 = int(fr%nout1, c_int32_t)
        cr%nout2 = int(fr%nout2, c_int32_t)
        cr%nsltn = int(fr%nsltn, c_int32_t)
        cr%quality = int(fr%quality, c_int32_t)
        cr%strike_avg = real(fr%strike_avg, c_double)
        cr%dip_avg = real(fr%dip_avg, c_double)
        cr%rake_avg = real(fr%rake_avg, c_double)
        cr%prob = real(fr%prob, c_double)
        cr%var_est = real(fr%var_est, c_double)
        cr%mfrac = real(fr%mfrac, c_double)
        cr%mavg = real(fr%mavg, c_double)
        cr%stdr = real(fr%stdr, c_double)
    end subroutine copy_result

    !> Set OpenMP thread count.
    subroutine cnchash_set_num_threads(nt) bind(C, name="cnchash_set_num_threads")
        integer(c_int32_t), value :: nt
        call hash_set_num_threads(int(nt, ik))
    end subroutine cnchash_set_num_threads

    !> Current OpenMP max thread count.
    function cnchash_get_max_threads() result(nt) bind(C, name="cnchash_get_max_threads")
        integer(c_int32_t) :: nt
        nt = int(hash_get_max_threads(), c_int32_t)
    end function cnchash_get_max_threads

    !> Version string (major.minor.patch as three ints).
    subroutine cnchash_version(ver_major, ver_minor, ver_patch) &
        bind(C, name="cnchash_version")
        integer(c_int32_t), intent(out) :: ver_major, ver_minor, ver_patch
        ver_major = int(1, c_int32_t)
        ver_minor = int(0, c_int32_t)
        ver_patch = int(0, c_int32_t)
    end subroutine cnchash_version

    !> Full polarity-only pipeline for one event.
    subroutine cnchash_run_event(p_azi_mc, p_the_mc, p_pol, p_qual, npol, nmc, &
                                 dang, maxout, badfrac, cangle, prob_max, &
                                 npolmin, max_agap, max_pgap, selection, res, &
                                 strike_all, dip_all, rake_all, faults_all, &
                                 slips_all) bind(C, name="cnchash_run_event")
        real(c_double), intent(in) :: p_azi_mc(*), p_the_mc(*)
        integer(c_int32_t), intent(in) :: p_pol(*), p_qual(*)
        integer(c_int32_t), value :: npol, nmc, maxout, selection
        real(c_double), value :: dang, badfrac, cangle, prob_max
        integer(c_int32_t), value :: npolmin, max_agap, max_pgap
        type(c_event_result_t), intent(out) :: res
        real(c_double), intent(out) :: strike_all(*), dip_all(*), rake_all(*)
        real(c_double), intent(out) :: faults_all(*), slips_all(*)
        type(focalmc_workspace_t) :: fws
        type(event_result_t) :: fr
        real(rk) :: s_all(maxout), d_all(maxout), r_all(maxout)
        real(rk) :: f_all(3, maxout), sl_all(3, maxout)
        integer(ik) :: n

        call run_event(fws, reshape(p_azi_mc(1:npol * nmc), [npol, nmc]), &
                       reshape(p_the_mc(1:npol * nmc), [npol, nmc]), &
                       int(p_pol(1:npol), ik), int(p_qual(1:npol), ik), &
                       int(npol, ik), int(nmc, ik), real(dang, rk), &
                       int(maxout, ik), real(badfrac, rk), real(cangle, rk), &
                       real(prob_max, rk), int(npolmin, ik), int(max_agap, ik), &
                       int(max_pgap, ik), int(selection, ik), fr, s_all, d_all, &
                       r_all, f_all, sl_all, .true.)

        call copy_result(fr, res)
        do n = 1, maxout
            strike_all(n) = real(s_all(n), c_double)
            dip_all(n) = real(d_all(n), c_double)
            rake_all(n) = real(r_all(n), c_double)
        end do
        ! faults_all/slips_all are (3, maxout) column-major.
        do n = 1, maxout
            faults_all(3 * (n - 1) + 1) = real(f_all(1, n), c_double)
            faults_all(3 * (n - 1) + 2) = real(f_all(2, n), c_double)
            faults_all(3 * (n - 1) + 3) = real(f_all(3, n), c_double)
            slips_all(3 * (n - 1) + 1) = real(sl_all(1, n), c_double)
            slips_all(3 * (n - 1) + 2) = real(sl_all(2, n), c_double)
            slips_all(3 * (n - 1) + 3) = real(sl_all(3, n), c_double)
        end do
    end subroutine cnchash_run_event

    !> Full polarity + S/P amplitude pipeline for one event.
    subroutine cnchash_run_event_amp(p_azi_mc, p_the_mc, p_pol, sp_amp, npsta, &
                                     nmc, dang, maxout, badfrac, qbadfac, &
                                     cangle, prob_max, npolmin, max_agap, &
                                     max_pgap, selection, res, strike_all, &
                                     dip_all, rake_all, faults_all, slips_all) &
        bind(C, name="cnchash_run_event_amp")
        real(c_double), intent(in) :: p_azi_mc(*), p_the_mc(*), sp_amp(*)
        integer(c_int32_t), intent(in) :: p_pol(*)
        integer(c_int32_t), value :: npsta, nmc, maxout, selection
        real(c_double), value :: dang, badfrac, qbadfac, cangle, prob_max
        integer(c_int32_t), value :: npolmin, max_agap, max_pgap
        type(c_event_result_t), intent(out) :: res
        real(c_double), intent(out) :: strike_all(*), dip_all(*), rake_all(*)
        real(c_double), intent(out) :: faults_all(*), slips_all(*)
        type(amplitude_workspace_t) :: aws
        type(event_result_t) :: fr
        real(rk) :: s_all(maxout), d_all(maxout), r_all(maxout)
        real(rk) :: f_all(3, maxout), sl_all(3, maxout)
        integer(ik) :: n

        call run_event_amp(aws, reshape(p_azi_mc(1:npsta * nmc), [npsta, nmc]), &
                           reshape(p_the_mc(1:npsta * nmc), [npsta, nmc]), &
                           int(p_pol(1:npsta), ik), real(sp_amp(1:npsta), rk), &
                           int(npsta, ik), int(nmc, ik), real(dang, rk), &
                           int(maxout, ik), real(badfrac, rk), real(qbadfac, rk), &
                           real(cangle, rk), real(prob_max, rk), int(npolmin, ik), &
                           int(max_agap, ik), int(max_pgap, ik), int(selection, ik), &
                           fr, s_all, d_all, r_all, f_all, sl_all, .true.)

        call copy_result(fr, res)
        do n = 1, maxout
            strike_all(n) = real(s_all(n), c_double)
            dip_all(n) = real(d_all(n), c_double)
            rake_all(n) = real(r_all(n), c_double)
            faults_all(3 * (n - 1) + 1) = real(f_all(1, n), c_double)
            faults_all(3 * (n - 1) + 2) = real(f_all(2, n), c_double)
            faults_all(3 * (n - 1) + 3) = real(f_all(3, n), c_double)
            slips_all(3 * (n - 1) + 1) = real(sl_all(1, n), c_double)
            slips_all(3 * (n - 1) + 2) = real(sl_all(2, n), c_double)
            slips_all(3 * (n - 1) + 3) = real(sl_all(3, n), c_double)
        end do
    end subroutine cnchash_run_event_amp

    !> Batch polarity-only pipeline (CSR-style flat inputs).
    subroutine cnchash_run_batch(nevent, npol_arr, nmc_arr, offsets, pick_offsets, &
                                 azi_flat, the_flat, pol_flat, qual_flat, dang, &
                                 maxout, badfrac, cangle, prob_max, npolmin, &
                                 max_agap, max_pgap, selection, results, strike_all, &
                                 dip_all, rake_all, faults_all, slips_all) &
        bind(C, name="cnchash_run_batch")
        integer(c_int32_t), value :: nevent
        integer(c_int32_t), intent(in) :: npol_arr(*), nmc_arr(*), offsets(*), pick_offsets(*)
        real(c_double), intent(in) :: azi_flat(*), the_flat(*)
        integer(c_int32_t), intent(in) :: pol_flat(*), qual_flat(*)
        real(c_double), value :: dang, badfrac, cangle, prob_max
        integer(c_int32_t), value :: maxout, npolmin, max_agap, max_pgap, selection
        type(c_event_result_t), intent(out) :: results(*)
        real(c_double), intent(out) :: strike_all(*), dip_all(*), rake_all(*)
        real(c_double), intent(out) :: faults_all(*), slips_all(*)
        type(event_result_t) :: fr(nevent)
        real(rk) :: s_all(maxout, nevent), d_all(maxout, nevent), r_all(maxout, nevent)
        real(rk) :: f_all(3, maxout, nevent), sl_all(3, maxout, nevent)
        integer(ik) :: ie, n

        call run_batch(int(nevent, ik), int(npol_arr(1:nevent), ik), &
                       int(nmc_arr(1:nevent), ik), int(offsets(1:nevent + 1), ik), &
                       int(pick_offsets(1:nevent + 1), ik), &
                       real(azi_flat(1:offsets(nevent + 1)), rk), &
                       real(the_flat(1:offsets(nevent + 1)), rk), &
                       int(pol_flat(1:pick_offsets(nevent + 1)), ik), &
                       int(qual_flat(1:pick_offsets(nevent + 1)), ik), &
                       real(dang, rk), int(maxout, ik), real(badfrac, rk), &
                       real(cangle, rk), real(prob_max, rk), int(npolmin, ik), &
                       int(max_agap, ik), int(max_pgap, ik), int(selection, ik), &
                       fr, s_all, d_all, r_all, f_all, sl_all)

        do ie = 1, nevent
            call copy_result(fr(ie), results(ie))
            do n = 1, maxout
                strike_all(n + maxout * (ie - 1)) = real(s_all(n, ie), c_double)
                dip_all(n + maxout * (ie - 1)) = real(d_all(n, ie), c_double)
                rake_all(n + maxout * (ie - 1)) = real(r_all(n, ie), c_double)
                faults_all(3 * (n - 1) + 1 + 3 * maxout * (ie - 1)) = real(f_all(1, n, ie), c_double)
                faults_all(3 * (n - 1) + 2 + 3 * maxout * (ie - 1)) = real(f_all(2, n, ie), c_double)
                faults_all(3 * (n - 1) + 3 + 3 * maxout * (ie - 1)) = real(f_all(3, n, ie), c_double)
                slips_all(3 * (n - 1) + 1 + 3 * maxout * (ie - 1)) = real(sl_all(1, n, ie), c_double)
                slips_all(3 * (n - 1) + 2 + 3 * maxout * (ie - 1)) = real(sl_all(2, n, ie), c_double)
                slips_all(3 * (n - 1) + 3 + 3 * maxout * (ie - 1)) = real(sl_all(3, n, ie), c_double)
            end do
        end do
    end subroutine cnchash_run_batch

    !> Batch polarity + S/P amplitude pipeline.
    subroutine cnchash_run_batch_amp(nevent, npsta_arr, nmc_arr, offsets, &
                                     pick_offsets, azi_flat, the_flat, pol_flat, &
                                     sp_flat, dang, maxout, badfrac, qbadfac, &
                                     cangle, prob_max, npolmin, max_agap, &
                                     max_pgap, selection, results, strike_all, &
                                     dip_all, rake_all, faults_all, slips_all) &
        bind(C, name="cnchash_run_batch_amp")
        integer(c_int32_t), value :: nevent
        integer(c_int32_t), intent(in) :: npsta_arr(*), nmc_arr(*), offsets(*), pick_offsets(*)
        real(c_double), intent(in) :: azi_flat(*), the_flat(*), sp_flat(*)
        integer(c_int32_t), intent(in) :: pol_flat(*)
        real(c_double), value :: dang, badfrac, qbadfac, cangle, prob_max
        integer(c_int32_t), value :: maxout, npolmin, max_agap, max_pgap, selection
        type(c_event_result_t), intent(out) :: results(*)
        real(c_double), intent(out) :: strike_all(*), dip_all(*), rake_all(*)
        real(c_double), intent(out) :: faults_all(*), slips_all(*)
        type(event_result_t) :: fr(nevent)
        real(rk) :: s_all(maxout, nevent), d_all(maxout, nevent), r_all(maxout, nevent)
        real(rk) :: f_all(3, maxout, nevent), sl_all(3, maxout, nevent)
        integer(ik) :: ie, n

        call run_batch_amp(int(nevent, ik), int(npsta_arr(1:nevent), ik), &
                           int(nmc_arr(1:nevent), ik), int(offsets(1:nevent + 1), ik), &
                           int(pick_offsets(1:nevent + 1), ik), &
                           real(azi_flat(1:offsets(nevent + 1)), rk), &
                           real(the_flat(1:offsets(nevent + 1)), rk), &
                           int(pol_flat(1:pick_offsets(nevent + 1)), ik), &
                           real(sp_flat(1:pick_offsets(nevent + 1)), rk), &
                           real(dang, rk), int(maxout, ik), real(badfrac, rk), &
                           real(qbadfac, rk), real(cangle, rk), real(prob_max, rk), &
                           int(npolmin, ik), int(max_agap, ik), int(max_pgap, ik), &
                           int(selection, ik), fr, s_all, d_all, r_all, f_all, sl_all)

        do ie = 1, nevent
            call copy_result(fr(ie), results(ie))
            do n = 1, maxout
                strike_all(n + maxout * (ie - 1)) = real(s_all(n, ie), c_double)
                dip_all(n + maxout * (ie - 1)) = real(d_all(n, ie), c_double)
                rake_all(n + maxout * (ie - 1)) = real(r_all(n, ie), c_double)
                faults_all(3 * (n - 1) + 1 + 3 * maxout * (ie - 1)) = real(f_all(1, n, ie), c_double)
                faults_all(3 * (n - 1) + 2 + 3 * maxout * (ie - 1)) = real(f_all(2, n, ie), c_double)
                faults_all(3 * (n - 1) + 3 + 3 * maxout * (ie - 1)) = real(f_all(3, n, ie), c_double)
                slips_all(3 * (n - 1) + 1 + 3 * maxout * (ie - 1)) = real(sl_all(1, n, ie), c_double)
                slips_all(3 * (n - 1) + 2 + 3 * maxout * (ie - 1)) = real(sl_all(2, n, ie), c_double)
                slips_all(3 * (n - 1) + 3 + 3 * maxout * (ie - 1)) = real(sl_all(3, n, ie), c_double)
            end do
        end do
    end subroutine cnchash_run_batch_amp

    !> Build (if needed) and expose the cached rotation grid.
    !> nrot is returned; b1/b2/b3 buffers (3, max_rot, column-major) are
    !> filled up to nrot entries.
    subroutine cnchash_get_rotation_grid(dang, nrot, b1, b2, b3, max_rot) &
        bind(C, name="cnchash_get_rotation_grid")
        real(c_double), value :: dang
        integer(c_int32_t), intent(out) :: nrot
        real(c_double), intent(out) :: b1(*), b2(*), b3(*)
        integer(c_int32_t), value :: max_rot
        integer(ik) :: i, n
        call ensure_rotation_grid_cached(real(dang, rk))
        n = min(grid_cache%nrot, int(max_rot, ik))
        nrot = int(n, c_int32_t)
        do i = 1, n
            b1(3 * (i - 1) + 1) = real(grid_cache%b1(i, 1), c_double)
            b1(3 * (i - 1) + 2) = real(grid_cache%b1(i, 2), c_double)
            b1(3 * (i - 1) + 3) = real(grid_cache%b1(i, 3), c_double)
            b2(3 * (i - 1) + 1) = real(grid_cache%b2(i, 1), c_double)
            b2(3 * (i - 1) + 2) = real(grid_cache%b2(i, 2), c_double)
            b2(3 * (i - 1) + 3) = real(grid_cache%b2(i, 3), c_double)
            b3(3 * (i - 1) + 1) = real(grid_cache%b3(i, 1), c_double)
            b3(3 * (i - 1) + 2) = real(grid_cache%b3(i, 2), c_double)
            b3(3 * (i - 1) + 3) = real(grid_cache%b3(i, 3), c_double)
        end do
    end subroutine cnchash_get_rotation_grid

    !> Maximum azimuthal/takeoff gaps (wraps GET_GAP).
    subroutine cnchash_get_gap(npol, p_azi, p_the, magap, mpgap) &
        bind(C, name="cnchash_get_gap")
        integer(c_int32_t), value :: npol
        real(c_double), intent(in) :: p_azi(*), p_the(*)
        real(c_double), intent(out) :: magap, mpgap
        real(rk) :: mg, mp
        call get_gap(int(npol, ik), real(p_azi(1:npol), rk), real(p_the(1:npol), rk), mg, mp)
        magap = real(mg, c_double)
        mpgap = real(mp, c_double)
    end subroutine cnchash_get_gap

    !> Misfit fraction and distribution ratio (wraps GET_MISF).
    subroutine cnchash_get_misfit(npol, p_azi, p_the, p_pol, p_qual, str, dip, &
                                  rake, mfrac, stdr) bind(C, name="cnchash_get_misfit")
        integer(c_int32_t), value :: npol
        real(c_double), intent(in) :: p_azi(*), p_the(*)
        integer(c_int32_t), intent(in) :: p_pol(*), p_qual(*)
        real(c_double), value :: str, dip, rake
        real(c_double), intent(out) :: mfrac, stdr
        real(rk) :: mf, sd
        call get_misfit(int(npol, ik), real(p_azi(1:npol), rk), &
                        real(p_the(1:npol), rk), int(p_pol(1:npol), ik), &
                        int(p_qual(1:npol), ik), real(str, rk), real(dip, rk), &
                        real(rake, rk), mf, sd)
        mfrac = real(mf, c_double)
        stdr = real(sd, c_double)
    end subroutine cnchash_get_misfit

    !> Misfit with S/P amplitude ratios (wraps GET_MISF_AMP).
    subroutine cnchash_get_misfit_amp(npol, p_azi, p_the, sp_amp, p_pol, str, &
                                      dip, rake, mfrac, mavg, stdr) &
        bind(C, name="cnchash_get_misfit_amp")
        integer(c_int32_t), value :: npol
        real(c_double), intent(in) :: p_azi(*), p_the(*), sp_amp(*)
        integer(c_int32_t), intent(in) :: p_pol(*)
        real(c_double), value :: str, dip, rake
        real(c_double), intent(out) :: mfrac, mavg, stdr
        real(rk) :: mf, ma, sd
        call get_misf_amp(int(npol, ik), real(p_azi(1:npol), rk), &
                          real(p_the(1:npol), rk), real(sp_amp(1:npol), rk), &
                          int(p_pol(1:npol), ik), real(str, rk), real(dip, rk), &
                          real(rake, rk), mf, ma, sd)
        mfrac = real(mf, c_double)
        mavg = real(ma, c_double)
        stdr = real(sd, c_double)
    end subroutine cnchash_get_misfit_amp

    !> Preferred mechanism(s) from a set of acceptable mechanisms.
    subroutine cnchash_mech_prob(nf, faults, slips, cangle, prob_max, nsltn, &
                                 str_avg, dip_avg, rak_avg, prob, rms_diff) &
        bind(C, name="cnchash_mech_prob")
        integer(c_int32_t), value :: nf
        real(c_double), intent(in) :: faults(*), slips(*)
        real(c_double), value :: cangle, prob_max
        integer(c_int32_t), intent(out) :: nsltn
        real(c_double), intent(out) :: str_avg(5), dip_avg(5), rak_avg(5)
        real(c_double), intent(out) :: prob(5), rms_diff(2, 5)
        real(rk) :: f(3, nf), s(3, nf)
        real(rk) :: sa(5), da(5), ra(5), pb(5), rd(2, 5)
        integer(ik) :: nsl, n, i
        do n = 1, nf
            do i = 1, 3
                f(i, n) = real(faults(3 * (n - 1) + i), rk)
                s(i, n) = real(slips(3 * (n - 1) + i), rk)
            end do
        end do
        call mech_prob(int(nf, ik), f, s, real(cangle, rk), real(prob_max, rk), &
                       nsl, sa, da, ra, pb, rd)
        nsltn = int(nsl, c_int32_t)
        str_avg = real(sa, c_double)
        dip_avg = real(da, c_double)
        rak_avg = real(ra, c_double)
        prob = real(pb, c_double)
        rms_diff = real(rd, c_double)
    end subroutine cnchash_mech_prob

    !> Minimum rotation angle between two mechanisms (wraps MECH_ROT).
    subroutine cnchash_mech_rot(n1, s1, n2, s2, rota) bind(C, name="cnchash_mech_rot")
        real(c_double), intent(in) :: n1(*), s1(*), n2(*), s2(*)
        real(c_double), intent(out) :: rota
        real(rk) :: a1(3), b1(3), a2(3), b2(3), r
        a1 = real(n1(1:3), rk)
        b1 = real(s1(1:3), rk)
        a2 = real(n2(1:3), rk)
        b2 = real(s2(1:3), rk)
        call mech_rot(a1, b1, a2, b2, r)
        rota = real(r, c_double)
    end subroutine cnchash_mech_rot

    !> Build a takeoff-angle table from a velocity model (no I/O).
    subroutine cnchash_build_velocity_table(z, alpha, npts, del1, del2, del3, &
                                            dep1, dep2, dep3, pmin, nump, &
                                            table, delttab, deptab, ndel_max, &
                                            ndep_max, ndel, ndep) &
        bind(C, name="cnchash_build_velocity_table")
        real(c_double), intent(in) :: z(*), alpha(*)
        integer(c_int32_t), value :: npts
        real(c_double), value :: del1, del2, del3, dep1, dep2, dep3, pmin
        integer(c_int32_t), value :: nump, ndel_max, ndep_max
        real(c_double), intent(out) :: table(*), delttab(*), deptab(*)
        integer(c_int32_t), intent(out) :: ndel, ndep
        real(rk) :: tbl(ndel_max, ndep_max)
        real(rk) :: dtab(ndel_max), etab(ndep_max)
        integer(ik) :: nd, np
        tbl = 0.0_rk
        call build_takeoff_table(real(z(1:npts), rk), real(alpha(1:npts), rk), &
                                 int(npts, ik), real(del1, rk), real(del2, rk), &
                                 real(del3, rk), real(dep1, rk), real(dep2, rk), &
                                 real(dep3, rk), real(pmin, rk), int(nump, ik), &
                                 tbl, dtab, etab, nd, np)
        ndel = int(nd, c_int32_t)
        ndep = int(np, c_int32_t)
        do np = 1, int(ndep, ik)
            do nd = 1, int(ndel, ik)
                table(nd + ndel_max * (np - 1)) = real(tbl(nd, np), c_double)
            end do
            deptab(np) = real(etab(np), c_double)
        end do
        do nd = 1, int(ndel, ik)
            delttab(nd) = real(dtab(nd), c_double)
        end do
    end subroutine cnchash_build_velocity_table

    !> Interpolate a takeoff angle from a table (wraps GET_TTS).
    subroutine cnchash_get_tts(table, delttab, deptab, ndel, ndep, del, qdep, &
                               tt, iflag) bind(C, name="cnchash_get_tts")
        real(c_double), intent(in) :: table(*), delttab(*), deptab(*)
        integer(c_int32_t), value :: ndel, ndep
        real(c_double), value :: del, qdep
        real(c_double), intent(out) :: tt
        integer(c_int32_t), intent(out) :: iflag
        real(rk) :: tbl(ndel, ndep), dtb(ndel), epb(ndep), t
        integer(ik) :: i, j, fl
        do j = 1, int(ndep, ik)
            do i = 1, int(ndel, ik)
                tbl(i, j) = real(table(i + ndel * (j - 1)), rk)
            end do
        end do
        do i = 1, int(ndel, ik)
            dtb(i) = real(delttab(i), rk)
        end do
        do j = 1, int(ndep, ik)
            epb(j) = real(deptab(j), rk)
        end do
        call get_tts(tbl, dtb, epb, int(ndel, ik), int(ndep, ik), real(del, rk), &
                     real(qdep, rk), t, fl)
        tt = real(t, c_double)
        iflag = int(fl, c_int32_t)
    end subroutine cnchash_get_tts

end module hash_c_api
