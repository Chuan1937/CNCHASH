!> Full per-event HASH pipeline and batch execution.
!>
!> Combines GET_GAP -> FOCALMC/FOCALAMP_MC -> MECH_PROB -> GET_MISF(_AMP)
!> -> quality rating into one callable, plus batch entries that
!> parallelize over events (with internal rotation parallelism
!> disabled to avoid nested OpenMP).
module hash_batch
    use hash_kinds, only: rk, ik
    use hash_focalmc, only: focalmc, focalmc_workspace_t
    use hash_amplitude, only: focalamp_mc, get_misf_amp, amplitude_workspace_t
    use hash_misfit, only: get_gap, get_misfit
    use hash_uncertainty, only: mech_prob
    use hash_runtime, only: hash_get_max_threads
    implicit none

    integer(ik), parameter, public :: max_solutions = 5

    type, public :: event_result_t
        integer(ik) :: success = 0
        integer(ik) :: nf = 0
        integer(ik) :: nout1 = 0
        integer(ik) :: nout2 = 0
        integer(ik) :: nsltn = 0
        integer(ik) :: quality = 5          ! 0=A,1=B,2=C,3=D,4=E,5=F
        real(rk) :: strike_avg(max_solutions) = 999.0_rk
        real(rk) :: dip_avg(max_solutions) = 99.0_rk
        real(rk) :: rake_avg(max_solutions) = 999.0_rk
        real(rk) :: prob(max_solutions) = 0.0_rk
        real(rk) :: var_est(2, max_solutions) = 99.0_rk
        real(rk) :: mfrac = 0.99_rk
        real(rk) :: mavg = 0.99_rk
        real(rk) :: stdr = 0.0_rk
    end type event_result_t

    type(focalmc_workspace_t) :: fws_serial, fws_thread
    type(amplitude_workspace_t) :: aws_serial, aws_thread
    !$omp threadprivate(fws_thread, aws_thread)

    private
    public :: run_event, run_event_amp, run_batch, run_batch_amp, rate_quality

contains

    !> Quality rating (letters) per HASH driver logic; returns 0..5.
    subroutine rate_quality(prob, var_avg, mfrac, stdr, mavg, has_amp, quality)
        real(rk), intent(in) :: prob, var_avg, mfrac, stdr, mavg
        logical, intent(in) :: has_amp
        integer(ik), intent(out) :: quality
        if (has_amp) then
            if (prob > 0.8_rk .and. var_avg <= 25.0_rk .and. mfrac <= 0.15_rk .and. &
                mavg <= 0.2_rk .and. stdr >= 0.5_rk) then
                quality = 0
            else if (prob > 0.6_rk .and. var_avg <= 35.0_rk .and. mfrac <= 0.2_rk .and. &
                     mavg <= 0.3_rk .and. stdr >= 0.4_rk) then
                quality = 1
            else if (prob > 0.5_rk .and. var_avg <= 45.0_rk .and. mfrac <= 0.3_rk .and. &
                     mavg <= 0.4_rk .and. stdr >= 0.3_rk) then
                quality = 2
            else
                quality = 3
            end if
        else
            if (prob > 0.8_rk .and. var_avg <= 25.0_rk .and. mfrac <= 0.15_rk .and. &
                stdr >= 0.5_rk) then
                quality = 0
            else if (prob > 0.6_rk .and. var_avg <= 35.0_rk .and. mfrac <= 0.2_rk .and. &
                     stdr >= 0.4_rk) then
                quality = 1
            else if (prob > 0.5_rk .and. var_avg <= 45.0_rk .and. mfrac <= 0.3_rk .and. &
                     stdr >= 0.3_rk) then
                quality = 2
            else
                quality = 3
            end if
        end if
    end subroutine rate_quality

    !> Full polarity-only pipeline for one event.
    subroutine run_event(fws, p_azi_mc, p_the_mc, p_pol, p_qual, npol, nmc, &
                         dang, maxout, badfrac, cangle, prob_max, npolmin, &
                         max_agap, max_pgap, selection, res, strike_all, &
                         dip_all, rake_all, faults_all, slips_all, use_omp)
        type(focalmc_workspace_t), intent(inout) :: fws
        real(rk), intent(in) :: p_azi_mc(:, :), p_the_mc(:, :)
        integer(ik), intent(in) :: p_pol(:), p_qual(:)
        integer(ik), intent(in) :: npol, nmc, maxout, selection
        real(rk), intent(in) :: dang, badfrac, cangle, prob_max
        integer(ik), intent(in) :: npolmin, max_agap, max_pgap
        type(event_result_t), intent(out) :: res
        real(rk), intent(out) :: strike_all(:), dip_all(:), rake_all(:)
        real(rk), intent(out) :: faults_all(:, :), slips_all(:, :)
        logical, intent(in) :: use_omp
        real(rk) :: magap, mpgap, var_avg
        integer(ik) :: nmismax, nextra

        res = event_result_t()
        res%strike_avg = 999.0_rk
        res%dip_avg = 99.0_rk
        res%rake_avg = 999.0_rk
        res%var_est = 99.0_rk

        if (npol < npolmin) then
            res%quality = 5
            return
        end if

        call get_gap(npol, p_azi_mc(:, 1), p_the_mc(:, 1), magap, mpgap)
        if (magap > real(max_agap, rk) .or. mpgap > real(max_pgap, rk)) then
            res%quality = 4
            return
        end if

        nmismax = max(int(real(npol, rk) * badfrac), 2)
        nextra = max(int(real(npol, rk) * badfrac * 0.5_rk), 2)

        call focalmc(fws, p_azi_mc, p_the_mc, p_pol, p_qual, npol, nmc, dang, &
                     maxout, nextra, nmismax, selection, res%nf, strike_all, &
                     dip_all, rake_all, faults_all, slips_all, use_omp)

        res%nout2 = min(maxout, res%nf)
        res%nout1 = min(maxout, res%nf)
        if (res%nf == 0) then
            res%quality = 5
            return
        end if

        call mech_prob(res%nout2, faults_all(:, 1:res%nout2), slips_all(:, 1:res%nout2), &
                       cangle, prob_max, res%nsltn, res%strike_avg, res%dip_avg, &
                       res%rake_avg, res%prob, res%var_est)

        if (res%nsltn > 0) then
            call get_misfit(npol, p_azi_mc(:, 1), p_the_mc(:, 1), p_pol, p_qual, &
                            res%strike_avg(1), res%dip_avg(1), res%rake_avg(1), &
                            res%mfrac, res%stdr)
        else
            res%mfrac = 0.99_rk
            res%stdr = 0.0_rk
        end if

        var_avg = (res%var_est(1, 1) + res%var_est(2, 1)) / 2.0_rk
        call rate_quality(res%prob(1), var_avg, res%mfrac, res%stdr, 0.0_rk, &
                          .false., res%quality)
        res%success = 1
    end subroutine run_event

    !> Full polarity + S/P amplitude pipeline for one event.
    subroutine run_event_amp(aws, p_azi_mc, p_the_mc, p_pol, sp_amp, npsta, &
                             nmc, dang, maxout, badfrac, qbadfac, cangle, &
                             prob_max, npolmin, max_agap, max_pgap, selection, &
                             res, strike_all, dip_all, rake_all, faults_all, &
                             slips_all, use_omp)
        type(amplitude_workspace_t), intent(inout) :: aws
        real(rk), intent(in) :: p_azi_mc(:, :), p_the_mc(:, :)
        integer(ik), intent(in) :: p_pol(:)
        real(rk), intent(in) :: sp_amp(:)
        integer(ik), intent(in) :: npsta, nmc, maxout, selection
        real(rk), intent(in) :: dang, badfrac, qbadfac, cangle, prob_max
        integer(ik), intent(in) :: npolmin, max_agap, max_pgap
        type(event_result_t), intent(out) :: res
        real(rk), intent(out) :: strike_all(:), dip_all(:), rake_all(:)
        real(rk), intent(out) :: faults_all(:, :), slips_all(:, :)
        logical, intent(in) :: use_omp
        real(rk) :: magap, mpgap, var_avg, qextra, qmismax
        integer(ik) :: npol, nspr, nmismax, nextra

        res = event_result_t()
        res%strike_avg = 999.0_rk
        res%dip_avg = 99.0_rk
        res%rake_avg = 999.0_rk
        res%var_est = 99.0_rk

        npol = count(p_pol /= 0)
        nspr = count(sp_amp /= 0.0_rk)

        if (npol < npolmin) then
            res%quality = 5
            return
        end if

        call get_gap(npsta, p_azi_mc(:, 1), p_the_mc(:, 1), magap, mpgap)
        if (magap > real(max_agap, rk) .or. mpgap > real(max_pgap, rk)) then
            res%quality = 4
            return
        end if

        nmismax = max(int(real(npol, rk) * 0.1_rk), 2)
        nextra = max(int(real(npol, rk) * 0.05_rk), 0)
        qmismax = real(nspr, rk) * qbadfac
        qextra = max(real(nspr, rk) * qbadfac * 0.5_rk, 2.0_rk)

        call focalamp_mc(aws, p_azi_mc, p_the_mc, sp_amp, p_pol, npsta, nmc, &
                         dang, maxout, nextra, nmismax, qextra, qmismax, &
                         selection, res%nf, strike_all, dip_all, rake_all, &
                         faults_all, slips_all, use_omp)

        res%nout2 = min(maxout, res%nf)
        res%nout1 = min(maxout, res%nf)
        if (res%nf == 0) then
            res%quality = 5
            return
        end if

        call mech_prob(res%nout2, faults_all(:, 1:res%nout2), slips_all(:, 1:res%nout2), &
                       cangle, prob_max, res%nsltn, res%strike_avg, res%dip_avg, &
                       res%rake_avg, res%prob, res%var_est)

        if (res%nsltn > 0) then
            call get_misf_amp(npsta, p_azi_mc(:, 1), p_the_mc(:, 1), sp_amp, &
                              p_pol, res%strike_avg(1), res%dip_avg(1), &
                              res%rake_avg(1), res%mfrac, res%mavg, res%stdr)
        else
            res%mfrac = 0.99_rk
            res%mavg = 0.99_rk
            res%stdr = 0.0_rk
        end if

        var_avg = (res%var_est(1, 1) + res%var_est(2, 1)) / 2.0_rk
        call rate_quality(res%prob(1), var_avg, res%mfrac, res%stdr, res%mavg, &
                          .true., res%quality)
        res%success = 1
    end subroutine run_event_amp

    !> Batch polarity-only pipeline.
    !> Inputs use a CSR-style flat layout:
    !>   offsets(nevent+1) - cumulative start of each event's data in the
    !>     flat arrays (picks laid out as (npol, nmc) per event)
    !>   npol_arr(nevent), nmc_arr(nevent)
    !> Output buffers are sized (maxout, nevent) and filled per event.
    subroutine run_batch(nevent, npol_arr, nmc_arr, offsets, azi_flat, the_flat, &
                         pol_flat, qual_flat, dang, maxout, badfrac, cangle, &
                         prob_max, npolmin, max_agap, max_pgap, selection, &
                         results, strike_all, dip_all, rake_all, faults_all, &
                         slips_all)
        integer(ik), intent(in) :: nevent
        integer(ik), intent(in) :: npol_arr(:), nmc_arr(:), offsets(:)
        real(rk), intent(in) :: azi_flat(:), the_flat(:)
        integer(ik), intent(in) :: pol_flat(:), qual_flat(:)
        real(rk), intent(in) :: dang, badfrac, cangle, prob_max
        integer(ik), intent(in) :: maxout, npolmin, max_agap, max_pgap, selection
        type(event_result_t), intent(out) :: results(:)
        real(rk), intent(out) :: strike_all(:, :), dip_all(:, :), rake_all(:, :)
        real(rk), intent(out) :: faults_all(:, :, :), slips_all(:, :, :)
        integer(ik) :: ie, npol, nmc, base, nthreads
        logical :: use_omp

        nthreads = hash_get_max_threads()
        ! Event-level parallelism when there are enough events; otherwise
        ! let each event parallelize internally over rotations.
        use_omp = nevent < 2 * nthreads

        if (use_omp) then
            do ie = 1, nevent
                npol = npol_arr(ie)
                nmc = nmc_arr(ie)
                base = offsets(ie)
                call run_event(fws_serial, &
                               reshape(azi_flat(base + 1:base + npol * nmc), [npol, nmc]), &
                               reshape(the_flat(base + 1:base + npol * nmc), [npol, nmc]), &
                               pol_flat(base + 1:base + npol), &
                               qual_flat(base + 1:base + npol), npol, nmc, dang, &
                               maxout, badfrac, cangle, prob_max, npolmin, &
                               max_agap, max_pgap, selection, results(ie), &
                               strike_all(:, ie), dip_all(:, ie), rake_all(:, ie), &
                               faults_all(:, :, ie), slips_all(:, :, ie), .true.)
            end do
        else
            !$omp parallel do schedule(dynamic) default(none) &
            !$omp shared(nevent, npol_arr, nmc_arr, offsets, azi_flat, the_flat, &
            !$omp         pol_flat, qual_flat, dang, maxout, badfrac, cangle, &
            !$omp         prob_max, npolmin, max_agap, max_pgap, selection, &
            !$omp         results, strike_all, dip_all, rake_all, faults_all, slips_all) &
            !$omp private(ie, npol, nmc, base)
            do ie = 1, nevent
                npol = npol_arr(ie)
                nmc = nmc_arr(ie)
                base = offsets(ie)
                call run_event(fws_thread, &
                               reshape(azi_flat(base + 1:base + npol * nmc), [npol, nmc]), &
                               reshape(the_flat(base + 1:base + npol * nmc), [npol, nmc]), &
                               pol_flat(base + 1:base + npol), &
                               qual_flat(base + 1:base + npol), npol, nmc, dang, &
                               maxout, badfrac, cangle, prob_max, npolmin, &
                               max_agap, max_pgap, selection, results(ie), &
                               strike_all(:, ie), dip_all(:, ie), rake_all(:, ie), &
                               faults_all(:, :, ie), slips_all(:, :, ie), .false.)
            end do
            !$omp end parallel do
        end if
    end subroutine run_batch

    !> Batch polarity + S/P amplitude pipeline (same layout as run_batch).
    subroutine run_batch_amp(nevent, npsta_arr, nmc_arr, offsets, azi_flat, &
                             the_flat, pol_flat, sp_flat, dang, maxout, &
                             badfrac, qbadfac, cangle, prob_max, npolmin, &
                             max_agap, max_pgap, selection, results, strike_all, &
                             dip_all, rake_all, faults_all, slips_all)
        integer(ik), intent(in) :: nevent
        integer(ik), intent(in) :: npsta_arr(:), nmc_arr(:), offsets(:)
        real(rk), intent(in) :: azi_flat(:), the_flat(:)
        integer(ik), intent(in) :: pol_flat(:)
        real(rk), intent(in) :: sp_flat(:)
        real(rk), intent(in) :: dang, badfrac, qbadfac, cangle, prob_max
        integer(ik), intent(in) :: maxout, npolmin, max_agap, max_pgap, selection
        type(event_result_t), intent(out) :: results(:)
        real(rk), intent(out) :: strike_all(:, :), dip_all(:, :), rake_all(:, :)
        real(rk), intent(out) :: faults_all(:, :, :), slips_all(:, :, :)
        integer(ik) :: ie, npsta, nmc, base, nthreads
        logical :: use_omp

        nthreads = hash_get_max_threads()
        use_omp = nevent < 2 * nthreads

        if (use_omp) then
            do ie = 1, nevent
                npsta = npsta_arr(ie)
                nmc = nmc_arr(ie)
                base = offsets(ie)
                call run_event_amp(aws_serial, &
                                   reshape(azi_flat(base + 1:base + npsta * nmc), [npsta, nmc]), &
                                   reshape(the_flat(base + 1:base + npsta * nmc), [npsta, nmc]), &
                                   pol_flat(base + 1:base + npsta), &
                                   sp_flat(base + 1:base + npsta), npsta, nmc, dang, &
                                   maxout, badfrac, qbadfac, cangle, prob_max, &
                                   npolmin, max_agap, max_pgap, selection, results(ie), &
                                   strike_all(:, ie), dip_all(:, ie), rake_all(:, ie), &
                                   faults_all(:, :, ie), slips_all(:, :, ie), .true.)
            end do
        else
            !$omp parallel do schedule(dynamic) default(none) &
            !$omp shared(nevent, npsta_arr, nmc_arr, offsets, azi_flat, the_flat, &
            !$omp         pol_flat, sp_flat, dang, maxout, badfrac, qbadfac, &
            !$omp         cangle, prob_max, npolmin, max_agap, max_pgap, selection, &
            !$omp         results, strike_all, dip_all, rake_all, faults_all, slips_all) &
            !$omp private(ie, npsta, nmc, base)
            do ie = 1, nevent
                npsta = npsta_arr(ie)
                nmc = nmc_arr(ie)
                base = offsets(ie)
                call run_event_amp(aws_thread, &
                                   reshape(azi_flat(base + 1:base + npsta * nmc), [npsta, nmc]), &
                                   reshape(the_flat(base + 1:base + npsta * nmc), [npsta, nmc]), &
                                   pol_flat(base + 1:base + npsta), &
                                   sp_flat(base + 1:base + npsta), npsta, nmc, dang, &
                                   maxout, badfrac, qbadfac, cangle, prob_max, &
                                   npolmin, max_agap, max_pgap, selection, results(ie), &
                                   strike_all(:, ie), dip_all(:, ie), rake_all(:, ie), &
                                   faults_all(:, :, ie), slips_all(:, :, ie), .false.)
            end do
            !$omp end parallel do
        end if
    end subroutine run_batch_amp

end module hash_batch
