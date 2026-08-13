!> Core focal-mechanism grid search (FOCALMC).
!>
!> Faithful port of HASH_complete/src/subs/fmech_subs.f FOCALMC, with:
!>   - all azimuth/takeoff -> Cartesian conversions hoisted out of the
!>     rotation loop (one pass over (npol, nmc), then only multiply/add)
!>   - optional OpenMP parallelization over rotations (independent per
!>     rotation). Batch mode passes use_omp=.false. when it parallelizes
!>     over events to avoid nested parallelism.
!>   - deterministic solution selection by default (first maxout in grid
!>     order); optional HASH-style random selection
!>   - no I/O, no allocations inside the hot loops
module hash_focalmc
    use hash_kinds, only: rk, ik
    use hash_geometry, only: to_cart, fpcoor
    use hash_rotation, only: rotation_grid_t
    use hash_runtime, only: ensure_rotation_grid_cached, grid_cache, &
                            hash_get_max_threads
#ifdef _OPENMP
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
#endif
    implicit none

    type, public :: focalmc_workspace_t
        real(rk), allocatable :: px(:, :)
        real(rk), allocatable :: py(:, :)
        real(rk), allocatable :: pz(:, :)
        integer(ik), allocatable :: fit0(:)
        integer(ik), allocatable :: fit(:)
        integer(ik), allocatable :: irotgood(:)
    end type focalmc_workspace_t

    private
    public :: focalmc, focalmc_prepare_workspace

contains

    !> Ensure workspace buffers match the requested sizes (reused across calls).
    subroutine focalmc_prepare_workspace(ws, npol, nmc, nrot)
        type(focalmc_workspace_t), intent(inout) :: ws
        integer(ik), intent(in) :: npol, nmc, nrot
        if (allocated(ws%px)) then
            ! nrot varies with dang, so the fit arrays must be checked too.
            if (size(ws%px, 1) == npol .and. size(ws%px, 2) == nmc .and. &
                size(ws%fit) == nrot) return
            deallocate (ws%px, ws%py, ws%pz, ws%fit0, ws%fit, ws%irotgood)
        end if
        allocate (ws%px(npol, nmc), ws%py(npol, nmc), ws%pz(npol, nmc))
        allocate (ws%fit0(nrot), ws%fit(nrot), ws%irotgood(nrot))
    end subroutine focalmc_prepare_workspace

    !> Grid search for acceptable focal mechanisms.
    !>
    !> Inputs:
    !>   p_azi_mc(npol,nmc) - azimuths (deg E of N) per observation and trial
    !>   p_the_mc(npol,nmc) - takeoff angles (deg, from vertical)
    !>   p_pol(npol)        - first motion, +1=up, -1=down
    !>   p_qual(npol)       - quality, 0=impulsive, 1=emergent
    !>   npol               - number of first motions
    !>   nmc                - number of trials
    !>   dang               - grid angle spacing (degrees)
    !>   maxout             - maximum number of fault planes to return
    !>   nextra             - additional misfits allowed above minimum
    !>   ntotal             - total number of allowed misfits
    !>   selection          - 0 = deterministic (first maxout), 1 = HASH random
    !>   use_omp            - parallelize over rotations (default .true.)
    !> Outputs:
    !>   nf                 - number of fault planes found
    !>   strike/dip/rake    - mechanism angles for each returned plane
    !>   faults(3,nf)       - fault normal vectors
    !>   slips(3,nf)        - slip vectors
    subroutine focalmc(ws, p_azi_mc, p_the_mc, p_pol, p_qual, npol, nmc, &
                       dang, maxout, nextra, ntotal, selection, nf, &
                       strike, dip, rake, faults, slips, use_omp)
        type(focalmc_workspace_t), intent(inout) :: ws
        real(rk), intent(in) :: p_azi_mc(:, :), p_the_mc(:, :)
        integer(ik), intent(in) :: p_pol(:), p_qual(:)
        integer(ik), intent(in) :: npol, nmc, maxout, nextra, ntotal, selection
        real(rk), intent(in) :: dang
        integer(ik), intent(out) :: nf
        real(rk), intent(out) :: strike(:), dip(:), rake(:)
        real(rk), intent(out) :: faults(:, :), slips(:, :)
        logical, intent(in), optional :: use_omp
        logical :: parallel
        integer(ik) :: nrot, i, m, irot, imaxout
        integer(ik) :: nmiss0min, nmissmin
        integer(ik) :: nmiss0max, nmissmax
        integer(ik) :: nfault
        real(rk) :: faultnorm(3), slip(3), s1, d1, r1
        integer(ik), allocatable :: nmiss01min(:)
        integer(ik), allocatable :: good_irot(:)
        integer(ik) :: jran, nsel
        real(rk) :: fran
        integer(ik), parameter :: lcg_im = 120050, lcg_ia = 2311, lcg_ic = 25367

        parallel = .true.
        if (present(use_omp)) parallel = use_omp
        ! A single OpenMP thread gets no benefit from a parallel region;
        ! run the serial path to avoid the fork/join overhead.
        if (parallel .and. hash_get_max_threads() <= 1) parallel = .false.

        call ensure_rotation_grid_cached(dang)

        nrot = grid_cache%nrot
        call focalmc_prepare_workspace(ws, npol, nmc, nrot)
        allocate (nmiss01min(0:npol), good_irot(nrot))

        ! Hoist all spherical -> Cartesian conversions out of the rotation loop.
        do m = 1, nmc
            do i = 1, npol
                call to_cart(p_the_mc(i, m), p_azi_mc(i, m), 1.0_rk, &
                             ws%px(i, m), ws%py(i, m), ws%pz(i, m))
            end do
        end do

        ws%irotgood = 0

        do m = 1, nmc
            nmiss0min = 999
            nmissmin = 999
            nmiss01min(:) = 999

            call count_misfits_trial(ws, m, npol, nrot, p_pol, p_qual, parallel, &
                                     nmiss0min, nmissmin, nmiss01min)

            ! Choose fit criteria (exact FOCALMC logic).
            if (nmiss0min == 0) then
                nmiss0max = ntotal
                nmissmax = ntotal
            else
                nmiss0max = ntotal
                nmissmax = npol
            end if
            if (nmiss0max < nmiss0min + nextra) nmiss0max = nmiss0min + nextra
            if (nmissmax < nmiss01min(nmiss0min) + nextra) &
                nmissmax = nmiss01min(nmiss0min) + nextra

            if (parallel) then
                !$omp parallel do schedule(static) default(none) &
                !$omp shared(grid_cache, ws, nrot, nmiss0max, nmissmax)
                do irot = 1, nrot
                    if (ws%fit0(irot) <= nmiss0max .and. ws%fit(irot) <= nmissmax) then
                        ws%irotgood(irot) = 1
                    end if
                end do
                !$omp end parallel do
            else
                do irot = 1, nrot
                    if (ws%fit0(irot) <= nmiss0max .and. ws%fit(irot) <= nmissmax) then
                        ws%irotgood(irot) = 1
                    end if
                end do
            end if
        end do

        ! Collect the acceptable rotations in grid order.
        nfault = 0
        do irot = 1, nrot
            if (ws%irotgood(irot) > 0) then
                nfault = nfault + 1
                good_irot(nfault) = irot
            end if
        end do

        ! Select output solutions.
        nf = 0
        imaxout = maxout
        if (nfault <= imaxout) then
            do i = 1, nfault
                nf = nf + 1
                call emit_mechanism(good_irot(i), faultnorm, slip, s1, d1, r1, &
                                    strike, dip, rake, faults, slips, nf)
            end do
        else if (selection == 0) then
            ! Deterministic: first maxout in grid order (matches CNCHASH Numba).
            do i = 1, imaxout
                nf = nf + 1
                call emit_mechanism(good_irot(i), faultnorm, slip, s1, d1, r1, &
                                    strike, dip, rake, faults, slips, nf)
            end do
        else
            ! HASH-style random selection (original FOCALMC rand(0) behavior,
            ! reproduced with a deterministic LCG so it is still reproducible).
            jran = 314159
            nsel = 0
            do while (nsel < imaxout)
                jran = mod(jran * lcg_ia + lcg_ic, lcg_im)
                fran = real(jran, rk) / real(lcg_im, rk)
                irot = nint(fran * real(nfault, rk) + 0.5_rk)
                if (irot < 1) irot = 1
                if (irot > nfault) irot = nfault
                if (good_irot(irot) <= 0) cycle
                nsel = nsel + 1
                nf = nf + 1
                call emit_mechanism(good_irot(irot), faultnorm, slip, s1, d1, r1, &
                                    strike, dip, rake, faults, slips, nf)
                good_irot(irot) = -1
            end do
        end if

    end subroutine focalmc

    !> Count polarity misfits over all rotations for one MC trial.
    subroutine count_misfits_trial(ws, m, npol, nrot, p_pol, p_qual, parallel, &
                                   nmiss0min, nmissmin, nmiss01min)
        type(focalmc_workspace_t), intent(inout) :: ws
        integer(ik), intent(in) :: m, npol, nrot
        integer(ik), intent(in) :: p_pol(:), p_qual(:)
        logical, intent(in) :: parallel
        integer(ik), intent(inout) :: nmiss0min, nmissmin
        integer(ik), intent(inout) :: nmiss01min(0:)
        integer(ik) :: irot, ista
        integer(ik) :: nm0, nm1, ipp, iq
        real(rk) :: pb1, pb3, pr, pa1, pa2, pa3

        if (parallel) then
            ! Transposed parallel path: every thread owns a contiguous
            ! rotation block and runs the same station-outer /
            ! rotation-inner loop as the serial path, so the inner loop
            ! is long, contiguous in memory and branch-free.
            ws%fit(1:nrot) = 0
            ws%fit0(1:nrot) = 0
            !$omp parallel default(none) &
            !$omp shared(grid_cache, ws, p_pol, p_qual, npol, nrot, m, &
            !$omp         nmiss01min, nmiss0min, nmissmin)
            block
                integer(ik) :: tid, nthreads, istart, istop, nchunk, il, ipp, iq
                integer(ik) :: nm0, nm1
                integer(ik) :: t_nmiss0min, t_nmissmin
                integer(ik), allocatable :: t_nmiss01min(:)
                integer(ik), allocatable :: t_fit(:), t_fit0(:)
                real(rk) :: pb1, pb3, pr, pa1, pa2, pa3
                allocate (t_nmiss01min(0:npol))
                t_nmiss0min = 999
                t_nmissmin = 999
                t_nmiss01min(:) = 999
                tid = omp_get_thread_num()
                nthreads = omp_get_num_threads()
                ! Static partition of the rotation range.
                istart = tid * (nrot / nthreads) + 1
                istop = (tid + 1) * (nrot / nthreads)
                if (tid == nthreads - 1) istop = nrot
                nchunk = istop - istart + 1
                allocate (t_fit(nchunk), t_fit0(nchunk))
                t_fit = 0
                t_fit0 = 0
                do ista = 1, npol
                    ipp = p_pol(ista)
                    iq = p_qual(ista)
                    pa1 = ws%px(ista, m)
                    pa2 = ws%py(ista, m)
                    pa3 = ws%pz(ista, m)
                    !$omp simd
                    do il = 1, nchunk
                        pb1 = grid_cache%b1(istart + il - 1, 1) * pa1 &
                            + grid_cache%b1(istart + il - 1, 2) * pa2 &
                            + grid_cache%b1(istart + il - 1, 3) * pa3
                        pb3 = grid_cache%b3(istart + il - 1, 1) * pa1 &
                            + grid_cache%b3(istart + il - 1, 2) * pa2 &
                            + grid_cache%b3(istart + il - 1, 3) * pa3
                        pr = pb1 * pb3
                        t_fit(il) = t_fit(il) &
                                    + merge(1_ik, 0_ik, (pr > 0.0_rk) .neqv. (ipp == 1))
                        t_fit0(il) = t_fit0(il) &
                                     + merge(1_ik, 0_ik, ((pr > 0.0_rk) .neqv. (ipp == 1)) &
                                               .and. (iq == 0))
                    end do
                end do
                do il = 1, nchunk
                    ws%fit(istart + il - 1) = t_fit(il)
                    ws%fit0(istart + il - 1) = t_fit0(il)
                    nm0 = t_fit0(il)
                    nm1 = t_fit(il)
                    if (nm0 < t_nmiss0min) t_nmiss0min = nm0
                    if (nm1 < t_nmissmin) t_nmissmin = nm1
                    if (nm1 < t_nmiss01min(nm0)) t_nmiss01min(nm0) = nm1
                end do
                !$omp critical(hash_focalmc_min)
                if (t_nmiss0min < nmiss0min) nmiss0min = t_nmiss0min
                if (t_nmissmin < nmissmin) nmissmin = t_nmissmin
                do ista = 0, npol
                    if (t_nmiss01min(ista) < nmiss01min(ista)) &
                        nmiss01min(ista) = t_nmiss01min(ista)
                end do
                !$omp end critical(hash_focalmc_min)
            end block
            !$omp end parallel
        else
            ! Transposed serial path: stations outer, rotations inner.
            ! The rotation loop is long (31032), contiguous in memory and
            ! branch-free, so it vectorizes well; misfit counts accumulate
            ! directly into ws%fit/ws%fit0 across stations.
            ws%fit(1:nrot) = 0
            ws%fit0(1:nrot) = 0
            do ista = 1, npol
                ipp = p_pol(ista)
                iq = p_qual(ista)
                pa1 = ws%px(ista, m)
                pa2 = ws%py(ista, m)
                pa3 = ws%pz(ista, m)
                !$omp simd
                do irot = 1, nrot
                    pb1 = grid_cache%b1(irot, 1) * pa1 &
                        + grid_cache%b1(irot, 2) * pa2 &
                        + grid_cache%b1(irot, 3) * pa3
                    pb3 = grid_cache%b3(irot, 1) * pa1 &
                        + grid_cache%b3(irot, 2) * pa2 &
                        + grid_cache%b3(irot, 3) * pa3
                    pr = pb1 * pb3
                    ws%fit(irot) = ws%fit(irot) &
                                   + merge(1_ik, 0_ik, (pr > 0.0_rk) .neqv. (ipp == 1))
                    ws%fit0(irot) = ws%fit0(irot) &
                                    + merge(1_ik, 0_ik, ((pr > 0.0_rk) .neqv. (ipp == 1)) &
                                              .and. (iq == 0))
                end do
            end do
            do irot = 1, nrot
                nm0 = ws%fit0(irot)
                nm1 = ws%fit(irot)
                if (nm0 < nmiss0min) nmiss0min = nm0
                if (nm1 < nmissmin) nmissmin = nm1
                if (nm1 < nmiss01min(nm0)) nmiss01min(nm0) = nm1
            end do
        end if
    end subroutine count_misfits_trial

    !> Write one mechanism (vectors + angles) into the output arrays.
    subroutine emit_mechanism(irot, faultnorm, slip, s1, d1, r1, &
                              strike, dip, rake, faults, slips, idx)
        integer(ik), intent(in) :: irot, idx
        real(rk), intent(inout) :: faultnorm(3), slip(3)
        real(rk), intent(out) :: s1, d1, r1
        real(rk), intent(inout) :: strike(:), dip(:), rake(:)
        real(rk), intent(inout) :: faults(:, :), slips(:, :)
        integer(ik) :: m
        faultnorm(1) = grid_cache%b3(irot, 1)
        faultnorm(2) = grid_cache%b3(irot, 2)
        faultnorm(3) = grid_cache%b3(irot, 3)
        slip(1) = grid_cache%b1(irot, 1)
        slip(2) = grid_cache%b1(irot, 2)
        slip(3) = grid_cache%b1(irot, 3)
        do m = 1, 3
            faults(m, idx) = faultnorm(m)
            slips(m, idx) = slip(m)
        end do
        call fpcoor(s1, d1, r1, faultnorm, slip, 2)
        strike(idx) = s1
        dip(idx) = d1
        rake(idx) = r1
    end subroutine emit_mechanism

end module hash_focalmc
