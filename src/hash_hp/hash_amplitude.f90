!> S/P amplitude-ratio constrained grid search (FOCALAMP_MC) and
!> misfit calculation (GET_MISF_AMP).
!>
!> Port of HASH_complete/src/subs/fmamp_subs.f. The P/S radiation
!> amplitude lookup tables are shared, read-only data built once per
!> (dang-independent) table configuration; the rotation grid comes from
!> the runtime cache.
module hash_amplitude
    use hash_kinds, only: rk, ik
    use hash_geometry, only: pi, deg_to_rad, to_cart, fpcoor, cross
    use hash_rotation, only: rotation_grid_t
    use hash_runtime, only: ensure_rotation_grid_cached, grid_cache
    implicit none

    integer(ik), parameter :: ntab = 180

    type, public :: amplitude_workspace_t
        real(rk), allocatable :: px(:, :)
        real(rk), allocatable :: py(:, :)
        real(rk), allocatable :: pz(:, :)
        real(rk), allocatable :: qmis(:)
        integer(ik), allocatable :: nmis(:)
        integer(ik), allocatable :: irotgood(:)
    end type amplitude_workspace_t

    private
    public :: focalamp_mc, get_misf_amp, amplitude_prepare_workspace

contains

    !> Ensure workspace buffers match the requested sizes.
    subroutine amplitude_prepare_workspace(ws, npol, nmc, nrot)
        type(amplitude_workspace_t), intent(inout) :: ws
        integer(ik), intent(in) :: npol, nmc, nrot
        if (allocated(ws%px)) then
            if (size(ws%px, 1) == npol .and. size(ws%px, 2) == nmc) return
            deallocate (ws%px, ws%py, ws%pz, ws%qmis, ws%nmis, ws%irotgood)
        end if
        allocate (ws%px(npol, nmc), ws%py(npol, nmc), ws%pz(npol, nmc))
        allocate (ws%qmis(nrot), ws%nmis(nrot), ws%irotgood(nrot))
    end subroutine amplitude_prepare_workspace

    !> Grid search using both P polarities and S/P amplitude ratios.
    !>
    !> Inputs:
    !>   p_azi_mc(npsta,nmc), p_the_mc(npsta,nmc) - geometry per trial
    !>   sp_amp(npsta)  - log10(S/P) amplitude ratios (0 = no data)
    !>   p_pol(npsta)   - polarities (+1/-1, 0 = no polarity)
    !>   npsta, nmc, dang, maxout
    !>   nextra, ntotal  - polarity misfit criteria
    !>   qextra, qtotal  - amplitude misfit criteria
    !>   selection       - 0 = deterministic, 1 = HASH random
    !> Outputs: nf, strike/dip/rake, faults(3,nf), slips(3,nf)
    subroutine focalamp_mc(ws, p_azi_mc, p_the_mc, sp_amp, p_pol, npsta, &
                           nmc, dang, maxout, nextra, ntotal, qextra, &
                           qtotal, selection, nf, strike, dip, rake, &
                           faults, slips, use_omp)
        type(amplitude_workspace_t), intent(inout) :: ws
        real(rk), intent(in) :: p_azi_mc(:, :), p_the_mc(:, :)
        real(rk), intent(in) :: sp_amp(:)
        integer(ik), intent(in) :: p_pol(:)
        integer(ik), intent(in) :: npsta, nmc, maxout, nextra, ntotal
        real(rk), intent(in) :: dang, qextra, qtotal
        integer(ik), intent(in) :: selection
        integer(ik), intent(out) :: nf
        real(rk), intent(out) :: strike(:), dip(:), rake(:)
        real(rk), intent(out) :: faults(:, :), slips(:, :)
        logical, intent(in), optional :: use_omp
        logical :: parallel
        real(rk), allocatable :: amptable(:, :, :), phitable(:, :), thetable(:)
        integer(ik) :: nrot, i, m, irot, ista, imaxout
        real(rk) :: nmis0min, qmis0min
        real(rk) :: nmis0max, qmis0max
        integer(ik) :: nfault, nadd, nsel
        real(rk) :: pb1, pb3, prod, p_proj1, p_proj2, p_proj3, plen
        real(rk) :: pp_b1, pp_b2, theta, phi_ang, p_amp, s_amp, sp_ratio
        real(rk) :: faultnorm(3), slip(3), s1, d1, r1
        integer(ik) :: ipl, idx
        real(rk) :: astep, bbb3, bbb1, s1c, s2c
        integer(ik), allocatable :: good_irot(:)
        integer(ik) :: jran
        real(rk) :: fran
        integer(ik), parameter :: lcg_im = 120050, lcg_ia = 2311, lcg_ic = 25367

        parallel = .true.
        if (present(use_omp)) parallel = use_omp

        call ensure_rotation_grid_cached(dang)
        nrot = grid_cache%nrot
        call amplitude_prepare_workspace(ws, npsta, nmc, nrot)
        allocate (good_irot(nrot))

        ! P/S radiation amplitude lookup tables (built once per call;
        ! small: 180*360 reals).
        astep = 1.0_rk / real(ntab, rk)
        allocate (amptable(2, ntab, 2 * ntab), phitable(2 * ntab + 1, 2 * ntab + 1), &
                  thetable(2 * ntab + 1))
        do i = 1, 2 * ntab + 1
            bbb3 = -1.0_rk + real(i - 1, rk) * astep
            thetable(i) = acos(max(-1.0_rk, min(1.0_rk, bbb3)))
            do m = 1, 2 * ntab + 1
                bbb1 = -1.0_rk + real(m - 1, rk) * astep
                phi_ang = atan2(bbb3, bbb1)
                if (phi_ang < 0.0_rk) phi_ang = phi_ang + 2.0_rk * pi
                phitable(i, m) = phi_ang
            end do
        end do
        do i = 1, 2 * ntab
            phi_ang = real(i - 1, rk) * pi * astep
            do m = 1, ntab
                theta = real(m - 1, rk) * pi * astep
                amptable(1, m, i) = abs(sin(2.0_rk * theta) * cos(phi_ang))
                s1c = cos(2.0_rk * theta) * cos(phi_ang)
                s2c = -cos(theta) * sin(phi_ang)
                amptable(2, m, i) = sqrt(s1c * s1c + s2c * s2c)
            end do
        end do

        ! Hoist spherical -> Cartesian conversions.
        do m = 1, nmc
            do i = 1, npsta
                call to_cart(p_the_mc(i, m), p_azi_mc(i, m), 1.0_rk, &
                             ws%px(i, m), ws%py(i, m), ws%pz(i, m))
            end do
        end do

        ws%irotgood = 0

        do m = 1, nmc
            nmis0min = 1.0e5_rk
            qmis0min = 1.0e5_rk

            call count_misfits_trial_amp(ws, m, npsta, nrot, p_pol, sp_amp, &
                                         parallel, amptable, phitable, thetable, &
                                         astep, nmis0min, qmis0min)

            nmis0max = real(ntotal, rk)
            if (nmis0max < nmis0min + real(nextra, rk)) nmis0max = nmis0min + real(nextra, rk)
            qmis0max = qtotal
            if (qmis0max < qmis0min + qextra) qmis0max = qmis0min + qextra

            ! Mark rotations meeting both criteria; if none, loosen the
            ! amplitude criterion to the minimum and retry (original logic).
            do
                nadd = 0
                do irot = 1, nrot
                    if (real(ws%nmis(irot), rk) <= nmis0max .and. ws%qmis(irot) <= qmis0max) then
                        ws%irotgood(irot) = 1
                        nadd = nadd + 1
                    end if
                end do
                if (nadd /= 0) exit
                qmis0min = 1.0e5_rk
                do irot = 1, nrot
                    if (real(ws%nmis(irot), rk) <= nmis0max .and. ws%qmis(irot) < qmis0min) then
                        qmis0min = ws%qmis(irot)
                    end if
                end do
                qmis0max = qtotal
                if (qmis0max < qmis0min + qextra) qmis0max = qmis0min + qextra
            end do
        end do

        ! Collect acceptable rotations in grid order.
        nfault = 0
        do irot = 1, nrot
            if (ws%irotgood(irot) > 0) then
                nfault = nfault + 1
                good_irot(nfault) = irot
            end if
        end do

        ! Select output solutions (deterministic or HASH random).
        nf = 0
        imaxout = maxout
        if (nfault <= imaxout) then
            do i = 1, nfault
                nf = nf + 1
                call emit_mechanism(good_irot(i), faultnorm, slip, s1, d1, r1, &
                                    strike, dip, rake, faults, slips, nf)
            end do
        else if (selection == 0) then
            do i = 1, imaxout
                nf = nf + 1
                call emit_mechanism(good_irot(i), faultnorm, slip, s1, d1, r1, &
                                    strike, dip, rake, faults, slips, nf)
            end do
        else
            jran = 314159
            nsel = 0
            do while (nsel < imaxout)
                jran = mod(jran * lcg_ia + lcg_ic, lcg_im)
                fran = real(jran, rk) / real(lcg_im, rk)
                idx = nint(fran * real(nfault, rk) + 0.5_rk)
                if (idx < 1) idx = 1
                if (idx > nfault) idx = nfault
                if (good_irot(idx) <= 0) cycle
                nsel = nsel + 1
                nf = nf + 1
                call emit_mechanism(good_irot(idx), faultnorm, slip, s1, d1, r1, &
                                    strike, dip, rake, faults, slips, nf)
                good_irot(idx) = -1
            end do
        end if

    end subroutine focalamp_mc

    !> Count polarity + amplitude misfits over all rotations for one trial.
    subroutine count_misfits_trial_amp(ws, m, npsta, nrot, p_pol, sp_amp, &
                                       parallel, amptable, phitable, thetable, &
                                       astep, nmis0min, qmis0min)
        type(amplitude_workspace_t), intent(inout) :: ws
        integer(ik), intent(in) :: m, npsta, nrot
        integer(ik), intent(in) :: p_pol(:)
        real(rk), intent(in) :: sp_amp(:)
        logical, intent(in) :: parallel
        real(rk), intent(in) :: amptable(:, :, :), phitable(:, :), thetable(:)
        real(rk), intent(in) :: astep
        real(rk), intent(inout) :: nmis0min, qmis0min
        integer(ik) :: irot, ista
        integer(ik) :: nm, ipol, ix, iy
        real(rk) :: qm, pb1, pb3, pr, pj1, pj2, pj3, pln, p1, p2
        real(rk) :: tht, pha, pam, sam, sratio

        if (parallel) then
            !$omp parallel default(none) private(ista) &
            !$omp shared(grid_cache, ws, p_pol, sp_amp, npsta, nrot, m, amptable, &
            !$omp         phitable, thetable, astep, nmis0min, qmis0min)
            block
                integer(ik) :: nm, ipol, ix, iy
                real(rk) :: qm, pb1, pb3, pr, pj1, pj2, pj3, pln, p1, p2
                real(rk) :: tht, pha, pam, sam, sratio
                real(rk) :: t_nmis0min, t_qmis0min
                t_nmis0min = 1.0e5_rk
                t_qmis0min = 1.0e5_rk
                !$omp do schedule(static)
                do irot = 1, nrot
                    nm = 0
                    qm = 0.0_rk
                    !$omp simd reduction(+:nm, qm)
                    do ista = 1, npsta
                        pb1 = grid_cache%b1(1, irot) * ws%px(ista, m) &
                            + grid_cache%b1(2, irot) * ws%py(ista, m) &
                            + grid_cache%b1(3, irot) * ws%pz(ista, m)
                        pb3 = grid_cache%b3(1, irot) * ws%px(ista, m) &
                            + grid_cache%b3(2, irot) * ws%py(ista, m) &
                            + grid_cache%b3(3, irot) * ws%pz(ista, m)
                        if (sp_amp(ista) /= 0.0_rk) then
                            pj1 = ws%px(ista, m) - pb3 * grid_cache%b3(1, irot)
                            pj2 = ws%py(ista, m) - pb3 * grid_cache%b3(2, irot)
                            pj3 = ws%pz(ista, m) - pb3 * grid_cache%b3(3, irot)
                            pln = sqrt(pj1 * pj1 + pj2 * pj2 + pj3 * pj3)
                            if (pln > 0.0_rk) then
                                pj1 = pj1 / pln
                                pj2 = pj2 / pln
                                pj3 = pj3 / pln
                            end if
                            p1 = grid_cache%b1(1, irot) * pj1 + grid_cache%b1(2, irot) * pj2 &
                                 + grid_cache%b1(3, irot) * pj3
                            p2 = grid_cache%b2(1, irot) * pj1 + grid_cache%b2(2, irot) * pj2 &
                                 + grid_cache%b2(3, irot) * pj3
                            ix = nint((pb3 + 1.0_rk) / astep) + 1
                            tht = thetable(ix)
                            ix = nint((p2 + 1.0_rk) / astep) + 1
                            iy = nint((p1 + 1.0_rk) / astep) + 1
                            pha = phitable(ix, iy)
                            ix = nint(pha / (pi * astep)) + 1
                            if (ix > 2 * ntab) ix = 1
                            iy = nint(tht / (pi * astep)) + 1
                            if (iy > ntab) iy = 1
                            pam = amptable(1, iy, ix)
                            sam = amptable(2, iy, ix)
                            if (pam == 0.0_rk) then
                                sratio = 4.0_rk
                            else if (sam == 0.0_rk) then
                                sratio = -2.0_rk
                            else
                                sratio = log10(4.9_rk * sam / pam)
                            end if
                            qm = qm + abs(sp_amp(ista) - sratio)
                        end if
                        if (p_pol(ista) /= 0) then
                            pr = pb1 * pb3
                            ipol = -1
                            if (pr > 0.0_rk) ipol = 1
                            if (ipol /= p_pol(ista)) nm = nm + 1
                        end if
                    end do
                    ws%qmis(irot) = qm
                    ws%nmis(irot) = nm
                    if (real(nm, rk) < t_nmis0min) t_nmis0min = real(nm, rk)
                    if (qm < t_qmis0min) t_qmis0min = qm
                end do
                !$omp end do
                !$omp critical(hash_amp_min)
                if (t_nmis0min < nmis0min) nmis0min = t_nmis0min
                if (t_qmis0min < qmis0min) qmis0min = t_qmis0min
                !$omp end critical(hash_amp_min)
            end block
            !$omp end parallel
        else
            do irot = 1, nrot
                nm = 0
                qm = 0.0_rk
                do ista = 1, npsta
                    pb1 = grid_cache%b1(1, irot) * ws%px(ista, m) &
                        + grid_cache%b1(2, irot) * ws%py(ista, m) &
                        + grid_cache%b1(3, irot) * ws%pz(ista, m)
                    pb3 = grid_cache%b3(1, irot) * ws%px(ista, m) &
                        + grid_cache%b3(2, irot) * ws%py(ista, m) &
                        + grid_cache%b3(3, irot) * ws%pz(ista, m)
                    if (sp_amp(ista) /= 0.0_rk) then
                        pj1 = ws%px(ista, m) - pb3 * grid_cache%b3(1, irot)
                        pj2 = ws%py(ista, m) - pb3 * grid_cache%b3(2, irot)
                        pj3 = ws%pz(ista, m) - pb3 * grid_cache%b3(3, irot)
                        pln = sqrt(pj1 * pj1 + pj2 * pj2 + pj3 * pj3)
                        if (pln > 0.0_rk) then
                            pj1 = pj1 / pln
                            pj2 = pj2 / pln
                            pj3 = pj3 / pln
                        end if
                        p1 = grid_cache%b1(1, irot) * pj1 + grid_cache%b1(2, irot) * pj2 &
                             + grid_cache%b1(3, irot) * pj3
                        p2 = grid_cache%b2(1, irot) * pj1 + grid_cache%b2(2, irot) * pj2 &
                             + grid_cache%b2(3, irot) * pj3
                        ix = nint((pb3 + 1.0_rk) / astep) + 1
                        tht = thetable(ix)
                        ix = nint((p2 + 1.0_rk) / astep) + 1
                        iy = nint((p1 + 1.0_rk) / astep) + 1
                        pha = phitable(ix, iy)
                        ix = nint(pha / (pi * astep)) + 1
                        if (ix > 2 * ntab) ix = 1
                        iy = nint(tht / (pi * astep)) + 1
                        if (iy > ntab) iy = 1
                        pam = amptable(1, iy, ix)
                        sam = amptable(2, iy, ix)
                        if (pam == 0.0_rk) then
                            sratio = 4.0_rk
                        else if (sam == 0.0_rk) then
                            sratio = -2.0_rk
                        else
                            sratio = log10(4.9_rk * sam / pam)
                        end if
                        qm = qm + abs(sp_amp(ista) - sratio)
                    end if
                    if (p_pol(ista) /= 0) then
                        pr = pb1 * pb3
                        ipol = -1
                        if (pr > 0.0_rk) ipol = 1
                        if (ipol /= p_pol(ista)) nm = nm + 1
                    end if
                end do
                ws%qmis(irot) = qm
                ws%nmis(irot) = nm
                if (real(nm, rk) < nmis0min) nmis0min = real(nm, rk)
                if (qm < qmis0min) qmis0min = qm
            end do
        end if
    end subroutine count_misfits_trial_amp

    !> Write one mechanism into the output arrays.
    subroutine emit_mechanism(irot, faultnorm, slip, s1, d1, r1, &
                              strike, dip, rake, faults, slips, idx)
        integer(ik), intent(in) :: irot, idx
        real(rk), intent(inout) :: faultnorm(3), slip(3)
        real(rk), intent(out) :: s1, d1, r1
        real(rk), intent(inout) :: strike(:), dip(:), rake(:)
        real(rk), intent(inout) :: faults(:, :), slips(:, :)
        integer(ik) :: m
        faultnorm(1) = grid_cache%b3(1, irot)
        faultnorm(2) = grid_cache%b3(2, irot)
        faultnorm(3) = grid_cache%b3(3, irot)
        slip(1) = grid_cache%b1(1, irot)
        slip(2) = grid_cache%b1(2, irot)
        slip(3) = grid_cache%b1(3, irot)
        do m = 1, 3
            faults(m, idx) = faultnorm(m)
            slips(m, idx) = slip(m)
        end do
        call fpcoor(s1, d1, r1, faultnorm, slip, 2)
        strike(idx) = s1
        dip(idx) = d1
        rake(idx) = r1
    end subroutine emit_mechanism

    !> Misfit fractions for a given mechanism with S/P ratios.
    !> Outputs: mfrac (polarity misfit fraction), mavg (mean |SP misfit|),
    !>          stdr (station distribution ratio).
    subroutine get_misf_amp(npol, p_azi, p_the, sp_ratio, p_pol, str_avg, &
                            dip_avg, rak_avg, mfrac, mavg, stdr)
        integer(ik), intent(in) :: npol
        real(rk), intent(in) :: p_azi(:), p_the(:)
        real(rk), intent(in) :: sp_ratio(:)
        integer(ik), intent(in) :: p_pol(:)
        real(rk), intent(in) :: str_avg, dip_avg, rak_avg
        real(rk), intent(out) :: mfrac, mavg, stdr
        real(rk) :: M(3, 3), a(3), b(3)
        real(rk) :: strike, dip, rake, qcount, scount, acount
        real(rk) :: azi, toff, wt, wo, pol
        real(rk) :: bb1(3), bb2(3), bb3(3)
        real(rk) :: sd, dd, rdg
        real(rk) :: p_a1, p_a2, p_a3, p_b1, p_b3
        real(rk) :: p_proj1, p_proj2, p_proj3, plen, pp_b1, pp_b2
        real(rk) :: phi_ang, theta_ang, p_amp, dot_val
        real(rk) :: s1c, s2c, s_amp, sp_rat
        integer(ik) :: k, in, jn

        strike = str_avg * deg_to_rad
        dip = dip_avg * deg_to_rad
        rake = rak_avg * deg_to_rad

        M(1, 1) = -sin(dip) * cos(rake) * sin(2 * strike) &
                  - sin(2 * dip) * sin(rake) * sin(strike) * sin(strike)
        M(2, 2) = sin(dip) * cos(rake) * sin(2 * strike) &
                  - sin(2 * dip) * sin(rake) * cos(strike) * cos(strike)
        M(3, 3) = sin(2 * dip) * sin(rake)
        M(1, 2) = sin(dip) * cos(rake) * cos(2 * strike) &
                  + 0.5_rk * sin(2 * dip) * sin(rake) * sin(2 * strike)
        M(2, 1) = M(1, 2)
        M(1, 3) = -cos(dip) * cos(rake) * cos(strike) &
                  - cos(2 * dip) * sin(rake) * sin(strike)
        M(3, 1) = M(1, 3)
        M(2, 3) = -cos(dip) * cos(rake) * sin(strike) &
                  + cos(2 * dip) * sin(rake) * cos(strike)
        M(3, 2) = M(2, 3)

        ! Moment tensor uses radians; the fault vectors use the original
        ! degree values (the original HASH misread radians as degrees here).
        sd = str_avg
        dd = dip_avg
        rdg = rak_avg
        call fpcoor(sd, dd, rdg, bb3, bb1, 1)
        call cross(bb3, bb1, bb2)

        mfrac = 0.0_rk
        qcount = 0.0_rk
        stdr = 0.0_rk
        scount = 0.0_rk
        mavg = 0.0_rk
        acount = 0.0_rk

        do k = 1, npol
            call to_cart(p_the(k), p_azi(k), 1.0_rk, p_a1, p_a2, p_a3)
            p_b1 = bb1(1) * p_a1 + bb1(2) * p_a2 + bb1(3) * p_a3
            p_b3 = bb3(1) * p_a1 + bb3(2) * p_a2 + bb3(3) * p_a3
            p_proj1 = p_a1 - p_b3 * bb3(1)
            p_proj2 = p_a2 - p_b3 * bb3(2)
            p_proj3 = p_a3 - p_b3 * bb3(3)
            plen = sqrt(p_proj1 * p_proj1 + p_proj2 * p_proj2 + p_proj3 * p_proj3)
            if (plen > 0.0_rk) then
                p_proj1 = p_proj1 / plen
                p_proj2 = p_proj2 / plen
                p_proj3 = p_proj3 / plen
            end if
            pp_b1 = bb1(1) * p_proj1 + bb1(2) * p_proj2 + bb1(3) * p_proj3
            pp_b2 = bb2(1) * p_proj1 + bb2(2) * p_proj2 + bb2(3) * p_proj3
            phi_ang = atan2(pp_b2, pp_b1)
            theta_ang = acos(max(-1.0_rk, min(1.0_rk, p_b3)))
            p_amp = abs(sin(2 * theta_ang) * cos(phi_ang))
            wt = sqrt(p_amp)

            if (p_pol(k) /= 0) then
                azi = p_azi(k) * deg_to_rad
                toff = p_the(k) * deg_to_rad
                a(1) = sin(toff) * cos(azi)
                a(2) = sin(toff) * sin(azi)
                a(3) = -cos(toff)
                do in = 1, 3
                    b(in) = 0.0_rk
                    do jn = 1, 3
                        b(in) = b(in) + M(in, jn) * a(jn)
                    end do
                end do
                dot_val = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
                if (dot_val < 0.0_rk) then
                    pol = -1.0_rk
                else
                    pol = 1.0_rk
                end if
                if (pol * real(p_pol(k), rk) < 0.0_rk) mfrac = mfrac + wt
                qcount = qcount + wt
                stdr = stdr + wt
                scount = scount + 1.0_rk
            end if

            if (sp_ratio(k) /= 0.0_rk) then
                s1c = cos(2 * theta_ang) * cos(phi_ang)
                s2c = -cos(theta_ang) * sin(phi_ang)
                s_amp = sqrt(s1c * s1c + s2c * s2c)
                if (p_amp == 0.0_rk) then
                    sp_rat = 4.0_rk
                else
                    sp_rat = log10(4.9_rk * s_amp / p_amp)
                end if
                mavg = mavg + abs(sp_ratio(k) - sp_rat)
                acount = acount + 1.0_rk
                stdr = stdr + wt
                scount = scount + 1.0_rk
            end if
        end do

        if (qcount > 0.0_rk) then
            mfrac = mfrac / qcount
        else
            mfrac = 0.0_rk
        end if
        if (acount > 0.0_rk) then
            mavg = mavg / acount
        else
            mavg = 0.0_rk
        end if
        if (scount > 0.0_rk) then
            stdr = stdr / scount
        else
            stdr = 0.0_rk
        end if

    end subroutine get_misf_amp

end module hash_amplitude
