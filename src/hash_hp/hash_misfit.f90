!> Misfit and station-coverage functions (GET_GAP, GET_MISF).
!>
!> Port of HASH_complete/src/subs/pol_subs.f.
!>
!> Notes on deliberate differences from the original:
!>   - Original GET_GAP implicitly typed magap/mpgap as INTEGER (implicit
!>     typing, 'm' in [i-n]); CNCHASH returns real gaps throughout so the
!>     Fortran and Python backends behave identically. The original
!>     truncation could only affect borderline max_agap/max_pgap checks.
module hash_misfit
    use hash_kinds, only: rk, ik
    use hash_geometry, only: pi, deg_to_rad, to_cart, fpcoor, cross
    implicit none

    private
    public :: get_gap, get_misfit, sort_asc

contains

    !> Find the maximum azimuthal and takeoff-angle gaps.
    !>
    !> Inputs: npol, p_azi(npol), p_the(npol)
    !> Outputs: magap (deg), mpgap (deg)
    subroutine get_gap(npol, p_azi, p_the, magap, mpgap)
        integer(ik), intent(in) :: npol
        real(rk), intent(in) :: p_azi(:), p_the(:)
        real(rk), intent(out) :: magap, mpgap
        real(rk), allocatable :: p2_azi(:), p2_the(:)
        integer(ik) :: k

        allocate (p2_azi(npol), p2_the(npol))

        do k = 1, npol
            if (p_the(k) > 90.0_rk) then
                p2_the(k) = 180.0_rk - p_the(k)
                p2_azi(k) = p_azi(k) - 180.0_rk
                if (p2_azi(k) < 0.0_rk) p2_azi(k) = p2_azi(k) + 360.0_rk
            else
                p2_the(k) = p_the(k)
                p2_azi(k) = p_azi(k)
            end if
        end do

        call sort_asc(npol, p2_azi)
        call sort_asc(npol, p2_the)

        magap = 0.0_rk
        mpgap = 0.0_rk
        do k = 2, npol
            if (p2_azi(k) - p2_azi(k - 1) > magap) magap = p2_azi(k) - p2_azi(k - 1)
            if (p2_the(k) - p2_the(k - 1) > mpgap) mpgap = p2_the(k) - p2_the(k - 1)
        end do
        if (p2_azi(1) - p2_azi(npol) + 360.0_rk > magap) &
            magap = p2_azi(1) - p2_azi(npol) + 360.0_rk
        if (90.0_rk - p2_the(npol) > mpgap) mpgap = 90.0_rk - p2_the(npol)
        if (p2_the(1) > mpgap) mpgap = p2_the(1)

    end subroutine get_gap

    !> Weighted fraction of misfit polarities and station distribution ratio.
    !>
    !> Inputs: npol, p_azi(npol), p_the(npol), p_pol(npol), p_qual(npol),
    !>         str_avg, dip_avg, rak_avg (mechanism, degrees)
    !> Outputs: mfrac (weighted misfit fraction), stdr (distribution ratio)
    subroutine get_misfit(npol, p_azi, p_the, p_pol, p_qual, str_avg, &
                          dip_avg, rak_avg, mfrac, stdr)
        integer(ik), intent(in) :: npol
        real(rk), intent(in) :: p_azi(:), p_the(:)
        integer(ik), intent(in) :: p_pol(:), p_qual(:)
        real(rk), intent(in) :: str_avg, dip_avg, rak_avg
        real(rk), intent(out) :: mfrac, stdr
        real(rk) :: M(3, 3), a(3), b(3)
        real(rk) :: strike, dip, rake, qcount, scount
        real(rk) :: azi, toff, wt, wo
        real(rk) :: bb1(3), bb2(3), bb3(3)
        real(rk) :: p_a1, p_a2, p_a3, p_b1, p_b3
        real(rk) :: p_proj1, p_proj2, p_proj3, plen, pp_b1, pp_b2
        real(rk) :: phi_ang, theta_ang, p_amp, dot_val, pol
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

        call fpcoor(strike, dip, rake, bb3, bb1, 1)
        call cross(bb3, bb1, bb2)

        mfrac = 0.0_rk
        qcount = 0.0_rk
        scount = 0.0_rk

        do k = 1, npol
            call to_cart(p_the(k), p_azi(k), 1.0_rk, p_a1, p_a2, p_a3)
            p_b1 = bb1(1) * p_a1 + bb1(2) * p_a2 + bb1(3) * p_a3
            p_b3 = bb3(1) * p_a1 + bb3(2) * p_a2 + bb3(3) * p_a3
            p_proj1 = p_a1 - p_b3 * bb3(1)
            p_proj2 = p_a2 - p_b3 * bb3(2)
            p_proj3 = p_a3 - p_b3 * bb3(3)
            plen = sqrt(p_proj1 * p_proj1 + p_proj2 * p_proj2 + p_proj3 * p_proj3)
            pp_b1 = bb1(1) * p_proj1 + bb1(2) * p_proj2 + bb1(3) * p_proj3
            pp_b2 = bb2(1) * p_proj1 + bb2(2) * p_proj2 + bb2(3) * p_proj3
            phi_ang = atan2(pp_b2, pp_b1)
            theta_ang = acos(max(-1.0_rk, min(1.0_rk, p_b3)))
            p_amp = abs(sin(2 * theta_ang) * cos(phi_ang))
            wt = sqrt(p_amp)

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
            if (p_qual(k) == 0) then
                wo = 1.0_rk
            else
                wo = 0.5_rk
            end if
            if (pol * real(p_pol(k), rk) < 0.0_rk) mfrac = mfrac + wt * wo
            qcount = qcount + wt * wo
            scount = scount + wo
        end do

        if (qcount > 0.0_rk) mfrac = mfrac / qcount
        if (scount > 0.0_rk) then
            stdr = qcount / scount
        else
            stdr = 0.0_rk
        end if

    end subroutine get_misfit

    !> In-place ascending heap sort (Numerical Recipes heapsort).
    subroutine sort_asc(n, ra)
        integer(ik), intent(in) :: n
        real(rk), intent(inout) :: ra(:)
        integer(ik) :: l, ir, i, j
        real(rk) :: rra
        if (n <= 1) return
        l = n / 2 + 1
        ir = n
        do
            if (l > 1) then
                l = l - 1
                rra = ra(l)
            else
                rra = ra(ir)
                ra(ir) = ra(1)
                ir = ir - 1
                if (ir == 1) then
                    ra(1) = rra
                    return
                end if
            end if
            i = l
            j = l + l
            do while (j <= ir)
                if (j < ir) then
                    if (ra(j) < ra(j + 1)) j = j + 1
                end if
                if (rra < ra(j)) then
                    ra(i) = ra(j)
                    i = j
                    j = j + j
                else
                    j = ir + 1
                end if
            end do
            ra(i) = rra
        end do
    end subroutine sort_asc

end module hash_misfit
