!> Uncertainty analysis: average mechanism, minimum rotation angle,
!> probability, and RMS angular differences (MECH_PROB).
!>
!> Port of HASH_complete/src/subs/uncert_subs.f (MECH_ROT, MECH_AVG,
!> MECH_PROB). Unlike the CNCHASH Numba approximation, this is a
!> faithful port: the full 4-combination minimum rotation search
!> (Kagan-style) is used everywhere, outlier removal uses the exact
!> maxrot criterion, and the orthogonalization uses the weighted
!> fract1 split with maxmisf = 0.01 deg.
module hash_uncertainty
    use hash_kinds, only: rk, ik
    use hash_geometry, only: pi, deg_to_rad, rad_to_deg, fpcoor, cross
    implicit none

    private
    public :: mech_rot, mech_avg, mech_prob

contains

    !> Minimum rotation angle between two mechanisms, trying the 4 nodal
    !> plane combinations. As in the original MECH_ROT, mech B vectors
    !> (n2, s2) are mutated to the best-matching orientation on return.
    subroutine mech_rot(n1, s1, n2, s2, rota)
        real(rk), intent(in) :: n1(3), s1(3)
        real(rk), intent(inout) :: n2(3), s2(3)
        real(rk), intent(out) :: rota
        real(rk) :: B1(3), B2(3)
        real(rk) :: n2t(3), s2t(3)
        real(rk) :: rotemp(4)
        real(rk) :: phi(3), theta(3), qdot(3)
        real(rk) :: dn(3, 3), scale(3), R(3), n1v(3), n2v(3)
        real(rk) :: qmins, qmindif, cval, sval, rnorm, dval
        integer(ik) :: iter, iout, iuse, k, i, j
        integer(ik) :: irot

        do iter = 1, 4
            if (iter < 3) then
                n2t = n2
                s2t = s2
            else
                n2t = s2
                s2t = n2
            end if
            if (iter == 2 .or. iter == 4) then
                n2t = -n2t
                s2t = -s2t
            end if

            call cross(n1, s1, B1)
            call cross(n2t, s2t, B2)

            dval = n1(1) * n2t(1) + n1(2) * n2t(2) + n1(3) * n2t(3)
            phi(1) = acos(max(-1.0_rk, min(1.0_rk, dval)))
            dval = s1(1) * s2t(1) + s1(2) * s2t(2) + s1(3) * s2t(3)
            phi(2) = acos(max(-1.0_rk, min(1.0_rk, dval)))
            dval = B1(1) * B2(1) + B1(2) * B2(2) + B1(3) * B2(3)
            phi(3) = acos(max(-1.0_rk, min(1.0_rk, dval)))

            if (phi(1) < 1.0e-4_rk .and. phi(2) < 1.0e-4_rk .and. phi(3) < 1.0e-4_rk) then
                rotemp(iter) = 0.0_rk
            else if (phi(1) < 1.0e-4_rk) then
                rotemp(iter) = rad_to_deg * phi(2)
            else if (phi(2) < 1.0e-4_rk) then
                rotemp(iter) = rad_to_deg * phi(3)
            else if (phi(3) < 1.0e-4_rk) then
                rotemp(iter) = rad_to_deg * phi(1)
            else
                ! Difference vectors; the rotation axis must be orthogonal
                ! to all three.
                dn(1, 1) = n1(1) - n2t(1)
                dn(2, 1) = n1(2) - n2t(2)
                dn(3, 1) = n1(3) - n2t(3)
                dn(1, 2) = s1(1) - s2t(1)
                dn(2, 2) = s1(2) - s2t(2)
                dn(3, 2) = s1(3) - s2t(3)
                dn(1, 3) = B1(1) - B2(1)
                dn(2, 3) = B1(2) - B2(2)
                dn(3, 3) = B1(3) - B2(3)
                do j = 1, 3
                    scale(j) = sqrt(dn(1, j)**2 + dn(2, j)**2 + dn(3, j)**2)
                    if (scale(j) > 0.0_rk) dn(:, j) = dn(:, j) / scale(j)
                end do
                qdot(3) = dn(1, 1) * dn(1, 2) + dn(2, 1) * dn(2, 2) + dn(3, 1) * dn(3, 2)
                qdot(2) = dn(1, 1) * dn(1, 3) + dn(2, 1) * dn(2, 3) + dn(3, 1) * dn(3, 3)
                qdot(1) = dn(1, 2) * dn(1, 3) + dn(2, 2) * dn(2, 3) + dn(3, 2) * dn(3, 3)

                ! Use the two largest difference vectors, unless parallel.
                iout = 0
                do i = 1, 3
                    if (qdot(i) > 0.9999_rk) iout = i
                end do
                if (iout == 0) then
                    qmins = 10000.0_rk
                    do i = 1, 3
                        if (scale(i) < qmins) then
                            qmins = scale(i)
                            iout = i
                        end if
                    end do
                end if

                k = 0
                do j = 1, 3
                    if (j /= iout) then
                        k = k + 1
                        if (k == 1) then
                            n1v = dn(:, j)
                        else
                            n2v = dn(:, j)
                        end if
                    end if
                end do

                call cross(n1v, n2v, R)
                rnorm = sqrt(R(1)**2 + R(2)**2 + R(3)**2)
                if (rnorm > 0.0_rk) R = R / rnorm

                dval = n1(1) * R(1) + n1(2) * R(2) + n1(3) * R(3)
                theta(1) = acos(max(-1.0_rk, min(1.0_rk, dval)))
                dval = s1(1) * R(1) + s1(2) * R(2) + s1(3) * R(3)
                theta(2) = acos(max(-1.0_rk, min(1.0_rk, dval)))
                dval = B1(1) * R(1) + B1(2) * R(2) + B1(3) * R(3)
                theta(3) = acos(max(-1.0_rk, min(1.0_rk, dval)))

                qmindif = 1000.0_rk
                iuse = 1
                do i = 1, 3
                    dval = abs(theta(i) - pi / 2.0_rk)
                    if (dval < qmindif) then
                        qmindif = dval
                        iuse = i
                    end if
                end do

                cval = cos(phi(iuse)) - cos(theta(iuse)) * cos(theta(iuse))
                sval = sin(theta(iuse)) * sin(theta(iuse))
                if (abs(sval) > 1.0e-10_rk) then
                    dval = max(-1.0_rk, min(1.0_rk, cval / sval))
                    rotemp(iter) = rad_to_deg * acos(dval)
                else
                    rotemp(iter) = rad_to_deg * phi(iuse)
                end if
            end if
        end do

        ! Minimum rotation over the 4 combos; mutate mech B to the best.
        rota = 180.0_rk
        irot = 1
        do iter = 1, 4
            if (abs(rotemp(iter)) < rota) then
                rota = abs(rotemp(iter))
                irot = iter
            end if
        end do
        if (irot >= 3) then
            call swap3(n2, s2)
        end if
        if (irot == 2 .or. irot == 4) then
            n2 = -n2
            s2 = -s2
        end if
    end subroutine mech_rot

    !> Average focal mechanism of a set (with plane matching and
    !> orthogonalization), per original MECH_AVG.
    subroutine mech_avg(nf, norm1, norm2, norm1_avg, norm2_avg)
        integer(ik), intent(in) :: nf
        real(rk), intent(in) :: norm1(:, :), norm2(:, :)
        real(rk), intent(out) :: norm1_avg(3), norm2_avg(3)
        real(rk) :: ref1(3), ref2(3), temp1(3), temp2(3)
        real(rk) :: rot_angle, ln1, ln2
        real(rk) :: d11, d22, a11, a22, avang1, avang2
        real(rk) :: dot1, misf, fract1, theta1, theta2
        real(rk), parameter :: maxmisf = 0.01_rk
        integer(ik) :: i, j, icount

        if (nf <= 1) then
            norm1_avg = norm1(:, 1)
            norm2_avg = norm2(:, 1)
            return
        end if

        norm1_avg = norm1(:, 1)
        norm2_avg = norm2(:, 1)
        ref1 = norm1(:, 1)
        ref2 = norm2(:, 1)

        do i = 2, nf
            temp1 = norm1(:, i)
            temp2 = norm2(:, i)
            rot_angle = 0.0_rk
            call mech_rot(ref1, ref2, temp1, temp2, rot_angle)
            norm1_avg = norm1_avg + temp1
            norm2_avg = norm2_avg + temp2
        end do
        ln1 = sqrt(sum(norm1_avg**2))
        ln2 = sqrt(sum(norm2_avg**2))
        if (ln1 > 0.0_rk) norm1_avg = norm1_avg / ln1
        if (ln2 > 0.0_rk) norm2_avg = norm2_avg / ln2

        ! RMS angular differences of the (unmatched) vectors, used to split
        ! the orthogonalization correction between the two planes.
        avang1 = 0.0_rk
        avang2 = 0.0_rk
        do i = 1, nf
            temp1 = norm1(:, i)
            temp2 = norm2(:, i)
            d11 = temp1(1) * norm1_avg(1) + temp1(2) * norm1_avg(2) + temp1(3) * norm1_avg(3)
            d22 = temp2(1) * norm2_avg(1) + temp2(2) * norm2_avg(2) + temp2(3) * norm2_avg(3)
            a11 = acos(max(-1.0_rk, min(1.0_rk, d11)))
            a22 = acos(max(-1.0_rk, min(1.0_rk, d22)))
            avang1 = avang1 + a11 * a11
            avang2 = avang2 + a22 * a22
        end do
        avang1 = sqrt(avang1 / real(nf, rk))
        avang2 = sqrt(avang2 / real(nf, rk))

        if ((avang1 + avang2) < 0.0001_rk) return

        fract1 = avang1 / (avang1 + avang2)
        do icount = 1, 100
            dot1 = norm1_avg(1) * norm2_avg(1) + norm1_avg(2) * norm2_avg(2) &
                   + norm1_avg(3) * norm2_avg(3)
            misf = 90.0_rk - acos(max(-1.0_rk, min(1.0_rk, dot1))) * rad_to_deg
            if (abs(misf) <= maxmisf) return
            theta1 = misf * fract1 * deg_to_rad
            theta2 = misf * (1.0_rk - fract1) * deg_to_rad
            do j = 1, 3
                d11 = norm1_avg(j)
                norm1_avg(j) = norm1_avg(j) - norm2_avg(j) * sin(theta1)
                norm2_avg(j) = norm2_avg(j) - d11 * sin(theta2)
            end do
            ln1 = sqrt(sum(norm1_avg**2))
            ln2 = sqrt(sum(norm2_avg**2))
            if (ln1 > 0.0_rk) norm1_avg = norm1_avg / ln1
            if (ln2 > 0.0_rk) norm2_avg = norm2_avg / ln2
        end do
    end subroutine mech_avg

    !> Preferred mechanism(s) from a set of acceptable mechanisms, with
    !> probability and RMS uncertainty (per original MECH_PROB).
    subroutine mech_prob(nf, norm1in, norm2in, cangle, prob_max, nsltn, &
                         str_avg, dip_avg, rak_avg, prob, rms_diff)
        integer(ik), intent(in) :: nf
        real(rk), intent(in) :: norm1in(:, :), norm2in(:, :)
        real(rk), intent(in) :: cangle, prob_max
        integer(ik), intent(out) :: nsltn
        real(rk), intent(out) :: str_avg(5), dip_avg(5), rak_avg(5)
        real(rk), intent(out) :: prob(5), rms_diff(2, 5)
        real(rk), allocatable :: norm1(:, :), norm2(:, :)
        real(rk) :: norm1_avg(3), norm2_avg(3), temp1(3), temp2(3)
        real(rk) :: rota(nf), rot_angle, d11, d22, a11, a22
        integer(ik) :: i, j, imult, icount, nfault, nc, imax
        real(rk) :: maxrot

        if (nf <= 1) then
            temp1 = norm1in(:, 1)
            temp2 = norm2in(:, 1)
            call fpcoor(str_avg(1), dip_avg(1), rak_avg(1), temp1, temp2, 2)
            prob(1) = 1.0_rk
            rms_diff(1, 1) = 0.0_rk
            rms_diff(2, 1) = 0.0_rk
            nsltn = 1
            return
        end if

        allocate (norm1(3, nf), norm2(3, nf))
        norm1 = norm1in
        norm2 = norm2in

        nfault = nf
        nc = nf
        do imult = 1, 5
            if (nc < 1) exit

            ! Repeatedly remove the mechanism with the largest angular
            ! difference from the running average, until all are within
            ! cangle of the average.
            do icount = 1, nf
                call mech_avg(nc, norm1, norm2, norm1_avg, norm2_avg)
                do i = 1, nc
                    temp1 = norm1(:, i)
                    temp2 = norm2(:, i)
                    rot_angle = 0.0_rk
                    call mech_rot(norm1_avg, norm2_avg, temp1, temp2, rot_angle)
                    rota(i) = rot_angle
                end do
                maxrot = 0.0_rk
                imax = 1
                do i = 1, nc
                    if (abs(rota(i)) > maxrot) then
                        maxrot = abs(rota(i))
                        imax = i
                    end if
                end do
                if (maxrot <= cangle) exit
                ! Move the outlier to the end.
                temp1 = norm1(:, imax)
                temp2 = norm2(:, imax)
                do j = imax, nc - 1
                    norm1(:, j) = norm1(:, j + 1)
                    norm2(:, j) = norm2(:, j + 1)
                end do
                norm1(:, nc) = temp1
                norm2(:, nc) = temp2
                nc = nc - 1
            end do

            prob(imult) = real(nc, rk) / real(nfault, rk)

            if (imult > 1 .and. prob(imult) < prob_max) exit

            ! Set up for next round: move the outliers to the front.
            do j = 1, nfault - nc
                norm1(:, j) = norm1(:, j + nc)
                norm2(:, j) = norm2(:, j + nc)
            end do
            nc = nfault - nc

            ! RMS angular differences between each input mechanism and the
            ! average, after matching planes.
            rms_diff(1, imult) = 0.0_rk
            rms_diff(2, imult) = 0.0_rk
            do i = 1, nfault
                temp1 = norm1in(:, i)
                temp2 = norm2in(:, i)
                rot_angle = 0.0_rk
                call mech_rot(norm1_avg, norm2_avg, temp1, temp2, rot_angle)
                d11 = temp1(1) * norm1_avg(1) + temp1(2) * norm1_avg(2) + temp1(3) * norm1_avg(3)
                d22 = temp2(1) * norm2_avg(1) + temp2(2) * norm2_avg(2) + temp2(3) * norm2_avg(3)
                a11 = acos(max(-1.0_rk, min(1.0_rk, d11)))
                a22 = acos(max(-1.0_rk, min(1.0_rk, d22)))
                rms_diff(1, imult) = rms_diff(1, imult) + a11 * a11
                rms_diff(2, imult) = rms_diff(2, imult) + a22 * a22
            end do
            rms_diff(1, imult) = rad_to_deg * sqrt(rms_diff(1, imult) / real(nfault, rk))
            rms_diff(2, imult) = rad_to_deg * sqrt(rms_diff(2, imult) / real(nfault, rk))

            call fpcoor(str_avg(imult), dip_avg(imult), rak_avg(imult), &
                        norm1_avg, norm2_avg, 2)
        end do

        nsltn = imult - 1
    end subroutine mech_prob

    !> Swap two 3-vectors.
    subroutine swap3(a, b)
        real(rk), intent(inout) :: a(3), b(3)
        real(rk) :: t(3)
        t = a
        a = b
        b = t
    end subroutine swap3

end module hash_uncertainty
