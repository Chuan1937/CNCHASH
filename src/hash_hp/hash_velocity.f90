!> 1D velocity-model ray tracing and takeoff-angle tables.
!>
!> Port of HASH_complete/src/subs/vel_subs.f (LAYERTRACE, MK_TABLE,
!> GET_TTS) with no I/O: velocity models enter as plain arrays.
!>
!> Deliberate deviations (documented):
!>   - GET_TTS depth-range check uses the actual table depth count instead
!>     of the original's uninitialized d(nd0) element (original bug).
!>   - asin argument clamped to [-1,1] to avoid NaN in double precision.
module hash_velocity
    use hash_kinds, only: rk, ik
    implicit none

    real(rk), parameter :: pi = 3.14159265358979323846_rk

    private
    public :: layer_trace, build_takeoff_table, get_tts

    integer(ik), parameter :: nray0 = 20001
    integer(ik), parameter :: max_model_pts = 2000

contains

    !> Ray trace through a single layer (original LAYERTRACE, double
    !> precision internally).
    !>
    !> Inputs: p (horizontal slowness), h (layer thickness),
    !>         utop, ubot (slowness at top/bottom), imth (method)
    !> Outputs: dx (range offset), dt (travel time), irtr (return code)
    subroutine layer_trace(p1, h1, utop1, ubot1, imth, dx1, dt1, irtr)
        real(rk), intent(in) :: p1, h1, utop1, ubot1
        integer(ik), intent(in) :: imth
        real(rk), intent(out) :: dx1, dt1
        integer(ik), intent(out) :: irtr
        real(rk) :: p, h, utop, ubot, u, y, q, qs, qr, b, etau, ex, dtau, dx, dt
        real(rk) :: zz

        p = p1
        h = h1
        utop = utop1
        ubot = ubot1

        if (h == 0.0_rk) then
            dx1 = 0.0_rk
            dt1 = 0.0_rk
            irtr = -1
            return
        end if

        u = utop
        y = u - p
        if (y <= 0.0_rk) then
            dx1 = 0.0_rk
            dt1 = 0.0_rk
            irtr = 0
            return
        end if

        q = y * (u + p)
        qs = sqrt(q)

        if (imth == 2) then
            y = u + qs
            if (p /= 0.0_rk) y = y / p
            qr = log(y)
        else if (imth == 3) then
            qr = atan2(qs, p)
        end if

        if (imth == 1) then
            b = -(utop**2 - ubot**2) / (2.0_rk * h)
        else if (imth == 2) then
            b = -((1.0_rk / utop) - (1.0_rk / ubot)) / h
        else
            b = -log(ubot / utop) / h
        end if

        if (b == 0.0_rk) then
            b = 1.0_rk / h
            etau = qs
            ex = p / qs
            irtr = 1
            go to 160
        end if

        if (imth == 1) then
            etau = -q * qs / 3.0_rk
            ex = -qs * p
        else if (imth == 2) then
            ex = qs / u
            etau = qr - ex
            if (p /= 0.0_rk) ex = ex / p
        else
            etau = qs - p * qr
            ex = qr
        end if

        u = ubot
        if (u <= p) then
            irtr = 2
            go to 160
        end if
        irtr = 1
        q = (u - p) * (u + p)
        qs = sqrt(q)

        if (imth == 1) then
            etau = etau + q * qs / 3.0_rk
            ex = ex + qs * p
        else if (imth == 2) then
            y = u + qs
            zz = qs / u
            etau = etau + zz
            if (p /= 0.0_rk) then
                y = y / p
                zz = zz / p
            end if
            qr = log(y)
            etau = etau - qr
            ex = ex - zz
        else
            qr = atan2(qs, p)
            etau = etau - qs + p * qr
            ex = ex - qr
        end if

160     dx = ex / b
        dtau = etau / b
        dt = dtau + p * dx
        dx1 = dx
        dt1 = dt
    end subroutine layer_trace

    !> Build a takeoff-angle table from a 1D velocity model (original
    !> MK_TABLE, no file I/O).
    !>
    !> Inputs:
    !>   z(1:npts), alpha(1:npts) - model depths (km) and P velocities (km/s)
    !>   npts              - number of model points
    !>   del1, del2, del3  - distance grid: start, end, step (km)
    !>   dep1, dep2, dep3  - depth grid: start, end, step (km)
    !>   pmin, nump        - ray-parameter grid start and count
    !> Outputs:
    !>   table(ndel,ndep)  - takeoff angles (degrees)
    !>   delttab(ndel)     - distance grid
    !>   deptab(ndep)      - depth grid
    !>   ndel, ndep        - grid sizes
    subroutine build_takeoff_table(z, alpha, npts, del1, del2, del3, dep1, &
                                   dep2, dep3, pmin, nump, table, delttab, &
                                   deptab, ndel, ndep)
        real(rk), intent(in) :: z(:), alpha(:)
        integer(ik), intent(in) :: npts
        real(rk), intent(in) :: del1, del2, del3, dep1, dep2, dep3, pmin
        integer(ik), intent(in) :: nump
        real(rk), intent(out) :: table(:, :)
        real(rk), intent(out) :: delttab(:), deptab(:)
        integer(ik), intent(out) :: ndel, ndep
        real(rk) :: slow(max_model_pts)
        real(rk) :: zz(max_model_pts), aa(max_model_pts)
        real(rk), allocatable :: deltab(:), tttab(:), ptab(:)
        real(rk), allocatable :: depxcor(:, :), depucor(:, :), deptcor(:, :)
        real(rk), allocatable :: xsave(:), tsave(:), psave(:), usave(:)
        real(rk) :: degrad, pmax, plongcut, pstep, p, x, t, h
        real(rk) :: dx, dt, xdeg, tmin, x2, t2, x1, t1, frac, scr1, angle
        real(rk) :: qtempdep2, dep, del, tt
        integer(ik) :: itab, idep, idel, i, j, np, imth, irtr, icount
        integer(ik) :: npts2, npts_old, npmax, i2, ncount, model_npts

        degrad = 180.0_rk / pi
        model_npts = npts
        if (model_npts < 2 .or. model_npts > max_model_pts) then
            ndel = 0
            ndep = 0
            return
        end if

        ! Depth grid.
        qtempdep2 = dep2 + dep3 / 20.0_rk
        ndep = int((qtempdep2 - dep1) / dep3) + 1
        do idep = 1, ndep
            dep = dep1 + dep3 * real(idep - 1, rk)
            deptab(idep) = dep
        end do
        ndel = int((del2 - del1) / del3) + 1
        allocate (deltab(nray0), tttab(nray0), ptab(nray0))
        allocate (xsave(nray0), tsave(nray0), psave(nray0), usave(nray0))
        allocate (depxcor(nray0, max(ndep, 1)))
        allocate (depucor(nray0, max(ndep, 1)))
        allocate (deptcor(nray0, max(ndep, 1)))

        ! Copy model and insert table depths as new model points.
        zz(1:model_npts) = z(1:model_npts)
        aa(1:model_npts) = alpha(1:model_npts)
        npts2 = model_npts
        npts_old = npts2
        do i = npts_old, 2, -1
            do idep = ndep, 1, -1
                if ((zz(i - 1) <= (deptab(idep) - 0.1_rk)) .and. &
                    (zz(i) >= (deptab(idep) + 0.1_rk))) then
                    npts2 = npts2 + 1
                    do j = npts2, i + 1, -1
                        zz(j) = zz(j - 1)
                        aa(j) = aa(j - 1)
                    end do
                    zz(i) = deptab(idep)
                    frac = (zz(i) - zz(i - 1)) / (zz(i + 1) - zz(i - 1))
                    aa(i) = aa(i - 1) + frac * (aa(i + 1) - aa(i - 1))
                end if
            end do
        end do
        do i = 1, npts2
            slow(i) = 1.0_rk / aa(i)
        end do
        pmax = slow(1)
        plongcut = slow(npts2)
        pstep = (pmax - pmin) / real(nump, rk)

        ! P-wave ray tracing.
        npmax = int((pmax + pstep / 2.0_rk - pmin) / pstep) + 1
        do np = 1, npmax
            p = pmin + pstep * real(np - 1, rk)
            ptab(np) = p
            x = 0.0_rk
            t = 0.0_rk
            imth = 3
            do idep = 1, ndep
                if (deptab(idep) == 0.0_rk) then
                    depxcor(np, idep) = 0.0_rk
                    deptcor(np, idep) = 0.0_rk
                    depucor(np, idep) = slow(1)
                else
                    depxcor(np, idep) = -999.0_rk
                    deptcor(np, idep) = -999.0_rk
                    depucor(np, idep) = -999.0_rk
                end if
            end do
            do i = 1, npts2 - 1
                if (zz(i) >= 9999.0_rk) then
                    deltab(np) = -999.0_rk
                    tttab(np) = -999.0_rk
                    go to 200
                end if
                h = zz(i + 1) - zz(i)
                if (h == 0.0_rk) cycle
                call layer_trace(p, h, slow(i), slow(i + 1), imth, dx, dt, irtr)
                x = x + dx
                t = t + dt
                if (irtr == 0 .or. irtr == 2) go to 105
                xdeg = x
                tmin = t
                do idep = 1, ndep
                    if (abs(zz(i + 1) - deptab(idep)) < 0.1_rk) then
                        depxcor(np, idep) = xdeg
                        deptcor(np, idep) = tmin
                        depucor(np, idep) = slow(i + 1)
                    end if
                end do
            end do
105         xdeg = 2.0_rk * x
            tmin = 2.0_rk * t
            deltab(np) = xdeg
            tttab(np) = tmin
200         continue
        end do

        ! Build the takeoff-angle table for each source depth.
        do idep = 1, ndep
            icount = 0
            x1 = -999.0_rk
            if (deptab(idep) == 0.0_rk) then
                i2 = np
                go to 223
            end if
            do i = 1, np
                x2 = depxcor(i, idep)
                if (x2 == -999.0_rk) exit
                if (x2 <= x1) exit
                t2 = deptcor(i, idep)
                icount = icount + 1
                xsave(icount) = x2
                tsave(icount) = t2
                psave(icount) = -ptab(i)
                usave(icount) = depucor(i, idep)
                x1 = x2
            end do
            i2 = i - 1
223         do i = i2, 1, -1
                if (depxcor(i, idep) == -999.0_rk) cycle
                if (deltab(i) == -999.0_rk) cycle
                x2 = deltab(i) - depxcor(i, idep)
                t2 = tttab(i) - deptcor(i, idep)
                icount = icount + 1
                xsave(icount) = x2
                tsave(icount) = t2
                psave(icount) = ptab(i)
                usave(icount) = depucor(i, idep)
            end do
            ncount = icount

            do idel = 1, ndel
                del = del1 + del3 * real(idel - 1, rk)
                delttab(idel) = del
                tt = 999.0_rk
                do i = 2, ncount
                    x1 = xsave(i - 1)
                    x2 = xsave(i)
                    if (x1 > del .or. x2 < del) go to 230
                    if (psave(i) > 0.0_rk .and. psave(i) < plongcut) go to 230
                    frac = (del - x1) / (x2 - x1)
                    t1 = tsave(i - 1) + frac * (tsave(i) - tsave(i - 1))
                    if (t1 < tt) then
                        tt = t1
                        scr1 = psave(i) / usave(i)
                        scr1 = max(-1.0_rk, min(1.0_rk, scr1))
                        angle = asin(scr1) * degrad
                        if (angle < 0.0_rk) then
                            angle = -angle
                        else
                            angle = 180.0_rk - angle
                        end if
                        table(idel, idep) = angle
                    end if
230                 continue
                end do
            end do
        end do

        if (delttab(1) == 0.0_rk) then
            do idep = 1, ndep
                table(1, idep) = 0.0_rk
            end do
        end if

    end subroutine build_takeoff_table

    !> Interpolate the takeoff angle from a table (original GET_TTS).
    !>
    !> Inputs: table(ndel,ndep), delttab, deptab, del (distance, km),
    !>         qdep (depth, km)
    !> Outputs: tt (takeoff angle, deg), iflag (-1 outside depth range,
    !>          0 interpolation, 1 extrapolation)
    subroutine get_tts(table, delttab, deptab, ndel, ndep, del, qdep, tt, iflag)
        real(rk), intent(in) :: table(:, :)
        real(rk), intent(in) :: delttab(:), deptab(:)
        integer(ik), intent(in) :: ndel, ndep
        real(rk), intent(in) :: del, qdep
        real(rk), intent(out) :: tt
        integer(ik), intent(out) :: iflag
        integer(ik) :: id1, id2, ix1, ix2, ix, id, ixbest1, ixbest2
        real(rk) :: xfrac, t1, t2, dfrac, xoffmin1, xoffmin2, xoff, xmid
        real(rk) :: tt1, tt2

        if (qdep < deptab(1) .or. qdep > deptab(ndep)) then
            iflag = -1
            tt = 999.0_rk
            return
        end if

        id1 = ndep - 1
        id2 = ndep
        do id = 2, ndep
            if (deptab(id) < qdep) cycle
            id1 = id - 1
            id2 = id
            exit
        end do

        ix1 = ndel - 1
        ix2 = ndel
        do ix = 2, ndel
            if (delttab(ix) < del) cycle
            ix1 = ix - 1
            ix2 = ix
            exit
        end do

        if (table(ix1, id1) /= 0.0_rk .and. table(ix1, id2) /= 0.0_rk .and. &
            table(ix2, id1) /= 0.0_rk .and. table(ix2, id2) /= 0.0_rk .and. &
            delttab(ix2) >= del) then
            iflag = 0
            xfrac = (del - delttab(ix1)) / (delttab(ix2) - delttab(ix1))
            t1 = table(ix1, id1) + xfrac * (table(ix2, id1) - table(ix1, id1))
            t2 = table(ix1, id2) + xfrac * (table(ix2, id2) - table(ix1, id2))
            dfrac = (qdep - deptab(id1)) / (deptab(id2) - deptab(id1))
            tt = t1 + dfrac * (t2 - t1)
            return
        end if

        ! Extrapolate.
        iflag = 1
        xoffmin1 = 999.0_rk
        xoffmin2 = 999.0_rk
        ixbest1 = -1
        ixbest2 = -1
        do ix = 2, ndel
            if (table(ix - 1, id1) /= 0.0_rk .and. table(ix, id1) /= 0.0_rk) then
                xmid = 0.5_rk * (delttab(ix - 1) + delttab(ix))
                xoff = abs(xmid - del)
                if (xoff < xoffmin1) then
                    xoffmin1 = xoff
                    ixbest1 = ix
                end if
            end if
            if (table(ix - 1, id2) /= 0.0_rk .and. table(ix, id2) /= 0.0_rk) then
                xmid = 0.5_rk * (delttab(ix - 1) + delttab(ix))
                xoff = abs(xmid - del)
                if (xoff < xoffmin2) then
                    xoffmin2 = xoff
                    ixbest2 = ix
                end if
            end if
        end do
        if (ixbest1 == -1 .or. ixbest2 == -1) then
            iflag = -1
            tt = 999.0_rk
            return
        end if
        xfrac = (del - delttab(ixbest1 - 1)) / (delttab(ixbest1) - delttab(ixbest1 - 1))
        t1 = table(ixbest1 - 1, id1)
        t2 = table(ixbest1, id1)
        tt1 = t1 + xfrac * (t2 - t1)
        xfrac = (del - delttab(ixbest2 - 1)) / (delttab(ixbest2) - delttab(ixbest2 - 1))
        t1 = table(ixbest2 - 1, id2)
        t2 = table(ixbest2, id2)
        tt2 = t1 + xfrac * (t2 - t1)
        dfrac = (qdep - deptab(id1)) / (deptab(id2) - deptab(id1))
        tt = tt1 + dfrac * (tt2 - tt1)
    end subroutine get_tts

end module hash_velocity
