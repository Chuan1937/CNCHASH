!> Rotation grid used by the focal-mechanism grid search.
!>
!> Port of the rotation-grid setup inside FOCALMC/FOCALAMP_MC
!> (HASH_complete/src/subs/fmech_subs.f lines 60-110). The grid is
!> cached per dang value; rebuilding happens only when dang changes.
!> All grid data is read-only once built, so the search loops can run
!> in parallel without locking.
module hash_rotation
    use hash_kinds, only: rk, ik
    use hash_geometry, only: deg_to_rad, cross
    implicit none

    type, public :: rotation_grid_t
        real(rk) :: dang = -1.0_rk
        integer(ik) :: nrot = 0
        ! Stored as (nrot, 3): each rotation's three components are
        ! contiguous in memory, so the search loops touch one cache line
        ! per rotation instead of three stride-nrot lines.
        real(rk), allocatable :: b1(:, :)
        real(rk), allocatable :: b2(:, :)
        real(rk), allocatable :: b3(:, :)
    end type rotation_grid_t

    private
    public :: build_rotation_grid, ensure_rotation_grid

contains

    !> Compute the total number of rotations for a given grid spacing.
    pure function count_rotations(dang) result(nrot)
        real(rk), intent(in) :: dang
        integer(ik) :: nrot
        integer(ik) :: ithe, iphi, max_the, max_zeta, numphi
        real(rk) :: the, rthe, fnumang, dphi

        max_the = int(90.1_rk / dang)
        max_zeta = int(179.9_rk / dang)
        nrot = 0
        do ithe = 0, max_the
            the = real(ithe, rk) * dang
            rthe = the * deg_to_rad
            fnumang = 360.0_rk / dang
            numphi = nint(fnumang * sin(rthe))
            if (numphi /= 0) then
                dphi = 360.0_rk / real(numphi, rk)
            else
                dphi = 10000.0_rk
            end if
            do iphi = 0, int(359.9_rk / dphi)
                nrot = nrot + (max_zeta + 1)
            end do
        end do
    end function count_rotations

    !> Build the rotation grid for a given angular spacing (degrees).
    subroutine build_rotation_grid(dang, grid)
        real(rk), intent(in) :: dang
        type(rotation_grid_t), intent(inout) :: grid
        integer(ik) :: max_the, max_zeta, irot, ithe, iphi, izeta
        integer(ik) :: numphi, nrot
        real(rk) :: the, rthe, costhe, sinthe, fnumang, dphi
        real(rk) :: phi, rphi, cosphi, sinphi
        real(rk) :: zeta, rzeta, coszeta, sinzeta
        real(rk) :: bb1(3), bb2(3), bb3(3)

        nrot = count_rotations(dang)
        if (allocated(grid%b1)) deallocate (grid%b1)
        if (allocated(grid%b2)) deallocate (grid%b2)
        if (allocated(grid%b3)) deallocate (grid%b3)
        allocate (grid%b1(nrot, 3), grid%b2(nrot, 3), grid%b3(nrot, 3))

        max_the = int(90.1_rk / dang)
        max_zeta = int(179.9_rk / dang)
        irot = 0
        do ithe = 0, max_the
            the = real(ithe, rk) * dang
            rthe = the * deg_to_rad
            costhe = cos(rthe)
            sinthe = sin(rthe)
            fnumang = 360.0_rk / dang
            numphi = nint(fnumang * sin(rthe))
            if (numphi /= 0) then
                dphi = 360.0_rk / real(numphi, rk)
            else
                dphi = 10000.0_rk
            end if
            do iphi = 0, int(359.9_rk / dphi)
                phi = real(iphi, rk) * dphi
                rphi = phi * deg_to_rad
                cosphi = cos(rphi)
                sinphi = sin(rphi)
                bb3(3) = costhe
                bb3(1) = sinthe * cosphi
                bb3(2) = sinthe * sinphi
                bb1(3) = -sinthe
                bb1(1) = costhe * cosphi
                bb1(2) = costhe * sinphi
                call cross(bb3, bb1, bb2)
                do izeta = 0, max_zeta
                    zeta = real(izeta, rk) * dang
                    rzeta = zeta * deg_to_rad
                    coszeta = cos(rzeta)
                    sinzeta = sin(rzeta)
                    irot = irot + 1
                    grid%b3(irot, 3) = bb3(3)
                    grid%b3(irot, 1) = bb3(1)
                    grid%b3(irot, 2) = bb3(2)
                    grid%b1(irot, 1) = bb1(1) * coszeta + bb2(1) * sinzeta
                    grid%b1(irot, 2) = bb1(2) * coszeta + bb2(2) * sinzeta
                    grid%b1(irot, 3) = bb1(3) * coszeta + bb2(3) * sinzeta
                    grid%b2(irot, 1) = bb2(1) * coszeta - bb1(1) * sinzeta
                    grid%b2(irot, 2) = bb2(2) * coszeta - bb1(2) * sinzeta
                    grid%b2(irot, 3) = bb2(3) * coszeta - bb1(3) * sinzeta
                end do
            end do
        end do
        grid%dang = dang
        grid%nrot = irot
    end subroutine build_rotation_grid

    !> Rebuild the grid only if dang differs from the cached value.
    subroutine ensure_rotation_grid(dang, grid)
        real(rk), intent(in) :: dang
        type(rotation_grid_t), intent(inout) :: grid
        if (grid%dang /= dang .or. grid%nrot == 0) then
            call build_rotation_grid(dang, grid)
        end if
    end subroutine ensure_rotation_grid

end module hash_rotation
