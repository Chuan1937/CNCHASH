!> Geometry utilities: spherical<->Cartesian transforms, fault-plane
!> coordinate conversions, vector primitives.
!>
!> Mathematical equivalents of HASH_complete/src/subs/util_subs.f
!> (TO_CAR, FPCOOR, CROSS), written in Modern Fortran with explicit
!> typing and no global state.
module hash_geometry
    use hash_kinds, only: rk, ik
    implicit none

    real(rk), parameter :: pi = 3.14159265358979323846_rk
    real(rk), parameter :: deg_to_rad = pi / 180.0_rk
    real(rk), parameter :: rad_to_deg = 180.0_rk / pi

    private
    public :: pi, deg_to_rad, rad_to_deg
    public :: to_cart, fpcoor, cross, normalize3

contains

    !> Transform spherical coordinates to Cartesian.
    !> theta = takeoff angle (deg, from vertical, up=0), phi = azimuth (deg E of N).
    !> (x=north, y=east, z=down)
    pure subroutine to_cart(theta, phi, r, x, y, z)
        real(rk), intent(in) :: theta, phi, r
        real(rk), intent(out) :: x, y, z
        real(rk) :: tr, pr
        tr = theta * deg_to_rad
        pr = phi * deg_to_rad
        z = -r * cos(tr)
        x = r * sin(tr) * cos(pr)
        y = r * sin(tr) * sin(pr)
    end subroutine to_cart

    !> Cross product v3 = v1 x v2.
    pure subroutine cross(v1, v2, v3)
        real(rk), intent(in) :: v1(3), v2(3)
        real(rk), intent(out) :: v3(3)
        v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
        v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
        v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
    end subroutine cross

    !> Normalize a 3-vector in place; no-op for zero-length input.
    pure subroutine normalize3(v)
        real(rk), intent(inout) :: v(3)
        real(rk) :: n
        n = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
        if (n > 0.0_rk) v = v / n
    end subroutine normalize3

    !> Convert between fault-plane coordinates and strike/dip/rake.
    !> idir = 1: (strike,dip,rake) -> (fnorm,slip)
    !> idir = 2: (fnorm,slip) -> (strike,dip,rake)
    !> Reference: Aki and Richards, p. 115. (x=north, y=east, z=down)
    pure subroutine fpcoor(strike, dip, rake, fnorm, slip, idir)
        real(rk), intent(inout) :: strike, dip, rake
        real(rk), intent(inout) :: fnorm(3), slip(3)
        integer(ik), intent(in) :: idir
        real(rk) :: phi, del, lam, a, clam, slam

        if (idir == 1) then
            phi = strike * deg_to_rad
            del = dip * deg_to_rad
            lam = rake * deg_to_rad
            fnorm(1) = -sin(del) * sin(phi)
            fnorm(2) = sin(del) * cos(phi)
            fnorm(3) = -cos(del)
            slip(1) = cos(lam) * cos(phi) + cos(del) * sin(lam) * sin(phi)
            slip(2) = cos(lam) * sin(phi) - cos(del) * sin(lam) * cos(phi)
            slip(3) = -sin(lam) * sin(del)
        else
            if ((1.0_rk - abs(fnorm(3))) <= 1.0e-7_rk) then
                del = 0.0_rk
                phi = atan2(-slip(1), slip(2))
                clam = cos(phi) * slip(1) + sin(phi) * slip(2)
                slam = sin(phi) * slip(1) - cos(phi) * slip(2)
                lam = atan2(slam, clam)
            else
                phi = atan2(-fnorm(1), fnorm(2))
                a = sqrt(fnorm(1) * fnorm(1) + fnorm(2) * fnorm(2))
                del = atan2(a, -fnorm(3))
                clam = cos(phi) * slip(1) + sin(phi) * slip(2)
                slam = -slip(3) / sin(del)
                lam = atan2(slam, clam)
                if (del > (0.5_rk * pi)) then
                    del = pi - del
                    phi = phi + pi
                    lam = -lam
                end if
            end if
            strike = phi * rad_to_deg
            if (strike < 0.0_rk) strike = strike + 360.0_rk
            dip = del * rad_to_deg
            rake = lam * rad_to_deg
            if (rake <= -180.0_rk) rake = rake + 360.0_rk
            if (rake > 180.0_rk) rake = rake - 360.0_rk
        end if
    end subroutine fpcoor

end module hash_geometry
