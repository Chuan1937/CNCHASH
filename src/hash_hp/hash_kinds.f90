!> Numerical kinds shared by all HASH-HP modules.
module hash_kinds
    use iso_fortran_env, only: real32, real64, int32, int64
    implicit none

    integer, parameter :: rk = real64
    integer, parameter :: ik = int32
    integer, parameter :: lk = int64
end module hash_kinds
