!> Runtime state: cached rotation grid and OpenMP thread control.
!>
!> The rotation grid is built once per dang value and then shared,
!> read-only, across all events and all threads. Building is guarded by
!> an OpenMP critical section (double-checked locking).
!>
!> Reentrancy note: this core is safe for OpenMP worker threads inside a
!> single native call, and for sequential calls. It is NOT reentrant for
!> independent concurrent API calls (the grid cache and the module-level
!> workspaces are process-global); the Python front-end serializes such
!> calls with a lock.
module hash_runtime
    use hash_kinds, only: rk, ik
    use hash_rotation, only: rotation_grid_t, build_rotation_grid
#ifdef _OPENMP
    use omp_lib, only: omp_set_num_threads, omp_get_max_threads, omp_in_parallel
#endif
    implicit none

    type(rotation_grid_t) :: grid_cache
    logical :: grid_valid = .false.

    private
    public :: ensure_rotation_grid_cached, grid_cache
    public :: hash_set_num_threads, hash_get_max_threads, hash_in_parallel_region

contains

    !> Rebuild the cached grid only if dang changed. Thread-safe.
    subroutine ensure_rotation_grid_cached(dang)
        real(rk), intent(in) :: dang
        if (.not. grid_valid .or. grid_cache%dang /= dang) then
            !$omp critical(hash_grid_build)
            if (.not. grid_valid .or. grid_cache%dang /= dang) then
                call build_rotation_grid(dang, grid_cache)
                grid_valid = .true.
            end if
            !$omp end critical(hash_grid_build)
        end if
    end subroutine ensure_rotation_grid_cached

    !> Set the OpenMP thread count for subsequent parallel regions.
    subroutine hash_set_num_threads(nt)
        integer(ik), intent(in) :: nt
#ifdef _OPENMP
        call omp_set_num_threads(max(1, nt))
#endif
    end subroutine hash_set_num_threads

    !> Current OpenMP max thread count (1 when built without OpenMP).
    function hash_get_max_threads() result(nt)
        integer(ik) :: nt
#ifdef _OPENMP
        nt = omp_get_max_threads()
#else
        nt = 1
#endif
    end function hash_get_max_threads

    !> True if called from inside an active parallel region.
    function hash_in_parallel_region() result(in_parallel)
        logical :: in_parallel
#ifdef _OPENMP
        in_parallel = omp_in_parallel()
#else
        in_parallel = .false.
#endif
    end function hash_in_parallel_region

end module hash_runtime
