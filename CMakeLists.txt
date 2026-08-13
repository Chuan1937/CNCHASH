# CNCHASH native backend (Modern Fortran + OpenMP)
#
# Builds libhashhp, the high-performance HASH-HP core. The backend is
# OPTIONAL for the Python package: cnchash falls back to the Numba
# backend when this library is unavailable.
#
# Options:
#   CNCHASH_FAST_MATH      - enable -ffast-math (default OFF; scientific
#                            correctness first)
#   CNCHASH_NATIVE_MARCH   - add -march=native (default ON for local
#                            performance; portable builds should set OFF)
#   CNCHASH_BUILD_TESTS    - build the native C test executable
cmake_minimum_required(VERSION 3.16)

project(hashhp LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

option(CNCHASH_FAST_MATH "Enable -ffast-math" OFF)
option(CNCHASH_NATIVE_MARCH "Enable -march=native" ON)
option(CNCHASH_BUILD_TESTS "Build native tests" OFF)

find_package(OpenMP QUIET)
if(NOT OpenMP_Fortran_FOUND)
    message(WARNING "OpenMP not found: building serial native backend")
endif()

file(GLOB HASHHP_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/hash_hp/*.f90")

add_library(hashhp SHARED ${HASHHP_SOURCES})

if(OpenMP_Fortran_FOUND)
    target_link_libraries(hashhp PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Intel")
    target_compile_options(hashhp PRIVATE -O3 -cpp)
    if(CNCHASH_FAST_MATH)
        target_compile_options(hashhp PRIVATE -ffast-math)
    endif()
    if(CNCHASH_NATIVE_MARCH AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64|aarch64")
        target_compile_options(hashhp PRIVATE -march=native)
    endif()
    # Thread-private workspace in a shared library needs initial-exec TLS.
    target_compile_options(hashhp PRIVATE -fPIC -ftls-model=initial-exec)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
    target_compile_options(hashhp PRIVATE -O3)
endif()

set_target_properties(hashhp PROPERTIES
    VERSION 1.0.0
    SOVERSION 1
    OUTPUT_NAME hashhp
    POSITION_INDEPENDENT_CODE ON
)

install(TARGETS hashhp LIBRARY DESTINATION lib)

if(CNCHASH_BUILD_TESTS)
    enable_testing()
    add_executable(c_api_smoke tests/c_api_smoke.c)
    target_link_libraries(c_api_smoke PRIVATE hashhp)
    add_test(NAME c_api_smoke COMMAND c_api_smoke)
endif()
