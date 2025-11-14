cmake_minimum_required(VERSION 3.15...4.0)
project(torus LANGUAGES CXX)
SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++17")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(Py_agama src/Py_agama.cpp)

set(SOURCE_FILES ${SOURCE_FILES}
    src/math_core.cpp
    src/math_fit.cpp
    src/math_fourier.cpp
    src/math_gausshermite.cpp
    src/math_geometry.cpp
    src/math_linalg.cpp
    src/math_ode.cpp
    src/math_optimization.cpp
    src/math_random.cpp
    src/math_sample.cpp
    src/math_specfunc.cpp
    src/math_sphharm.cpp
    src/math_spline.cpp
    src/particles_io.cpp
    src/actions_harmonicoscillator.cpp
    src/actions_newisochrone.cpp
    src/actions_newgenfnc.cpp
    src/actions_newtorus.cpp
    src/actions_spherical.cpp
    src/actions_staeckel.cpp
    src/coord.cpp
    src/cubature.cpp
    src/df_base.cpp
    src/df_disk.cpp
    src/df_ergodic.cpp
    src/df_factory.cpp
    src/df_halo.cpp
    src/df_interpolated.cpp
    src/df_spherical.cpp
    src/galaxymodel_base.cpp
    src/galaxymodel_densitygrid.cpp
    src/galaxymodel_jeans.cpp
    src/galaxymodel_losvd.cpp
    src/galaxymodel_selfconsistent.cpp
    src/galaxymodel_spherical.cpp
    src/galaxymodel_velocitysampler.cpp
    src/orbit.cpp
    src/map.cpp
    src/orbit_lyapunov.cpp
    src/obs_base.cpp
    src/obs_dust.cpp
    src/obs_los.cpp
    src/potential_analytic.cpp
    src/potential_bars.cpp
    src/potential_base.cpp
    src/potential_composite.cpp
    src/potential_cylspline.cpp
    src/potential_dehnen.cpp
    src/potential_disk.cpp
    src/potential_factory.cpp
    src/potential_ferrers.cpp
    src/potential_galpot.cpp
    src/potential_king.cpp
    src/potential_interpolators.cpp
    src/potential_multipole.cpp
    src/potential_perfect_ellipsoid.cpp
    src/potential_spheroid.cpp
    src/potential_utils.cpp
    src/utils.cpp
    src/utils_config.cpp
)
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
add_library(lib STATIC ${SOURCE_FILES})
target_link_libraries(Py_agama PUBLIC lib)
SET(PYBIND_INCLUDE_DIR ${pybind11_INCLUDE_DIR})
SET(PYTHON_INCLUDE_DIR ${pybind11_INCLUDE_DIRS})
include_directories(${MAKE_CURRENT_LIST_DIR})
include_directories(${PYBIND_INCLUDE_DIR})
include_directories(${PYTHON_INCLUDE_DIR})
message(STATUS, "Looking for GsL(Required)")
find_package(GSL REQUIRED)
target_link_libraries(Py_agama PUBLIC GSL::gsl GSL::gslcblas)
message(STATUS, "Looking for OpenMP")
find_package(OpenMP)
if(OPENMP_FOUND)
	target_link_libraries(Py_agama PUBLIC OpenMP::OpenMP_CXX)
	message(STATUS, "Found OpenMP")
else()
	message(STATUS, "Could not find OpenMP. This will slow down performance")
endif()
