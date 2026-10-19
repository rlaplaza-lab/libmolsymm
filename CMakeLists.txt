cmake_minimum_required(VERSION 3.20)
project(molsymm VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(MOLSYMM_BUILD_PYTHON "Build nanobind Python module" OFF)
option(MOLSYMM_BUILD_FORTRAN "Build Fortran module" OFF)
option(MOLSYMM_BUILD_TESTS "Build tests" ON)
option(MOLSYMM_BUILD_SHARED "Build shared library" ON)

# scikit-build-core: ship a single extension module with a static core.
if(DEFINED SKBUILD)
  set(MOLSYMM_BUILD_PYTHON ON CACHE BOOL "" FORCE)
  set(MOLSYMM_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(MOLSYMM_BUILD_FORTRAN OFF CACHE BOOL "" FORCE)
  set(MOLSYMM_BUILD_SHARED OFF CACHE BOOL "" FORCE)
endif()

include(FetchContent)
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  # Tarball avoids flaky git clones against gitlab.com under CI load.
  FetchContent_Declare(
    Eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2
    URL_HASH SHA256=b4c198460eba6f28d34894e3a5710998818515104d6e74e5cc331ce31e46e626
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  )
  set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(Eigen)
endif()

set(MOLSYMM_SOURCES
  src/masses.cpp
  src/molecule.cpp
  src/operations.cpp
  src/equivalence.cpp
  src/detect.cpp
  src/character_table.cpp
  src/subgroups.cpp
  src/symmetrize.cpp
  src/displacements.cpp
  src/c_api.cpp
)

if(MOLSYMM_BUILD_SHARED)
  add_library(molsymm SHARED ${MOLSYMM_SOURCES})
else()
  add_library(molsymm STATIC ${MOLSYMM_SOURCES})
endif()

add_library(molsymm::molsymm ALIAS molsymm)

target_include_directories(molsymm
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(molsymm PUBLIC Eigen3::Eigen)
# Treat Eigen headers as system includes so -Werror does not fire inside them.
get_target_property(_molsymm_eigen_inc Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
if(_molsymm_eigen_inc)
  target_include_directories(molsymm SYSTEM PUBLIC ${_molsymm_eigen_inc})
endif()
target_compile_features(molsymm PUBLIC cxx_std_20)
# Needed when a static core is linked into the Python extension module.
set_target_properties(molsymm PROPERTIES POSITION_INDEPENDENT_CODE ON)

function(molsymm_set_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE
      -Wall -Wextra -Wpedantic -Wshadow -Werror
    )
  endif()
endfunction()

molsymm_set_warnings(molsymm)

if(NOT DEFINED SKBUILD)
  include(GNUInstallDirs)
  install(TARGETS molsymm
    EXPORT MolsymmTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
  install(DIRECTORY include/molsymm DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(EXPORT MolsymmTargets
    FILE MolsymmTargets.cmake
    NAMESPACE molsymm::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Molsymm
  )

  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/MolsymmConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MolsymmConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/MolsymmConfig.cmake"
    @ONLY
  )
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/MolsymmConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MolsymmConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Molsymm
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/molsymm.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/molsymm.pc"
    @ONLY
  )
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/molsymm.pc"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
  )
endif()

if(MOLSYMM_BUILD_TESTS)
  enable_testing()
  add_executable(molsymm_tests tests/test_main.cpp)
  target_link_libraries(molsymm_tests PRIVATE molsymm)
  molsymm_set_warnings(molsymm_tests)
  target_compile_definitions(molsymm_tests PRIVATE
    MOLSYMM_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data"
  )
  add_test(NAME molsymm_tests COMMAND molsymm_tests)
endif()

if(NOT DEFINED SKBUILD)
  add_executable(analyze_xyz examples/analyze_xyz.cpp)
  target_link_libraries(analyze_xyz PRIVATE molsymm)
  molsymm_set_warnings(analyze_xyz)
endif()

if(MOLSYMM_BUILD_PYTHON)
  find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
  find_package(nanobind CONFIG QUIET)
  if(NOT nanobind_FOUND)
    FetchContent_Declare(
      nanobind
      GIT_REPOSITORY https://github.com/wjakob/nanobind.git
      GIT_TAG v2.2.0
      GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(nanobind)
  endif()
  nanobind_add_module(molsymm_py python/molsymm.cpp)
  target_link_libraries(molsymm_py PRIVATE molsymm)
  set_target_properties(molsymm_py PROPERTIES OUTPUT_NAME molsymm)
  # Nanobind headers trip -Wshadow/-Wpedantic; keep first-party warning set on C++ targets only.
  if(MOLSYMM_BUILD_TESTS)
    # Some interpreters (e.g. conda) ship an older libstdc++ via RPATH that cannot
    # load a library built with the system g++. Preload the compiler's libstdc++.
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.so.6
      OUTPUT_VARIABLE _molsymm_libstdcxx
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    get_filename_component(_molsymm_libstdcxx_real "${_molsymm_libstdcxx}" REALPATH)
    add_test(
      NAME molsymm_python_test
      COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_python.py
    )
    set_tests_properties(molsymm_python_test PROPERTIES
      ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:molsymm_py>;MOLSYMM_TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests/data;LD_PRELOAD=${_molsymm_libstdcxx_real}"
    )
  endif()
  if(DEFINED SKBUILD)
    install(TARGETS molsymm_py LIBRARY DESTINATION . RUNTIME DESTINATION . ARCHIVE DESTINATION .)
  endif()
endif()

if(MOLSYMM_BUILD_FORTRAN)
  enable_language(Fortran)
  add_library(molsymm_fortran STATIC fortran/molsymm.f90)
  target_link_libraries(molsymm_fortran PUBLIC molsymm)
  set_target_properties(molsymm_fortran PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_mod
  )
  target_include_directories(molsymm_fortran
    INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/fortran_mod
  )
  if(MOLSYMM_BUILD_TESTS)
    add_executable(molsymm_fortran_test tests/test_fortran.f90)
    target_link_libraries(molsymm_fortran_test PRIVATE molsymm_fortran)
    add_test(NAME molsymm_fortran_test COMMAND molsymm_fortran_test)
  endif()
endif()
