cmake_minimum_required(VERSION 3.16)
project(orion_openfhe CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Point CMAKE_PREFIX_PATH (or OpenFHE_DIR) at an OpenFHE installation.
# An OpenFHE built with the Intel HEXL extension needs the HEXL package as well.
find_package(HEXL QUIET)
find_package(OpenFHE REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenFHE_CXX_FLAGS}")

add_library(orion_openfhe SHARED src/orion_openfhe.cpp)
# OpenFHE's headers only enable their OpenMP code paths (and thread limits) with PARALLEL.
target_compile_definitions(orion_openfhe PRIVATE PARALLEL)
target_include_directories(orion_openfhe PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/../c_api
  ${OPENMP_INCLUDES}
  ${OpenFHE_INCLUDE}
  ${OpenFHE_INCLUDE}/third-party/include
  ${OpenFHE_INCLUDE}/core
  ${OpenFHE_INCLUDE}/pke
  ${OpenFHE_INCLUDE}/binfhe)
target_link_directories(orion_openfhe PRIVATE ${OpenFHE_LIBDIR} ${OPENMP_LIBRARIES})
target_link_libraries(orion_openfhe PRIVATE ${OpenFHE_SHARED_LIBRARIES})
set_target_properties(orion_openfhe PROPERTIES
  OUTPUT_NAME "openfhe-linux" PREFIX "" SUFFIX ".so"
  BUILD_RPATH "${OpenFHE_LIBDIR}" INSTALL_RPATH "${OpenFHE_LIBDIR}")
