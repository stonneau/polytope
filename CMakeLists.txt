cmake_minimum_required(VERSION 2.6)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/cpack.cmake)

if ( UNIX )  
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")
endif ()

SET(PROJECT_NAME polytope)
SET(PROJECT_DESCRIPTION
   "Stability margin implementation with polytopes"
 )
 SET(PROJECT_URL "")


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${PROJECT_SOURCE_DIR}/cmake2/")
set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin/")

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

set(EXTERNAL_LIBRARY_DIR "${PROJECT_SOURCE_DIR}/external/lib")

if ( MSVC )  
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS")
endif()

SETUP_PROJECT()

# Declare Headers
SET(${PROJECT_NAME}_HEADERS
    include/polytope/config.h
    include/polytope/stability_margin.h
  )

find_package(Eigen3)
if (EIGEN3_INCLUDE_DIR)
        message(STATUS "Found Eigen3")
else()
        message(STATUS "Eigen3 not found looking for Eigen2")
        find_package(Eigen2 REQUIRED)
endif()

find_package(CDD REQUIRED)

add_subdirectory (src/polytope)
add_subdirectory (test)

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()
