cmake_minimum_required(VERSION 2.6)

include_directories("${SRC_DIR}")
include_directories("${INCLUDE_DIR}")
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${CDD_INCLUDE_DIR}")

PROJECT(polytope)


if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )
	
set(SOURCES
    ${INCLUDE_DIR}/polytope/config.h
    ${INCLUDE_DIR}/polytope/stability_margin.h
    stability_margin.cpp
)

add_library(polytope SHARED ${SOURCES})

TARGET_LINK_LIBRARIES(polytope ${CDD_LIBRARY})

SET_TARGET_PROPERTIES(polytope PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
SET_TARGET_PROPERTIES(polytope PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
