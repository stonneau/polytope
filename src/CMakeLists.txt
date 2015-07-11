cmake_minimum_required(VERSION 2.6)

include_directories("${SRC_DIR})
include_directories("${SRC_TOOLS_DIR})
include_directories("${EIGEN3_INCLUDE_DIR}")
include_directories("${CDD_INCLUDE_DIR}")

PROJECT(cddwrap)


if ( MSVC )
	SET(CMAKE_DEBUG_POSTFIX d)
endif ( MSVC )
	
set(SOURCES
    DynamicStability.h
    DynamicStability.cpp
    MatrixDefs.h MatrixDefs.cpp
    MatrixDefsInternal.h
)


LINK_DIRECTORIES(${CDD_LIBRARY_DIR})

#add_library(cddwrap ${SOURCES})

#TARGET_LINK_LIBRARIES(cddwrap ${CDD_LIBRARY})

#SET_TARGET_PROPERTIES(cddwrap PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
#SET_TARGET_PROPERTIES(cddwrap PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")


add_executable(cddwrap ${SOURCES})

TARGET_LINK_LIBRARIES(cddwrap ${CDD_LIBRARY})
