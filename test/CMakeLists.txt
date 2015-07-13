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
    main.cpp
)

add_executable(polytopetest ${SOURCES})

TARGET_LINK_LIBRARIES(polytopetest polytope)
#~ TARGET_LINK_LIBRARIES(polytopetest polytope ${SRC_DIR}/../external/cddlib-094b/lib-src/libcdd.a)
