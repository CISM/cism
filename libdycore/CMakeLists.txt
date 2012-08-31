# libdycore build

FILE(GLOB BISICLES_SOURCES *.cpp BISICLES/*.cpp)
FILE(GLOB BISICLES_HEADERS *.H BISICLES/*.H)

add_library(DyCoreToGlimmer ${BISICLES_SOURCES} ${BISICLES_HEADERS})

# Need include directories from Trilinos but also mod files from glimmer

include_directories (${BISICLES_INTERFACE_DIR})

MESSAGE("CISM_BISICLES_DIR: " ${CISM_BISICLES_DIR})


