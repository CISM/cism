# libdycore build

FILE(GLOB DYCORE_SOURCES *.cpp *.cpp)
FILE(GLOB DYCORE_HEADERS *.H *.H)

#IK, 8/6/13: added logic to add BISICLES/* files only if 
#simple_bisicles is being built.  This is necessary 
#if there is the option of other external dycores besides bisicles.
IF (CISM_BUILD_SIMPLE_BISICLES) 
  FILE(GLOB BISICLES_SOURCES *.cpp BISICLES/*.cpp)
  FILE(GLOB BISICLES_HEADERS *.H BISICLES/*.H)
ENDIF()
# IK, 8/6/13: added FELIX directory sources/headers
IF (CISM_BUILD_SIMPLE_FELIX) 
  FILE(GLOB FELIX_SOURCES *.cpp FELIX/*.cpp)
  FILE(GLOB FELIX_HEADERS *.H FELIX/*.H)
ENDIF() 

add_library(DyCoreToGlimmer ${DYCORE_SOURCES} ${DYCORE_HEADERS} 
                            ${BISICLES_SOURCES} ${BISICLES_HEADERS}
                            ${FELIX_SOURCES} ${FELIX_HEADERS})

# Need include directories from Trilinos but also mod files from glimmer

#IK, 8/6/13: added include dirs from Albany/FELIX glimmerInterface source code 
#This is where felix_driver.cpp/.H are located.  
include_directories (${BISICLES_INTERFACE_DIR} ${FELIX_INTERFACE_DIR})

#IK, 8/6/13: is CISM_BISICLES_DIR set somewhere, or remnant from earlier
#version of code? 
MESSAGE("CISM_BISICLES_DIR: " ${CISM_BISICLES_DIR})


