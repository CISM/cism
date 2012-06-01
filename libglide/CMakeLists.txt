# Get all F90 files 
FILE(GLOB SOURCES *.F90 fgmresD.f90)
# Then remove this unused file
LIST(REMOVE_ITEM SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/relaxed.F90)

# Make a library from these, forcing it to be after the modules are built
add_library(glidecmake ${SOURCES})
add_dependencies(glidecmake glimmermodules)

# Specify include directories needed to find header files, where
# mod files are in BINARY tree, not source tree
include_directories (../ ${GLIMMER_BINARY_DIR}/libglimmer
    ${GLIMMER_BINARY_DIR}/libglimmer-solve ../libglimmer-trilinos)
