# First create a library of all the sources that produce modules
# that are used in other files. These are forced to be built before
# the rest of the F90 files by the "add_dependencies" line.
set(MODSOURCES
    glimmer_config.F90
    glimmer_global.F90
    glimmer_log.F90
    parallel_mpi.F90
    glimmer_filenames.F90
    glimmer_coordinates.F90
    glimmer_paramets.F90
    glimmer_physcon.F90
    glimmer_ncdf.F90
    glimmer_utils.F90
    glimmer_scales.F90
    glimmer_map_types.F90
    glimmer_map_init.F90
    glimmer_map_CFproj.F90
    kelvin.F90
    profile.F90
   )
add_library(glimmermodules ${MODSOURCES})
# End creation of library of modules

# Now create libglimmer. It doesn't seem to hurt to GLOB all the
# F90 files, which redundantly grrabs all the enumerated MODSOURCES above.

# Get all F90 files (would be better to list explicitly then to "GLOB")
FILE(GLOB SOURCES *.F90 writestats.c)

# Then remove these three that really aren't part of the build.
LIST(REMOVE_ITEM SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/parallel_single.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/parallel_slap.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/unittest_glimmer_deriv.F90 )

# Make a glmmercmake library from these, after glimmermodules
add_library(glimmercmake ${SOURCES})
add_dependencies(glimmercmake glimmermodules)

# Specify include directories needed to find header files
# .mod files end up in BINARY tree, not source tree.
include_directories (../ ${GLIMMER_BINARY_DIR}/libglimmer ${GLIMMER_NETCDF_DIR}/include)
