#IK, 8/6/13: modified this CMakeLists.txt for simple_felix 
# simple_felix build
# Get libraries for link line from Trilinos build information
# 

IF (NOT ${NO_TRILINOS})
  set(CISM_TRILINOS_LIBS glimmercismcpp ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
     ${Trilinos_EXTRA_LD_FLAGS} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
ELSE()
  set(CISM_TRILINOS_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )
ENDIF()

#message("")
#message("  CMake detected the following libraries for linking Fortran with C++ compiler:")
#message("     ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ")

# Need include directories from Trilinos but also mod files from glimmer
include_directories (${GLIMMER_BINARY_DIR}/include  ${PYTHON_INC_DIR}
                     ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

#IK, 8/6/13: the following line ensures Albany Felix and glimmerInterface
#libraries get linked.
#Add error checking if not found ${GLIMMER_ALBANY_DIR}/FELIX or
#${GLIMMER_ALBANY_DIR}/glimmerInterface
link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}
                  ${GLIMMER_DYCORE_DIR} ${GLIMMER_ALBANY_DIR}
                  ${CISM_HDF5_LIB_DIR} ${PYTHON_LIB_DIR} )

# These are local source files needed to make  simple_felix  executable
add_executable(simple_felix  simple_felix.F90 simple_forcing.F90 testsfg.F90)

MESSAGE("GLIMMER_ALBANY_DIR: " ${GLIMMER_ALBANY_DIR})
MESSAGE("CISM_HDF5_LIBS: " ${CISM_HDF5_LIBS})

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)
target_link_libraries(simple_felix
      ${GLIMMER_FMAIN}
      glimmercismfortran
      DyCoreToGlimmer
      libFelixCISMInterface.a
      libFELIX.a
      ${PYTHON_LIBS}
      ${GLIMMER_NETCDF_LIBS}
      ${CISM_HDF5_LIBS}
      ${GLIMMER_MPI_LIBS}
      ${CISM_TRILINOS_LIBS}
      ${GLIMMER_EXTRA_LIBS}
      )



# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable  simple_felix  should appear in dir: build/example-drivers/simple_felix/src")
MESSAGE("")
