# cism_driver and CISM front-end build
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

link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}
                  ${GLIMMER_DYCORE_DIR}
		  ${BISICLES_INTERFACE_DIR}/${BISICLES_LIB_SUBDIR}
                  ${CISM_HDF5_LIB_DIR} ${PYTHON_LIB_DIR} )

# These are local source files needed to make the cism_driver executable and CISM front-end
add_executable(cism_driver cism_driver.F90 cism_front_end.F90 cism_cesm_interface.F90
               cism_internal_dycore_interface.F90 cism_external_dycore_interface.F90 
               simple_forcing.F90 testsfg.F90)

MESSAGE("GLIMMER_BISICLES_DIR: " ${GLIMMER_BISICLES_DIR})
MESSAGE("CISM_HDF5_LIBS: " ${CISM_HDF5_LIBS})

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)

IF (${CISM_BUILD_SIMPLE_BISICLES} OR ${CISM_BUILD_SIMPLE_FELIX})
  set(CISM_USE_EXTERNAL_DYCORE ON)
  ELSE()
  set(CISM_USE_EXTERNAL_DYCORE OFF)
ENDIF()
MESSAGE("CISM_USE_EXTERNAL_DYCORE: " ${CISM_USE_EXTERNAL_DYCORE})
IF (NOT ${CISM_USE_EXTERNAL_DYCORE})
  link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}
                  ${CISM_HDF5_LIB_DIR} ${PYTHON_LIB_DIR} )

  target_link_libraries(cism_driver
      glimmercismfortran
      ${PYTHON_LIBS}
      ${GLIMMER_NETCDF_LIBS}
      ${CISM_HDF5_LIBS}
      ${GLIMMER_MPI_LIBS}
      ${CISM_TRILINOS_LIBS}
      ${GLIMMER_EXTRA_LIBS}
  )
ELSEIF (${CISM_BUILD_SIMPLE_BISICLES})
MESSAGE("GLIMMER_DYCORE_DIR: " ${GLIMMER_DYCORE_DIR})
  link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}
                  ${GLIMMER_DYCORE_DIR}
		  ${BISICLES_INTERFACE_DIR}/${BISICLES_LIB_SUBDIR}
                  ${CISM_HDF5_LIB_DIR} ${PYTHON_LIB_DIR} )
  target_link_libraries(cism_driver
      ${GLIMMER_FMAIN}
      glimmercismfortran
      DyCoreToGlimmer
      libBisicles.a
      libChomboLibs.a
      ${PYTHON_LIBS}
      ${GLIMMER_NETCDF_LIBS}
      ${CISM_HDF5_LIBS}
      ${GLIMMER_MPI_LIBS}
      ${CISM_TRILINOS_LIBS}
      ${GLIMMER_EXTRA_LIBS}
      )
ENDIF()

# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable  simple_bisicles  should appear in dir: build/example-drivers/simple_bisicles/src")
MESSAGE("")
