# simple_bisicles build

IF (${CISM_USE_TRILINOS})
  LIST(INSERT CISM_TRILINOS_LIBS 0 glimmercismcpp)
ENDIF()

# Need include directories from Trilinos but also mod files from glimmer
include_directories (${GLIMMER_BINARY_DIR} ${PYTHON_INC_DIR}
                     ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS}
                  ${GLIMMER_DYCORE_DIR}
		  ${BISICLES_INTERFACE_DIR}/${BISICLES_LIB_SUBDIR}
                  ${CISM_HDF5_LIB_DIR} ${PYTHON_LIB_DIR} )

# These are local source files needed to make  simple_bisicles  executable
add_executable(simple_bisicles  simple_bisicles.F90
               ../../simple_glide/src/simple_forcing.F90
               ../../simple_glide/src/testsfg.F90)

IF (CISM_STATIC_LINKING)
  SET_TARGET_PROPERTIES(simple_bisicles PROPERTIES LINK_SEARCH_START_STATIC 1)
  SET_TARGET_PROPERTIES(simple_bisicles PROPERTIES LINK_SEARCH_END_STATIC 1)
ENDIF()

IF (CISM_FORCE_FORTRAN_LINKER)
  SET_PROPERTY(TARGET simple_bisicles PROPERTY LINKER_LANGUAGE Fortran)
ENDIF()

MESSAGE("GLIMMER_BISICLES_DIR: " ${GLIMMER_BISICLES_DIR})
MESSAGE("CISM_HDF5_LIBS: " ${CISM_HDF5_LIBS})

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)
target_link_libraries(simple_bisicles
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

# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable  simple_bisicles  should appear in dir: build/example-drivers/simple_bisicles/src")
MESSAGE("")
