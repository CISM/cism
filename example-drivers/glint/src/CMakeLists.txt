# glint_example build
# Get libraries for link line from Trilinos build information
# 

IF (${CISM_USE_TRILINOS})
  set(CISM_TRILINOS_LIBS glimmercismcpp ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
     ${Trilinos_EXTRA_LD_FLAGS} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
ELSE()
  set(CISM_TRILINOS_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )
ENDIF()

#message("")
#message("  CMake detected the following libraries for linking Fortran with C++ compiler:")
#message("     ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ")

# Need include directories from Trilinos but also mod files from glimmer
include_directories (${GLIMMER_BINARY_DIR}/include)

# These are local source files needed to make glint_example  executable
add_executable(glint_example glint_example.F90 
               ${GLIMMER_BINARY_DIR}/fortran_autogen_srcs/glint_io.F90)

IF (CISM_STATIC_LINKING)
  SET_TARGET_PROPERTIES(glint_example PROPERTIES LINK_SEARCH_START_STATIC 1)
  SET_TARGET_PROPERTIES(glint_example PROPERTIES LINK_SEARCH_END_STATIC 1)
ENDIF()

IF (CISM_FORCE_FORTRAN_LINKER)
  SET_PROPERTY(TARGET glint_example PROPERTY LINKER_LANGUAGE Fortran)
ENDIF()

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)
target_link_libraries(glint_example
      ${GLIMMER_FMAIN}
      glimmercismfortran
      ${GLIMMER_NETCDF_LIBS}
      ${GLIMMER_GPTL_LIB}
      ${CISM_TRILINOS_LIBS}
      ${GLIMMER_EXTRA_LIBS})

# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable glint_example should appear in dir: build/example-drivers/glint/src")
MESSAGE("")
