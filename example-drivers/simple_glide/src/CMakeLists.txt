# Get libraries for link line from Trilinos build information
set(TRILINOS_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
   ${Trilinos_EXTRA_LD_FLAGS} ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )

# Need include directories from Trilinos but also mod files from glimmer
include_directories ( ./ ${GLIMMER_BINARY_DIR}/libglimmer
${GLIMMER_BINARY_DIR}/libglide ${GLIMMER_BINARY_DIR}/libglimmer-solve
${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

# These are local source files needed to make  sgcmake  execultable
add_executable(sgcmake  simple_glide.F90 simple_forcing.F90 testsfg.F90)

# Executable depends on several glimmer libraries and Trilinos,
# and potentially an f90main.o file (on jaguar with PGI at least)
target_link_libraries(sgcmake  
      ${GLIMMER_FMAIN}
      glidecmake glimmercmake glimmersolvecmake
      glimmersolveslapcmake glimmertrilinoscmake
      ${TRILINOS_LIBS}) 

# Helpful(?) message near end of configuration step
MESSAGE("\n   Execute 'make -j 8'")
MESSAGE("   Executable  sgcmake  should appear in dir: build/example-drivers/simple_glide/src")
MESSAGE("")
