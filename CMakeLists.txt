# CMAKE File for GLIMMER building against an installed Trilinos

cmake_minimum_required(VERSION 2.7)

# Error check up front
IF (NOT DEFINED GLIMMER_TRILINOS_DIR)
  MESSAGE(FATAL_ERROR "\nGLIMMER Error: cmake must define GLIMMER_TRILINOS_DIR:
  (-D GLIMMER_TRILINOS_DIR:FILEPATH=<trilinos_install_prefix>)!")
ENDIF()

# Get Trilinos as one entity
SET(CMAKE_PREFIX_PATH ${GLIMMER_TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)

MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
MESSAGE("End of Trilinos details\n")

SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

ENABLE_LANGUAGE(Fortran) 

#IF (NOT DEFINED GLIMMER_NETCDF_DIR)
#  MESSAGE(FATAL_ERROR "\nGLIMMER Error: cmake must define GLIMMER_NETCDF_DIR:
#  (-D GLIMMER_NETCDF_DIR=<netcdf_install_prefix>)!")
#ENDIF()


PROJECT(GLIMMER)


# End of setup and error checking
#  NOTE: PROJECT command checks for compilers, so this statement
#        is moved AFTER setting CMAKE_CXX_COMPILER from Trilinos

## Use CMAKE_CXX_FLAGS   CMAKE_Fortran_FLAGS   to override  Trilinos flags
## USe GLIMMER_CXX_FLAGS GLIMMER_Fortran_FLAGS to append to Trilinos flags

IF (NOT CMAKE_CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS} )
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLIMMER_CXX_FLAGS}")

IF (NOT CMAKE_Fortran_FLAGS)
  SET(CMAKE_Fortran_FLAGS ${Trilinos_Fortran_COMPILER_FLAGS} )
ENDIF()
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GLIMMER_Fortran_FLAGS}")

## GLIMMER requires Trilinos to include the Piro build 
MESSAGE("-- Looking for Piro within Trilinos:")
LIST(FIND Trilinos_PACKAGE_LIST  Piro Piro_List_ID)
IF (Piro_List_ID GREATER -1)
  MESSAGE("-- Looking for Piro: -- found! Continuing.")
ELSE()
  MESSAGE(FATAL_ERROR "-- Looking for Piro: -- not found! Rebuild Trilinos with Piro package.")
ENDIF()

######

include_directories(${GLIMMER_NETCDF_DIR}/include)
link_directories(${GLIMMER_NETCDF_DIR}/lib)

ADD_DEFINITIONS(-DGLIMMER_MPI -DTRILINOS)

add_subdirectory(libglimmer-trilinos)
add_subdirectory(libglimmer-solve)
add_subdirectory(libglimmer)
add_subdirectory(libglide)

add_subdirectory(example-drivers)
