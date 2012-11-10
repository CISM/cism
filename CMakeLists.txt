# CMAKE File for GLIMMER building against an installed Trilinos

cmake_minimum_required(VERSION 2.8.6)

OPTION (NO_TRILINOS "Toggle to Not use Trilinos: defaults to OFF (Trilinos on by default)" OFF)
OPTION (CISM_MPI_MODE "Toggle to Configure with MPI: defaults to ON" ON)
OPTION (CISM_SERIAL_MODE "Toggle to Configure in Serial mode: defaults to OFF " OFF)
OPTION (CISM_BUILD_SIMPLE_GLIDE "Toggle to build simple_glide, ON by default)" ON)
OPTION (CISM_BUILD_SIMPLE_BISICLES "Toggle to build simple_bisicles, off by default)" OFF)
OPTION (CISM_BUILD_EXTRA_EXECUTABLES "Toggle to other executables, off by default)" OFF)
OPTION (CISM_USE_GPTL_INSTRUMENTATION "Toggle to use GPTL instrumentation, on by default " ON)



IF (NOT DEFINED(CISM_BINARY_DIR))
   SET(CISM_BINARY_DIR ${CMAKE_BINARY_DIR})
ENDIF()

MESSAGE("Building in: ${CISM_BINARY_DIR}")

IF (NOT ${NO_TRILINOS})
  OPTION(GLIMMER_TRILINOS_DIR "Required path to installed Trilinos")
  OPTION(GLIMMER_NETCDF_DIR "Required path to installed Netcdf")

  # Error check up front
  IF (NOT DEFINED GLIMMER_TRILINOS_DIR)
    MESSAGE(FATAL_ERROR "\nGLIMMER Error: cmake must define GLIMMER_TRILINOS_DIR:
    (-D GLIMMER_TRILINOS_DIR=<trilinos_install_prefix>)!")
  ENDIF()

  # Get Trilinos as one entity
  SET(CMAKE_PREFIX_PATH ${GLIMMER_TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
  FIND_PACKAGE(Trilinos REQUIRED)

  IF (${Trilinos_VERSION} VERSION_LESS  10.8.0)
    MESSAGE(FATAL_ERROR "Trilinos version 10.8 or newer required!")
  ENDIF()

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
  IF (NOT DEFINED CMAKE_CXX_COMPILER) 
    SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
  ENDIF()
  IF (NOT DEFINED CMAKE_C_COMPILER) 
    SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
  ENDIF()
  IF (NOT DEFINED CMAKE_Fortran_COMPILER)
    SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
  ENDIF()
ENDIF()

ENABLE_LANGUAGE(Fortran) 

MESSAGE(">> GLIMMER_NETCDF_DIR set to : ${GLIMMER_NETCDF_DIR}")

IF (NOT DEFINED GLIMMER_NETCDF_DIR)
  MESSAGE(FATAL_ERROR "\nGLIMMER Error: cmake must define GLIMMER_NETCDF_DIR:
  (-D GLIMMER_NETCDF_DIR=<netcdf_install_prefix>)!")
ENDIF()


FIND_PATH(CISM_NETCDFF_FOUND libnetcdff.a ${GLIMMER_NETCDF_DIR}/lib)

#MESSAGE(" CISM_NETCDFF_FOUND = ${CISM_NETCDFF_FOUND}")
IF (${CISM_NETCDFF_FOUND} STREQUAL ${GLIMMER_NETCDF_DIR}/lib )
   SET(GLIMMER_NETCDF_LIBS "netcdff;netcdf" CACHE STRING "Netcdf Library Names(s)")
ELSE()
   SET(GLIMMER_NETCDF_LIBS "netcdf" CACHE STRING "Netcdf Library Names(s)") 
ENDIF()
MESSAGE(">> GLIMMER_NETCDF_LIBS Library(s) set to : ${GLIMMER_NETCDF_LIBS}")

IF (DEFINED CISM_MPI_BASE_DIR)
  IF (NOT DEFINED CISM_MPI_LIB_DIR)
    SET(CISM_MPI_LIB_DIR ${CISM_MPI_BASE_DIR}/lib)
  ENDIF()
  IF (NOT DEFINED CISM_MPI_INC_DIR)
    SET(CISM_MPI_INC_DIR ${CISM_MPI_BASE_DIR}/include)
  ENDIF()
ENDIF()
MESSAGE(">> CISM_MPI_LIB_DIR set to : ${CISM_MPI_LIB_DIR}")
MESSAGE(">> CISM_MPI_INC_DIR set to : ${CISM_MPI_INC_DIR}")


SET(CMAKE_Fortran_MODULE_DIRECTORY ${CISM_BINARY_DIR}/fortran_mod_files)


PROJECT(GLIMMER)

MESSAGE(">> Calling cmake-scripts/autogenerate-script")
EXECUTE_PROCESS(COMMAND ./autogenerate-script
                WORKING_DIRECTORY ${GLIMMER_SOURCE_DIR}/cmake-scripts
                OUTPUT_FILE ${CISM_BINARY_DIR}/autogenerate.log)
MESSAGE(">> see ${CISM_BINARY_DIR}/autogenerate.log")

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

IF (NOT ${NO_TRILINOS}) 
   ## GLIMMER requires Trilinos to include the Piro build 
   MESSAGE("-- Looking for Piro within Trilinos:")
   MESSAGE("Trilinos pkg list: " ${Trilinos_PACKAGE_LIST})
   LIST(FIND Trilinos_PACKAGE_LIST  Piro Piro_List_ID)
   IF (Piro_List_ID GREATER -1)
     MESSAGE("-- Looking for Piro: -- found! Continuing.")
   ELSE()
     MESSAGE(FATAL_ERROR "-- Looking for Piro: -- not found! Rebuild Trilinos with Piro package.")
   ENDIF()

   ADD_DEFINITIONS(-DGLIMMER_MPI -DTRILINOS)
   MESSAGE(">> Trilinos build: adding  -DTRILINOS  to compile lines")

   include_directories(${GLIMMER_NETCDF_DIR}/include)
   link_directories(${GLIMMER_NETCDF_DIR}/lib) 

ELSE()

   IF (${CISM_MPI_MODE})
      ADD_DEFINITIONS(-DGLIMMER_MPI)
   ENDIF()
   include_directories(${GLIMMER_NETCDF_DIR}/include ${CISM_MPI_INC_DIR})
   link_directories(${GLIMMER_NETCDF_DIR}/lib ${CISM_MPI_LIB_DIR} ${CISM_SCI_LIB_DIR} 
  		 ${GLIMMER_DYCORE_DIR} ${CISM_HDF5_LIB_DIR} ${BISICLES_INTERFACE_DIR}) 

ENDIF()

IF (CISM_BUILD_SIMPLE_BISICLES)
  MESSAGE("\nBISICLES build configured.\n")  

  SET(GLIMMER_DYCORE_DIR ${GLIMMER_SOURCE_DIR}/libdycore)
  SET(GLIMMER_TO_DYCORE ${GLIMMER_DYCORE_DIR}/glimmer_to_dycore.F90) 
  SET(GLIMMER_BISICLES_DIR ${GLIMMER_DYCORE_DIR}/BISICLES)
  
  INCLUDE_DIRECTORIES(${GLIMMER_DYCORE_DIR})
ENDIF()

IF ( CISM_USE_GPTL_INSTRUMENTATION AND 
    (NOT CISM_BUILD_SIMPLE_BISICLES) AND
    (NOT CISM_SERIAL_MODE))
 OPTION(GLIMMER_GPTL_DIR "Optional path to installed gptl library")
 IF (GLIMMER_GPTL_DIR)
   message(">> GPTL Library reqested, installed in ${GLIMMER_GPTL_DIR}")
   ADD_DEFINITIONS(-DCESMTIMERS)
   include_directories(${GLIMMER_GPTL_DIR})
   link_directories(${GLIMMER_GPTL_DIR})
   SET(GLIMMER_GPTL_LIB gptl)
 ELSE()
   message(">>GPTL Library not requested: can set GLIMMER_GPTL_DIR to enable")
 ENDIF()
ENDIF()

#####Build All Fortran Sources #####
file(GLOB FORTRANSOURCES
     libglimmer-solve/SLAP/*.f
     libglimmer-solve/*.F90
     libglimmer/*.F90 libglimmer/writestats.c
     libglide/*.F90 libglide/fgmresD.f90 
     ${GLIMMER_TO_DYCORE})

LIST(REMOVE_ITEM FORTRANSOURCES
   ${GLIMMER_SOURCE_DIR}/libglimmer/unittest_glimmer_deriv.F90 
   ${GLIMMER_SOURCE_DIR}/libglide/relaxed.F90
   ${GLIMMER_SOURCE_DIR}/libglide/test_setup.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_commandline.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_config.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_filenames.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_integrate.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_ts.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_writestats.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/nc2config.F90 )

IF ((NOT ${NO_TRILINOS}) OR CISM_MPI_MODE)
   LIST(REMOVE_ITEM FORTRANSOURCES
      ${GLIMMER_SOURCE_DIR}/libglimmer/parallel_slap.F90   
      ${GLIMMER_SOURCE_DIR}/libglimmer/glimmer_horiz_bcs_serial.F90)

ELSE(${CISM_SERIAL_MODE})
   LIST(REMOVE_ITEM FORTRANSOURCES
      ${GLIMMER_SOURCE_DIR}/libglimmer/parallel_mpi.F90
      ${GLIMMER_SOURCE_DIR}/libglimmer/glimmer_horiz_bcs_parallel.F90)
ENDIF()

### (For CESM) Remove source files with names already in GLIMMER_SOURCEMOD_DIR
OPTION(GLIMMER_SOURCEMOD_DIR 
      "Path to SourceMod directory of F90 files to replace Glimmer files")

#MESSAGE("Fortran Source Files: ${FORTRANSOURCES}")

FILE(GLOB SOURCEMODFILES ${GLIMMER_SOURCEMOD_DIR}/*.F90)

# MESSAGE("glimmer_sourcemod_dir: " ${GLIMMER_SOURCEMOD_DIR})
# MESSAGE("Fortran Mod Files: ${SOURCEMODFILES}")

FOREACH( MODFILE ${SOURCEMODFILES})
  STRING(FIND ${MODFILE} / index REVERSE)
  MATH(EXPR index ${index}+1)
  STRING(SUBSTRING ${MODFILE} ${index} -1 filename)
  FOREACH( SOURCEFILE ${FORTRANSOURCES})
    STRING(REGEX MATCH ${filename} match_found ${SOURCEFILE})

    IF(match_found)
      MESSAGE("--SourceMod: removing ${SOURCEFILE} in favor of ${MODFILE}")
      LIST(REMOVE_ITEM FORTRANSOURCES ${SOURCEFILE})
   ENDIF()
  ENDFOREACH()
ENDFOREACH()
###


SET(CMAKE_INCLUDE_PATH ${CISM_BINARY_DIR}/include)


#message("Fortransource ${FORTRANSOURCES}")

# include-dirs: Binary: all the .mod files; SOURCE: config.inc; 
#               libglide: glide_mask.inc

include_directories (${CISM_BINARY_DIR} ${GLIMMER_SOURCE_DIR}
                     ${GLIMMER_SOURCE_DIR}/libglide)


SET(LIBRARY_OUTPUT_PATH ${CISM_BINARY_DIR}/lib)

add_library(glimmercismfortran ${FORTRANSOURCES})


#####Build C++ Sources #####

IF (NOT ${NO_TRILINOS})
  add_subdirectory(libglimmer-trilinos)
ENDIF()

#####Build simple_glide executable (unless turned off) #####
OPTION(GLIMMER_NO_EXECUTABLE "Set to  ON  to just build libraries (default:OFF)"  OFF)


IF (CISM_BUILD_SIMPLE_GLIDE AND (NOT GLIMMER_NO_EXECUTABLE))
  add_subdirectory(example-drivers/simple_glide/src)
ENDIF()


IF (CISM_BUILD_SIMPLE_BISICLES)
  add_subdirectory(example-drivers/simple_bisicles/src)
  add_subdirectory(libdycore)
  SET(GLIMMER_EXTERNAL_DYCORE_DIR ${GLIMMER_SOURCE_DIR}/libdycore)
  SET(GLIMMER_BISICLES_DIR ${GLIMMER_EXTERNAL_DYCORE_DIR}/BISICLES)
  include_directories (${GLIMMER_EXTERNAL_DYCORE_DIR}
                       ${GLIMMER_BISICLES_DIR})

  message("glimmer src dir: ${GLIMMER_SOURCE_DIR}")
  message("glimmer extern dycore dir: ${GLIMMER_EXTERNAL_DYCORE_DIR}")
  get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
  message("inc_dirs = ${inc_dirs}")

ENDIF()


IF(CISM_BUILD_EXTRA_EXECUTABLES)
  add_subdirectory(example-drivers/eis/src)
  message("Building extra executables: ")
ENDIF()

#COMMAND (environment)


