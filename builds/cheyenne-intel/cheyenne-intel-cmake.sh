#!/bin/sh

# After this script completes, type: make -j 8
# If rebuilding, type 'make clean' before running 'make -j 8'

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#   $ source $CISM/builds/linux-gnu-cism/linux-gnu-cism-cmake $CISM
# where $CISM is the path to the top level cism directory.
if [ $# -eq 0 ]
then
    cism_top="../.." 
else
    cism_top=${1%/}
fi

echo CISM: "${cism_top}"


module purge
module load ncarenv/1.2
module load intel/17.0.1
module load mkl/2017.0.1
module load mpt/2.15f
module load netcdf-mpi/4.4.1.1
module load ncarcompilers/0.4.1
module load pnetcdf/1.8.0
module load cmake/3.7.2
module load python/2.7.13
module load numpy/1.12.0
module load netcdf4-python/1.2.7

# remove old build data:
rm -f ./CMakeCache.txt
rm -rf ./CMakeFiles

echo
echo "Doing CMake Configuration step"

# Note: the compilation flags were taken from the defaults for a CESM build on
# cheyenne-intel (using cime at 84aafd5). Some of these options are probably
# unnecessary for a standalone cism build, but I am keeping things consistent
# with the CESM build for simplicity.

# CISM_USE_GPTL_INSTRUMENTATION -- ON by default, set to OFF to not use GPTL instrumentation.

cmake \
  -D CISM_BUILD_CISM_DRIVER:BOOL=ON \
  -D CISM_ENABLE_BISICLES=OFF \
  -D CISM_ENABLE_FELIX=OFF \
\
  -D CISM_USE_TRILINOS:BOOL=OFF \
  -D CISM_MPI_MODE:BOOL=ON \
  -D CISM_SERIAL_MODE:BOOL=OFF \
\
  -D CISM_USE_GPTL_INSTRUMENTATION:BOOL="${CISM_USE_GPTL_INSTRUMENTATION:=ON}" \
  -D CISM_COUPLED:BOOL=OFF \
  -D CISM_USE_CISM_FRONT_END:BOOL=OFF \
\
  -D CISM_NETCDF_DIR="$NETCDF" \
  -D CISM_GPTL_DIR= "utils/libgptl" \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
\
  -D CMAKE_CXX_COMPILER=mpiicpc \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_Fortran_COMPILER=mpif90 \
\
  -D CMAKE_Fortran_FLAGS:STRING="-qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -qopt-report -O2 -debug minimal " \
  -D CMAKE_C_FLAGS:STRING="-qno-opt-dynamic-align -fp-model precise -std=gnu99 -qopt-report -O2 -debug minimal " \
  -D CMAKE_CXX_FLAGS:STRING="-qno-opt-dynamic-align -fp-model precise -std=gnu99 -qopt-report -O2 -debug minimal " \
  "${cism_top}"

