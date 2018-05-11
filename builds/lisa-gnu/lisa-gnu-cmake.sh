# Run this script by typing: source lisa-gnu-cmake
# After this script completes, type: make -j 8
# If rebuilding, type 'make clean' before running 'make -j 8'

# This cmake configuration script is set up to perform a parallel build

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#  source $CISM/builds/lisa-gnu/lisa-gnu-cmake $CISM
# where $CISM is the path to the top level cism directory.
if [ $# -eq 0 ]
then
    cism_top="../.." 
else
    cism_top=${1%/}
fi

echo CISM: "${cism_top}"

module purge
module load surfsara

# This results in the following module list
# 1) libgfortran/32/1(default)   5) oldwheezy/1.0(default)  
# 2) stdenv/1.3(default)         6) moab/default            
# 3) compilerwrappers            7) surfsara/1.1(default)   
# 4) licenses/1.0(default)

# Note: netcdf was compiled by hand to match the used compiler version with
# the same module set as above
##netCDF-Fortran 4.4.4:
# ./configure --prefix=/home/hgoelzer/prog/netcdf/build
# make install
##netCDF configure 4.5.0
# ./configure --prefix=/home/hgoelzer/prog/netcdf/build --disable-netcdf-4 --disable-dap
# make install

NETCDF_ROOT="/home/hgoelzer/prog/netcdf/build"

# Note: openmpi paths are specified manually as loading the corresponding module 
# leads to side effects.

MPI_INC_DIR="/sara/sw/openmpi-gnu-1.6.5-gf4.7/include"
MPI_LIB_DIR="/sara/sw/openmpi-gnu-1.6.5-gf4.7/lib"

# remove old build data:
rm -f ./CMakeCache.txt
rm -rf ./CMakeFiles

echo
echo "Doing CMake Configuration step"

cmake \
  -D CISM_BUILD_CISM_DRIVER:BOOL=ON \
  -D CISM_ENABLE_BISICLES=OFF \
  -D CISM_ENABLE_FELIX=OFF \
\
  -D CISM_USE_TRILINOS:BOOL=${CISM_USE_TRILINOS:=OFF} \
  -D CISM_MPI_MODE:BOOL=ON \
  -D CISM_SERIAL_MODE:BOOL=OFF \
\
  -D CISM_USE_GPTL_INSTRUMENTATION:BOOL=OFF \
  -D CISM_COUPLED:BOOL=OFF \
\
  -D CISM_GNU=ON \
\
  -D CISM_TRILINOS_DIR=$CISM_TRILINOS_DIR \
\
  -D CISM_NETCDF_DIR=$NETCDF_ROOT \
\
  -D CMAKE_Fortran_FLAGS="-g -O2 -ffree-line-length-none -fPIC -fno-range-check" \
\
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_Fortran_COMPILER=mpif90 \
\
  -D CISM_EXTRA_LIBS:STRING="-lblas" \
\
  -D CISM_MPI_INC_DIR=$MPI_INC_DIR \
  -D CISM_MPI_LIB_DIR=$MPI_LIB_DIR \
\
  -D CMAKE_VERBOSE_MAKEFILE=OFF \
  "${cism_top}"

# Note: last argument above  "../.." or ${cism_top} is path to top cism directory
