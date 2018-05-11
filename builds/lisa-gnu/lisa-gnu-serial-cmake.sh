# Run this script by typing: source lisa-gnu-serial-cmake
# After this script completes, type: make 
# If rebuilding, type 'make clean' before running 'make'

# This cmake configuration script is set up to perform a serial build

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#  source $CISM/builds/lisa-gnu/lisa-gnu-serial-cmake $CISM
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

# Note: netcdf was compiled by hand to match the used compiler version with
# the same module set as above
##netCDF-Fortran 4.4.4:
# ./configure --prefix=/home/hgoelzer/prog/netcdf/build
# make install
##netCDF configure 4.5.0
# ./configure --prefix=/home/hgoelzer/prog/netcdf/build --disable-netcdf-4 --disable-dap
# make install

NETCDF_ROOT="/home/hgoelzer/prog/netcdf/build"

# remove old build data:
rm -f ./CMakeCache.txt
rm -rf ./CMakeFiles

echo
echo "Doing CMake Configuration step"

cmake \
  -D CISM_USE_TRILINOS:BOOL=OFF \
  -D CISM_COUPLED:BOOL=OFF \
  -D CISM_MPI_MODE:BOOL=OFF \
  -D CISM_SERIAL_MODE:BOOL=ON \
  -D CISM_BUILD_CISM_DRIVER=ON \
  -D CISM_BUILD_SIMPLE_GLIDE:BOOL=OFF \
  -D CISM_BUILD_SIMPLE_BISICLES:BOOL=OFF \
  -D CISM_BUILD_GLINT_EXAMPLE:BOOL=OFF \
  -D CISM_USE_GPTL_INSTRUMENTATION:BOOL=OFF \
  -D CISM_USE_DEFAULT_IO:BOOL=OFF \
  -D CISM_USE_CISM_FRONT_END:BOOL=ON \
\
  -D CISM_GNU=ON \
\
  -D CISM_NETCDF_DIR=$NETCDF_ROOT \
\
  -D CMAKE_Fortran_FLAGS="-g -O2 -ffree-line-length-none -fPIC -fno-range-check" \
\
  -D CMAKE_CXX_COMPILER=g++ \
  -D CMAKE_C_COMPILER=gcc \
  -D CMAKE_Fortran_COMPILER=gfortran \
\
  -D CISM_EXTRA_LIBS:STRING="-lblas" \
\
  -D CMAKE_VERBOSE_MAKEFILE=ON \
  "${cism_top}"

# Note: last argument above  "../.." or ${cism_top} is path to top cism directory

