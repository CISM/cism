# Run this script by typing: source cartesius-intel-serial-cmake
# After this script completes, type: make 
# If rebuilding, type 'make clean' before running 'make'

# This cmake configuration script is set up to perform a serial build

# Set path to top cism directory
# Note, this is an easy way to build out of source.
# In directory you want to build in, run:
#  source $CISM/builds/cartesius-intel/cartesius-intel-serial-cmake $CISM
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
module load hdf5/serial/intel/1.8.10-patch1
module load netcdf
module load mkl
module load cmake
module load lapack/mkl
module load blas/netlib/intel

## Resulting module list 
# 1) compilerwrappers                7) hdf5/serial/intel/1.8.10-patch1        
# 2) c/intel/15.0.0(default)         8) netcdf/serial/intel/4.1.3(default)     
# 3) fortran/intel/15.0.0(default)   9) mkl/11.2(default)                      
# 4) mpi/impi/5.0.3.048(default)    10) cmake/3.5.2(default)                   
# 5) bull                           11) lapack/mkl(default)                    
# 6) surfsara                       12) blas/netlib/intel/2011.04.19(default)  

## where 'module load surfsara' will give you the following 
# 1) compilerwrappers          3) fortran/intel/15.0.0(default)   5) bull      
# 2) c/intel/15.0.0(default)   4) mpi/impi/5.0.3.048(default)     6) surfsara  

# remove old build data:
rm -f ./CMakeCache.txt
rm -fr ./CMakeFiles

echo
echo "Doing CMake Configuration step"

# A few non-intuitive things:
#
# - CISM_FORCE_FORTRAN_LINKER: without this, cmake tries to use a C++ linker, which doesn't work
#
# - CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES: (this is a note that applies to the
#   parallel build with trilinos, and may or may not apply to this serial
#   build): if this is on (the default), some libraries are included on the link
#   line which can't be found (e.g., hdf5). This may be related to the fact that
#   trilinos on yellowstone is old, and/or the fact that cmake wants to use a
#   C++ linker but we're telling it to use a fortran linker.

cmake \
  -D CISM_USE_TRILINOS:BOOL=OFF \
  -D CISM_COUPLED:BOOL=OFF \
  -D CISM_MPI_MODE:BOOL=OFF \
  -D CISM_SERIAL_MODE:BOOL=ON \
  -D CISM_BUILD_SIMPLE_GLIDE:BOOL=ON \
  -D CISM_BUILD_SIMPLE_BISICLES:BOOL=OFF \
  -D CISM_BUILD_GLINT_EXAMPLE:BOOL=OFF \
  -D CISM_BUILD_CISM_DRIVER:BOOL=ON \
  -D CISM_USE_GPTL_INSTRUMENTATION:BOOL=OFF \
  -D CISM_USE_DEFAULT_IO:BOOL=OFF \
  -D CISM_USE_CISM_FRONT_END:BOOL=OFF \
\
  -D CISM_NETCDF_DIR=/hpc/sw/netcdf-fortran-4.2-intel-seq \
  -D CISM_FORCE_FORTRAN_LINKER:BOOL=ON \
  -D CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES:BOOL=ON \
  -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
\
  -D CMAKE_CXX_COMPILER=icpc \
  -D CMAKE_C_COMPILER=icc \
  -D CMAKE_Fortran_COMPILER=ifort \
\
  -D CMAKE_Fortran_FLAGS:STRING="-fp-model source -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -xHost -O2 -L/  /hpc/sw/blas-2011.04.19-intel/lib/libblas.a" \
  -D CMAKE_C_FLAGS:STRING="-O2 -fp-model precise -xHost" \
  -D CMAKE_CXX_FLAGS:STRING="-O2 -fp-model precise -xHost" \
  "${cism_top}"

# Note: last argument above  "../.." or ${cism_top} is path to top cism directory
