#!/bin/sh -f
# Master build script for hopper, last updated 10/25/2012 with v1583 
# build the code in the 4 ways currently supported and submits some test jobs
# there are fewer tests run here than for jaguar since the allocation amount is smaller
# this script may not execute successfully if there are old builds or module settings laying around
 
# (1) execute from the main seacism directory
# (2) set the next two commands 

#add logic at the top to decide which versions to build 
setenv TEST_DIR "$SCRATCH/cism_tests"
setenv CODE_DIR "$HOME/seacism"
cd $CODE_DIR
#mkdir -pv $TEST_DIR/configure_output

# even if these are set in your env you need these when running the script
echo 'set the pgi env'
module unload cmake netcdf python
module swap PrgEnv-gnu PrgEnv-pgi; module load cmake/2.8.7 python netcdf-hdf5parallel/4.2.0 subversion

# NEEDED AFTER A FRESH CHECKOUT
#echo 'bootstrap'
#./bootstrap

#SERIAL BUILD WITH AUTOCONF
echo 'make distclean'
make distclean 
echo 'configure serial autoconf build'
./configure-scripts/hopper-config-serial >& conf_serial.out
echo 'make serial'
make >& serial_build.out 

echo 'copy serial executable to test directory'
cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_serial

# PARALLEL BUILD WITH AUTOCONF PGI
echo 'make distclean'
make distclean
echo 'configure parallel autoconf build'
./configure-scripts/hopper-config >& conf_auto_pgi.out
echo 'make parallel'
make >& auto_pgi_build.out 

echo 'copy autconf parallel executable to test directory'
cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide

# set up directories for building cmake executables
rm -rf xe6-pgi
rm -rf xe6-gnu
mkdir xe6-pgi
mkdir xe6-gnu
cp cmake-scripts/hopper-pgi-cmake xe6-pgi
cp cmake-scripts/hopper-gnu-cmake xe6-gnu

# PARALLEL BUILD WITH CMAKE PGI
cd $CODE_DIR/xe6-pgi
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure pgi cmake build'
./hopper-pgi-cmake >& conf_cmake_pgi.out
echo 'make parallel pgi'
make -j 4 >& cmake_pgi_build.out

echo 'copy pgi parallel executable to test directory'
cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_pgi

# PARALLEL BUILD WITH CMAKE GNU
echo 'change to gnu env'
module unload cmake netcdf-hdf5parallel/4.2.0 python
module swap PrgEnv-pgi PrgEnv-gnu; module load cmake/2.8.7 python netcdf-hdf5parallel/4.2.0

cd $CODE_DIR/xe6-gnu
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure gnu cmake build'
./hopper-gnu-cmake >& conf_gnu.out
echo 'make parallel gnu'
make -j 4 >& cmake_gnu_build.out

echo 'copy gnu parallel executable to test directory'
cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_gnu

# execute tests on hopper
# TODO the small jobs need to be combined into one hopjob submission to get through the queue

# simplest case, runs all builds and on a range of small processor counts 
cd $TEST_DIR/dome30
qsub hopjob

# large but not challenging case, to test large processor counts, not yet configured for hopper
#cd $TEST_DIR/dome500
#qsub hopjob

# ISMIP test case A - not activated until BC set
#cd $TEST_DIR/ismip-hom-a
#qsub hopjob

# ISMIP test case C - not activated until BC set
#cd $TEST_DIR/ismip-hom-c
#qsub hopjob

# confined shelf to periodic BC
cd $TEST_DIR/confined-shelf
qsub hopjob

# circular shelf to periodic BC
cd $TEST_DIR/circular-shelf
qsub hopjob

# smaller GIS case to test realistic ice sheet configuration
cd $TEST_DIR/gis_10km
qsub hopjob

# high resolution GIS case to test realistic ice sheet configuration and longer time series, no yet configured for hopper
#cd $TEST_DIR/gis_5km
#qsub hopjob

