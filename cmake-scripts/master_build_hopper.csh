#!/bin/csh -f
# Master build script for hopper, last updated 9/28/2012. 
# build the code in the 4 ways currently supported and submits some test jobs
# there are fewer tests run here than for jaguar since the allocation amount is smaller
# this script may not execute successfully if there are old builds or module settings laying around
 
# (1) execute from the main seacism directory
# (2) set the next two commands 

setenv TEST_DIR "$SCRATCH/cism_tests"
setenv CODE_DIR "$HOME/seacism"
cd $CODE_DIR

echo 'set the pgi env'
module unload cmake netcdf python
module swap PrgEnv-gnu PrgEnv-pgi; module load cmake/2.8.7 python netcdf/4.1.2 subversion

# NEEDED AFTER A FRESH CHECKOUT
#echo 'bootstrap'
#./bootstrap

#SERIAL BUILD WITH AUTOCONF
echo 'make distclean'
make distclean
echo 'configure serial autoconf build'
./configure-scripts/xk6-config-serial >& conf_serial.out
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
echo 'configure pgi cmake build'
cd xe6-pgi
./hopper-pgi-cmake >& conf_cmake_pgi.out
echo 'make parallel pgi'
make -j 4 >& cmake_pgi_build.out

echo 'copy pgi parallel executable to test directory'
cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_pgi

# PARALLEL BUILD WITH CMAKE GNU
echo 'change to gnu env'
module unload cmake netcdf python
module swap PrgEnv-pgi PrgEnv-gnu; module load cmake/2.8.7 python netcdf

cd ../xe6-gnu
echo 'configure gnu cmake build'
./hopper-gnu-cmake >& conf_gnu.out
echo 'make parallel gnu'
make -j 4 >& cmake_gnu_build.out

echo 'copy gnu parallel executable to test directory'
cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_gnu

# execute tests on hopper
# this part has not yet been made ready for hopper 9/28/2012
# TODO the small jobs need to be combined into one ijob submission to get through the queue

# simplest case, runs all builds and on a range of small processor counts 
#cd $TEST_DIR/dome30
#qsub ijob

# large but not challenging case, to test large processor counts
#cd $TEST_DIR/dome500
#qsub ijob

# ISMIP test case A
#cd $TEST_DIR/ismip-hom-a
#qsub ijob

# ISMIP test case C
#cd $TEST_DIR/ismip-hom-c
#qsub ijob

# confined shelf to periodic BC
#cd $TEST_DIR/confined-shelf
#qsub ijob

# circular shelf to periodic BC
#cd $TEST_DIR/circular-shelf
#qsub ijob

# smaller GIS case to test realistic ice sheet configuration
#cd $TEST_DIR/gis_10km
#qsub ijob

# high resolution GIS case to test realistic ice sheet configuration and longer time series
#cd $TEST_DIR/gis_5km
#qsub ijob

