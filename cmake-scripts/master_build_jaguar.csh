#!/bin/csh -f

# user needs to set these two commands
setenv TEST_DIR "/tmp/work/4ue/SEACISM/cism_tests"
setenv CODE_DIR "$HOME/seacism"
cd $CODE_DIR

# NEEDED AFTER A FRESH CHECKOUT
echo 'bootstrap'
./bootstrap

module load netcdf-hdf5parallel/4.2.0 python subversion

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
./configure-scripts/xk6-config >& conf_auto_pgi.out
echo 'make parallel'
make >& auto_pgi_build.out 

echo 'copy autconf parallel executable to test directory'
cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide

# PARALLEL BUILD WITH CMAKE PGI
echo 'configure pgi cmake build'
cd xk6-pgi
./jaguar-pgi-cmake >& conf_cmake_pgi.out
echo 'make parallel pgi'
make -j 8 >& cmake_pgi_build.out

module load cmake
echo 'copy pgi parallel executable to test directory'
cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_pgi

# PARALLEL BUILD WITH CMAKE GNU
echo 'change to gnu env'
module swap PrgEnv-pgi PrgEnv-gnu; module load cmake python netcdf-hdf5parallel/4.2.0

cd ../xk6-gnu
echo 'configure gnu cmake build'
./jaguar-gnu-cmake >& conf_gnu.out
echo 'make parallel gnu'
make -j 8 >& cmake_gnu_build.out

echo 'copy gnu parallel executable to test directory'
cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_gnu

# execute tests on jaguar
# TODO the small jobs need to be combined into one ijob submission to get through the queue

# simplest case, runs all builds and on a range of small processor counts 
cd $TEST_DIR/dome30
qsub ijob

# large but not challenging case, to test large processor counts
cd $TEST_DIR/dome500
qsub ijob

# ISMIP test case A
cd $TEST_DIR/ismip-hom-a
qsub ijob

# ISMIP test case C
#cd $TEST_DIR/ismip-hom-c
#qsub ijob

# confined shelf to periodic BC
cd $TEST_DIR/confined-shelf
qsub ijob

# circular shelf to periodic BC
cd $TEST_DIR/circular-shelf
qsub ijob

# smaller GIS case to test realistic ice sheet configuration
cd $TEST_DIR/gis_10km
qsub ijob

# high resolution GIS case to test realistic ice sheet configuration and longer time series
cd $TEST_DIR/gis_5km
qsub ijob

