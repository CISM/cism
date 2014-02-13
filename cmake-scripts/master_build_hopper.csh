#!/bin/csh
# Master build script for hopper, last updated 11/09/2012 with v1617 
# build the code in the 4 ways currently supported and submits some test jobs
# there are fewer tests run here than for jaguar since the allocation amount is smaller
 
# (1) execute from the main seacism directory
# (2) set the next two commands 

#add logic at the top to decide which versions to build 
setenv TEST_DIR "$GSCRATCH/higher-order"
setenv CODE_DIR "$HOME/PISCEES/trunk"
cd $CODE_DIR
# setting to 0 mens don't build that version
setenv build_autoconf 1
setenv build_cmake 1
# flags set for regression and performance suites
setenv REG_TEST 1
setenv PERF_TEST 0

#mkdir -pv $TEST_DIR/configure_output

# even if these are set in your env you need these when running the script
echo 'set the pgi env'
module unload cmake netcdf python netcdf-hdf5parallel
module swap PrgEnv-gnu PrgEnv-pgi
module load cmake python netcdf-hdf5parallel/4.2.0 subversion usg-default-modules/1.0
module load boost/1.49.0

# 0 is a successful build
setenv build_no 0

# NEEDED AFTER A FRESH CHECKOUT
echo 'bootstrap'
./bootstrap

if ($build_autoconf == 1 ) then
echo 'build with autoconf'
#SERIAL BUILD WITH AUTOCONF
echo 'make distclean'
make distclean 
echo 'configure serial autoconf build'
./configure-scripts/hopper-config-serial >& conf_serial.out
echo 'make serial'
make >& serial_build.out 

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy serial executable to test directory'
 cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_serial
else
 echo "autoconf parallel build failed, no executable"
 @ build_no = 1
endif

# PARALLEL BUILD WITH AUTOCONF PGI
echo 'make distclean'
make distclean
echo 'configure parallel autoconf build'
./configure-scripts/hopper-config-cesmtimers >& conf_auto_pgi.out
echo 'make parallel'
make >& auto_pgi_build.out 

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy autoconf parallel executable to test directory'
 cp -f $CODE_DIR/example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide
else
 echo "autoconf parallel build failed, no executable"
 @ build_no = 1
endif

else # build with autoconf option
 echo 'not building autoconf code option'
endif # build with autoconf option

echo $build_no
echo 'flag to build cmake option:' $build_cmake
if ($build_cmake == 1 ) then
# set up directories for building cmake executables
rm -rf xe6-pgi
rm -rf xe6-gnu
mkdir xe6-pgi
mkdir xe6-gnu
cp cmake-scripts/hopper-pgi-cmake-cesmtimers xe6-pgi
cp cmake-scripts/hopper-gnu-cmake-cesmtimers xe6-gnu

# PARALLEL BUILD WITH CMAKE PGI
cd $CODE_DIR/xe6-pgi
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure pgi cmake build'
./hopper-pgi-cmake-cesmtimers >& conf_cmake_pgi.out
echo 'make parallel pgi'
make -j 4 >& cmake_pgi_build.out

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy pgi parallel executable to test directory'
 cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_pgi
else
 echo "cmake pgi build failed, no executable"
 @ build_no = 1
endif

cd $CODE_DIR
make distclean

echo $build_no
# PARALLEL BUILD WITH CMAKE GNU
echo 'change to gnu env'
module unload cmake
module unload PrgEnv-cray PrgEnv-gnu PrgEnv-intel PrgEnv-pathscale PrgEnv-pgi
module unload hdf5 hdf5-parallel netcdf-hdf5parallel netcdf
module unload python cray-shmem cray-mpich2
module unload boost

#module --silent purge

module load PrgEnv-gnu/4.1.40
module swap gcc gcc/4.7.2

module load modules/3.2.6.6
module load cmake/2.8.7
module load python/2.7.1
module load cray-shmem/6.0.1
module load cray-mpich/6.0.1
module load torque/4.2.3.h5_notcpretry
module load netcdf-hdf5parallel/4.3.0
module load boost

module list

cd $CODE_DIR/xe6-gnu
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake
echo 'configure gnu cmake build'
./hopper-gnu-cmake-cesmtimers >& conf_gnu.out
echo 'make parallel gnu'
make -j 4 >& cmake_gnu_build.out

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy gnu parallel executable to test directory'
 cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_gnu
else
 echo "cmake gnu build failed, no executable"
 @ build_no = 1
endif

else # build with cmake option
 echo 'not building cmake code option'
endif # build with cmake option

echo $build_no
# execute tests on hopper
if ($build_no == 1 ) then
  echo "no job sumbitted, build/builds failed"
else
 echo 'submitting jobs to compute nodes'

  if ($REG_TEST == 0 ) then
    echo "no regression suite jobs submitted"
  else
  # simplest case, runs all builds and on a range of small processor counts
    echo 'submitting regression jobs to compute nodes'
    echo 'Go to carver.nersc.gov to complete Visualization and Verification (LIVV)'

  #diagnostic dome test case
    cd $TEST_DIR/reg_test/dome30/diagnostic
    qsub hopjob

  #evolving dome test case
    cd $TEST_DIR/reg_test/dome30/evolving
    qsub hopjob

  #confined shelf to periodic BC
    cd $TEST_DIR/reg_test/confined-shelf
    qsub hopjob

  #circular shelf to periodic BC
    cd $TEST_DIR/reg_test/circular-shelf
    qsub hopjob
  
  #ISMIP test case A 
    cd $TEST_DIR/reg_test/ismip-hom-a/80km
    qsub hopjob

  #ISMIP test case A 
    cd $TEST_DIR/reg_test/ismip-hom-a/20km
    qsub hopjob

  #ISMIP test case C - not operational until BC set
    #cd $TEST_DIR/reg_test/ismip-hom-c/80km
    #qsub hopjob

  #smaller GIS case to test realistic ice sheet configuration
    cd $TEST_DIR/reg_test/gis_10km
    qsub hopjob
  endif

# large but not challenging case, to test large processor counts, not yet configured for hopper
  if ($PERF_TEST == 0 ) then
    echo "no performance suite jobs submitted"
  else
    echo 'submitting performance jobs to compute nodes'
    echo 'Go to carver.nersc.gov to complete Visualization and Verification (LIVV)'

  #dome 30 test case
    cd $TEST_DIR/perf_test/dome30
    qsub hopjob

  #dome 60 test case
    cd $TEST_DIR/perf_test/dome60
    qsub hopjob

  #dome 120 test case
    cd $TEST_DIR/perf_test/dome120
    qsub hopjob

  #dome 240 test case
    cd $TEST_DIR/perf_test/dome240
    qsub hopjob

  #dome 500 test case
    cd $TEST_DIR/perf_test/dome500
    qsub hopjob

  #dome 1000 test case
    cd $TEST_DIR/perf_test/dome1000
    qsub hopjob
  endif

# high resolution GIS case to test realistic ice sheet configuration and longer time series, current setup gives
#convergence problems
#cd $TEST_DIR/gis_5km
#qsub hopjob

endif
