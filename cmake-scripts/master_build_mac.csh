#!/bin/csh
# Master build script for mac laptops. Last updated 8/2013 
# For now, only supports parallel build with Trilinos using gnu and cmake. 
# Only a subset of the small, standard tests are run, on both 1 and 4 procs.

# (1) execute from the main seacism directory
# (2) set the next two commands 
#add logic at the top to decide which versions to build 

# change these to mathc your env (Steves by default)
setenv TEST_DIR "~/work/modeling/cism/seacism-livv/tests/higher-order"
setenv CODE_DIR "~/work/modeling/cism/seacism-livv"

cd $CODE_DIR
setenv build_no 0
setenv build_cmake 0

echo $build_no
echo 'flag to build cmake option:' $build_cmake

if ($build_cmake == 1 ) then    
# make directories for building cmake executables
rm -rf mac-gnu
mkdir mac-gnu
cp cmake-scripts/mac-cmake mac-gnu/

echo $build_no

# PARALLEL BUILD WITH CMAKE GNU

cd $CODE_DIR/mac-gnu
echo 'clean out the build dir'
rm -rf CMakeCache.txt CMakeFiles fortran_mod_files lib libglimmer-trilinos autogenerate.log cmake_install.cmake

echo 'configure gnu cmake build'
./mac-cmake >& conf_gnu.out

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
# execute tests

if ($build_no == 1 ) then
  echo "no job sumbitted, build/builds failed"
else

# simplest case, runs build and on a range of small processor counts 
 echo 'submitting jobs to compute nodes'

#diagnostic dome test case
cd $TEST_DIR/reg_test/dome30/diagnostic
sh macjob

#evolving dome test case
#cd $TEST_DIR/reg_test/dome30/evolving
#sh macjob

## ISMIP test case A - not operational until BC set
cd $TEST_DIR/reg_test/ismip-hom-a/80km
sh macjob

## ISMIP test case A - not operational until BC set
cd $TEST_DIR/reg_test/ismip-hom-a/20km
sh macjob

## ISMIP test case C - not operational until BC set
#cd $TEST_DIR/reg_test/ismip-hom-c/80km
#sh macjob

## confined shelf to periodic BC
cd $TEST_DIR/reg_test/confined-shelf
sh macjob

## circular shelf to periodic BC
cd $TEST_DIR/reg_test/circular-shelf
sh macjob

## GIS 10km test, 1 time step
#cd $TEST_DIR/reg_test/gis_10km
#sh macjob

endif
