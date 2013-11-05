#!/bin/csh
# Master build script for mac laptops. Last updated 2/28/2013 by SFP.
# This is a hacked version of Kate's original script for use on Hopper.
# For now, only supports parallel build with Trilinos using gnu and cmake. 
# Only a subset of the small, standard tests are run, on both 1 and 4 procs.
 
# (1) execute from the main seacism directory
# (2) set the next two commands 

#add logic at the top to decide which versions to build 

# PARALLEL BUILD WITH CMAKE GNU

# !!! moved this to .bashrc !!!
# setenv TEST_DIR "/USERS/$USER/work/modeling/cism/seacism-oceans11/tests/higher-order"
setenv build_problem 0

echo
echo 'To use this script, type: csh mac-gnu-build-and-test.csh'
echo
echo 'For a quick test (dome only), type: csh mac-gnu-build-and-test.csh quick-test'

echo
echo 'See the LIVV documentation for instructions on setting up the test directory (TEST_DIR).'
echo
echo 'The following environment variables must be set: TEST_DIR, GLIMMER_TRILINOS_DIR'
echo 'Examples (place in .cshrc or .bashrc):'
echo 'csh, tcsh:  setenv GLIMMER_TRILINOS_DIR "/Users/$USER/Trilinos/gcc-build/install"'
echo 'bash:       export GLIMMER_TRILINOS_DIR="/Users/$USER/Trilinos/gcc-build/install"'
echo
echo 'Setting TEST_DIR to the standard location (bash):'
echo 'export TEST_DIR="/Users/$USER/work/modeling/cism/seacism-oceans11/tests/higher-order"'

# PARALLEL BUILD WITH CMAKE GNU

echo
echo "Configuring and building in directory: " $PWD
echo 

echo 'Configuring gnu cmake build...'
source ./mac-gnu-cmake >& conf_gnu.out
echo 'Making parallel gnu...'
make -j 4 >& cmake_gnu_build.out

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'copy gnu parallel executable to test directory'
 cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_gnu
else
 echo "cmake gnu build failed, no executable"
 @ build_problem = 1
endif

if ($build_problem == 1 ) then
  echo "No job sumbitted -- cmake build failed."
else  # execute tests:
 echo 'Submitting test jobs to compute nodes...'

 setenv run_all_tests 1
 if ($1 == "quick-test") setenv run_all_tests 0 

 #diagnostic dome test case
 cd $TEST_DIR/reg_test/dome30/diagnostic
 sh macjob

 if ($run_all_tests == 1) then

  #evolving dome test case
  cd $TEST_DIR/reg_test/dome30/evolving
  sh macjob

  # confined shelf to periodic BC
  cd $TEST_DIR/reg_test/confined-shelf
  sh macjob

  # circular shelf to periodic BC
  cd $TEST_DIR/reg_test/circular-shelf
  sh macjob

  # ISMIP test case A, 80 km 
  cd $TEST_DIR/reg_test/ismip-hom-a/80km
  sh macjob

  # ISMIP test case A, 20 km 
  cd $TEST_DIR/reg_test/ismip-hom-a/20km
  sh macjob

  ## ISMIP test case C - not operational 
  #cd $TEST_DIR/reg_test/ismip-hom-c/80km
  #sh macjob
 endif

 echo
 echo "Test runs completed."
 echo 

 cd $TEST_DIR/livv
 bash mac_VV_new.bash from-script $1

 echo
 echo "If there were errors finding ncl, add the ncl installation directory to your PATH in ~/.bashrc."
 echo

endif
