#!/bin/csh

# 1/29/2014 DMR -- This script has been modified to run in the builds/titan-pgi directory, using
# titan-pgi-cmake.

# Master build script for mac laptops. Last updated 2/28/2013 by SFP.
# This is a hacked version of Kate's original script for use on Hopper.
# For now, only supports parallel build with Trilinos using gnu and cmake. 
# Only a subset of the small, standard tests are run, on both 1 and 4 procs.
 
# (1) execute from the builds/mac-gnu subdirectory of CISM

#add logic at the top to decide which versions to build 

# PARALLEL BUILD WITH CMAKE 

# !!! moved this to .bashrc !!!
# setenv TEST_DIR "/USERS/$USER/work/modeling/cism/seacism-oceans11/tests/higher-order"

setenv TEST_DIR /lustre/atlas/scratch/$USER/cli062/higher-order

if (! -d $TEST_DIR) mkdir -p $TEST_DIR

setenv TEST_SUITE_DEFAULT_DIR /ccs/proj/cli062/test_suite

setenv build_problem 0

set COMPILER_NAME = pgi
set PLATFORM_NAME = titan

# set PLATFORM_NAME = $1
# set COMPILER_NAME = $2

set CMAKE_SCRIPT = $PLATFORM_NAME'-'$COMPILER_NAME'-cmake'
set CMAKE_CONF_OUT = 'conf_'$COMPILER_NAME'.out'
set CMAKE_BUILD_OUT = 'cmake_'$COMPILER_NAME'_build.out'
#set CISM_RUN_SCRIPT = $PLATFORM_NAME'job' 
set CISM_RUN_SCRIPT = 'ijob' 
#set CISM_VV_SCRIPT = $PLATFORM_NAME'_VV.bash'
set CISM_VV_SCRIPT = 'rhea_VV.bash'

echo
echo 'To use this script, type: csh '$PLATFORM_NAME'-'$COMPILER_NAME'-build-and-test.csh'
echo
echo 'For a quick test (dome only), type: csh '$PLATFORM_NAME'-'$COMPILER_NAME'-build-and-test.csh quick-test'

echo
echo 'See the LIVV documentation for instructions on setting up the test directory (TEST_DIR).'
echo
#echo 'The following environment variables must be set: TEST_DIR, GLIMMER_TRILINOS_DIR'
#echo 'Examples (place in .cshrc or .bashrc):'
#echo 'csh, tcsh:  setenv GLIMMER_TRILINOS_DIR "/Users/$USER/Trilinos/gcc-build/install"'
#echo 'bash:       export GLIMMER_TRILINOS_DIR="/Users/$USER/Trilinos/gcc-build/install"'
echo
echo 'Setting TEST_DIR to the location: '
echo 'TEST_DIR =' $TEST_DIR

# PARALLEL BUILD WITH CMAKE

echo
echo "Configuring and building in directory: " $PWD
echo 

echo 'Configuring '$COMPILER_NAME' cmake build...'
source ./$CMAKE_SCRIPT >& $CMAKE_CONF_OUT
echo 'Making parallel '$COMPILER_NAME'...'
make -j 8 >& $CMAKE_BUILD_OUT

if ( -e example-drivers/simple_glide/src/simple_glide ) then
 echo 'Copying '$COMPILER_NAME' parallel executable to test directory'
 cp -f example-drivers/simple_glide/src/simple_glide $TEST_DIR/simple_glide_$COMPILER_NAME
else
 echo "cmake '$COMPILER_NAME' build failed, no executable"
 @ build_problem = 1
endif

if ($build_problem == 1 ) then
  echo "No job submitted -- cmake build failed."
else  # execute tests:
 

 # Make copy of test suite in $TEST_DIR:

 echo "Copying default reg_test and LIVV to $TEST_DIR"
 pushd . > /dev/null
 cp $TEST_SUITE_DEFAULT_DIR/reg_test_default.tar $TEST_DIR/reg_test_default.tar
 cd $TEST_DIR
 tar xf reg_test_default.tar
 rm reg_test_default.tar
 popd > /dev/null

 cp -rf ../../tests/higher-order/livv $TEST_DIR


 echo 'Submitting test jobs to compute nodes.'

 setenv run_all_tests 1
 if ($1 == "quick-test") setenv run_all_tests 0 

 #diagnostic dome test case
 cd $TEST_DIR/reg_test/dome30/diagnostic
 qsub $CISM_RUN_SCRIPT


 if ($run_all_tests == 1) then

  #evolving dome test case
  cd $TEST_DIR/reg_test/dome30/evolving
  qsub $CISM_RUN_SCRIPT

  # confined shelf to periodic BC
  cd $TEST_DIR/reg_test/confined-shelf
  qsub $CISM_RUN_SCRIPT

  # circular shelf to periodic BC
  cd $TEST_DIR/reg_test/circular-shelf
  qsub $CISM_RUN_SCRIPT

  # ISMIP test case A, 80 km 
  cd $TEST_DIR/reg_test/ismip-hom-a/80km
  qsub $CISM_RUN_SCRIPT

  # ISMIP test case A, 20 km 
  cd $TEST_DIR/reg_test/ismip-hom-a/20km
  qsub $CISM_RUN_SCRIPT

  ## ISMIP test case C - not operational 
  #cd $TEST_DIR/reg_test/ismip-hom-c/80km
  #qsub $CISM_RUN_SCRIPT
 endif

 echo
 echo "Test Suite jobs started -- using qstat to monitor."
 echo 

 set still_running = 1
 set counter = 0
 set timeout_error = 0

 set run_list = "dome_30_test dome_30_evolve conf_shelf circ_shelf ishoma_80 ishoma_20"
 while ($still_running)
  set ls_out = `qstat | grep $USER`

  set found = 0 
  foreach cur ($run_list)
   foreach elem ($ls_out)
    if ("$cur" == "$elem") then
     if (($counter % 5) == 0) echo "Still running: $cur" 
     set found = 1
    endif
    # if ($found == 1) break 
   end 
  end
  if ($found == 0) then 
   echo "All jobs completed."
   set still_running = 0
  else 
   sleep 60
  endif
  @ counter = $counter + 1
  if ($counter == 120) then
   set still_running = 0
   set timeout_error = 1
   echo "Timeout error -- jobs are taking too long. Exiting script."
  endif
  if (($counter % 5) == 0) echo "Minutes: $counter"
 end

 if ($timeout_error == 0) then
  echo "Total minutes: $counter"
  echo

  echo "Call disabled to: $CISM_VV_SCRIPT, which is located in:" 
  echo "$TEST_DIR/livv"
  echo
  echo "Perform this step on rhea after the Test Suite jobs have completed."
  # cd $TEST_DIR/livv
  # bash $CISM_VV_SCRIPT from-script $1
 endif

 echo
 # echo "If there were errors finding ncl, add the ncl installation directory to your PATH in ~/.bashrc."
 echo

endif
