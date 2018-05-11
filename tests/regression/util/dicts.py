"""A set of dictionaries for the regression tests.

These dictionaries will determine what tests are run and what options to pass
to each test. See descriptions of each dictionary.
"""

# SPECIFIC TESTS DICTIONARIES
# ===========================
# These dictionaries describe the NON DEFAULT options to pass to the test run
# scripts found within the CISM tests directory. These modified test cases are
# only called if the performance option is specified.
# 
# Dictionaries for use in the MAIN TEST DICTIONARY (at bottom). They should
# consist of key-value pairs like:
#   key: anything descriptive -- used only internally
#   value: a string containing the options to pass to the script. 
# See the dictionaries below for an example.
# NOTE: these options should not include the -e/--executable, --hpc, 
#       -o/--output-dir, -m/--modifier, or -s/--setup-only options because the
#       build_and_test script will include those automatically.

# dome test
# --------------------------
# for tests with -n N for N < 16
dome_perf_small = { 
        's1': '--scale 1 -n 1',
        's2': '--scale 2 -n 1',
        'p1': '-n 2',
        'p2': '-n 4',
        'b1': '--scale 1 -n 4',
        'p3': '-n 8',
        'b2': '--scale 2 -n 16',
        }

# for tests with -n N for N > 16
dome_perf_large = {
        's3': '--scale 2 -n 256',
        's4': '--scale 3 -n 256',
        'b3': '--scale 3 -n 64',
        'b4': '--scale 4 -n 256',
        }

# shelf tests
# ----------------------------
# NOTE: empty dict because no performance testing for confined shelf. Leaving
#       here for possible future expansion. 
shelfConfined_perf_small = {}

# NOTE: empty dict because no performance testing for circular shelf. Leaving
#       here for possible future expansion. 
shelfCircular_perf_small = {}

# NOTE: empty dict because no performance testing for ISMIP-HOM. Leaving here
#       for possible future expansion. 
ismip_perf_small = {}

# NOTE: empty dict because no performance testing for stream. Leaving here for
#       possible future expansion. 
stream_perf_small = {}

# NOTE: this is a purposely empty dict to be used for repeated default tests
#       where only one performance dictionary is needed (e.g., ISMIP-HOM). 
keep_empty = {}

# MAIN TEST DICTIONARY
# ====================
# This is the main dictionary that describes what tests to run.
# Each dictionary item should consist of key-value pairs like:
#   key: 'path_to_test_from_$CISM/tests SIZE(optional) CASE'
#       example: 'dome'
#       NOTE: key can be a space separated list with the first entry the path to
#             the test directory, the last entry is the specific test case, and 
#             the rest of the list used to define uniqueness. This is useful for 
#             tests that have multiple run scripts like shelf or ISMIP-HOM.
#           example: 'ismip-hom 20 a'
#       
#   value: tuple of (run_script, perf_dict) where run_script is the test run
#          script that can be found within the directory specified by the key 
#          and necessary command line options, and perf_dict is a dictionary
#          containing the options to pass to therun_script when doing 
#          performance testing (as described above)
#        example: ('runISMIP_HOM.py -r a --size 20', dome_perf_small)
#
#NOTE: foe ISMIP-HOM, experiment f, there is only one size (100km) and instead 
#      experiments are run with different silp rations. CISM only runs the no-slip 
#      (ratio = 0) experiment, and as such, the uniquness part of the key and the 
#      size option in the run command reflects the slip ratio not the domain size 
#      (like in the other ISMIP-HOM tests).
test_dict = {
        'dome dome': ('runDome.py', dome_perf_small),
        'shelf shelf-confined': ('runShelfConfined.py', shelfConfined_perf_small),
        'shelf shelf-circular': ('runShelfCircular.py', shelfConfined_perf_small),
        'ismip-hom 20 ismip-hom-a': ('runISMIP_HOM.py -r a --size 20', ismip_perf_small),
        'ismip-hom 20 ismip-hom-c': ('runISMIP_HOM.py -r c --size 20', keep_empty),
        'ismip-hom 80 ismip-hom-a': ('runISMIP_HOM.py -r a --size 80', keep_empty),
        'ismip-hom 80 ismip-hom-c': ('runISMIP_HOM.py -r c --size 80', keep_empty),
        'ismip-hom 0 ismip-hom-f': ('runISMIP_HOM.py -r f --size 0', keep_empty),
        'stream stream': ('runStream.py', stream_perf_small),
        }

perf_dict = {
        'dome dome': ('runDome.py', dome_perf_large),
        }

# HPC PLATFORM DICTIONARIES
# =========================
# There should be a dictionary for each supported HPC platform which specifies 
# the default batch scheduler options to use. Additionally, there is an
# optional 'RUN_CMD' key whose value will replace the run command ('aprun' is
# the default).
hopper_dict = {
        'PBS_A': 'm1795',
        'PBS_q': 'regular',
        'PBS_N': 'reg_test_all',
        'PBS_RES': 'mppwidth',
        'RES_NUM': '24',
        'PBS_walltime': '01:00:00',
        }


titan_dict = {
        'PBS_A': 'cli106ice',
        'PBS_q': 'batch',
        'PBS_N': 'reg_test_all',
        'PBS_RES': 'nodes',
        'RES_NUM': '1',
        'PBS_walltime': '01:00:00',
        }

cheyenne_dict = {
        'PBS_A': 'P93300601',
        'PBS_q': 'regular',
        'PBS_N': 'reg_test_all',
        'PBS_RES': 'select',
        'RES_NUM': '2:ncpus=36:mpiprocs=36',
        'PBS_walltime': '01:00:00',
        'RUN_CMD': 'mpiexec_mpt',
        }


# MAIN HPC DICTIONARY
# ===================
# Collection of all the HPC platform dictionaries. 
hpc_dict = {
        'titan': titan_dict,
        'hopper': hopper_dict,
        'cheyenne': cheyenne_dict,
        }
