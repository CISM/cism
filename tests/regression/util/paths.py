import os
import sys
import errno
import fnmatch
import argparse
 
# small helper function so argparse will understand unsigned integers
def unsigned_int(x):
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type!"+ 
                "Should be an integer greater than or equal to zero.")
    return x

#NOTE: this parser is used to find the directory names needed for the reg_test
#      tree structure. The DOF and PROCESSOR directory levels are determined by
#      the individual test options defined in the util.dicts module
run_parser = argparse.ArgumentParser(prog='RUN',add_help=False)
run_parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
run_parser.add_argument('--scale', type=unsigned_int, default=0, 
        help="Change the degrees of freedom.")
run_parser.add_argument('--sizes', type=unsigned_int, metavar='KM',
        help="Change the domain size.")


def recursive_glob(tree, pattern):
    matches = []
    for base, dirs, files in os.walk(tree):
        goodfiles = fnmatch.filter(files, pattern)
        matches.extend(os.path.join(base, f) for f in goodfiles)
    return matches


def file_modifier_list(args):
    mod_list = []
    if args.tmod:
        mod_list.append(args.tmod)
    return mod_list


def make_absolute(args):
    args.out_dir = os.path.abspath(args.out_dir) 
    args.cism_dir = os.path.abspath(args.cism_dir) 
    args.build_dir = os.path.abspath(args.build_dir) 
    return args


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def case_run_directory(case_dir, run_args):
    case_run_dir = os.path.join(case_dir, 's'+str(run_args.scale), 'p'+str(run_args.parallel))
    if run_args.sizes is not None:
        case_run_dir += os.sep+'z'+str(run_args.sizes)
    return case_run_dir

def mkdir_test(args, test_dict):
    """
    Sets up a regression testing data directory that looks like:
        reg_test
        |-- PLATFORM-COMPILER
            |-- ICE_MODEL
                |-- TEST
                    |-- CASE
                        |-- DOF (degrees of freedom)
                            |-- PROCESSORS
                                |-- [OPTIONAL TEST SPECIFIC DIRS]
                                    |-- files.ext
    """

    #TODO: ice_model will be a indicator for which particular ice-sheet model was used.
    #      right now, this is a rather moot point as BATS only works for CISM-GLISADE
    ice_model = "CISM_glissade"
    data_dir = os.path.join(args.out_dir, args.platform, ice_model)

    # make the output directories softly
    mkdir_p(args.out_dir)
    mkdir_p(data_dir)
   
    # clean up run files, if the exist
    # NOTE: this is needed because runnit.hpc uses recursive_globs over the created
    #       run files in order to determine what run files are new. With force, all 
    #       run files will appear to be old. 
    if args.force:
        all_run_files = recursive_glob(data_dir,"*.run")
        for rf in all_run_files:
            os.remove(rf)


    for case in test_dict:
        case_split = str.split(case," ")
        case_dir = os.path.normpath(os.path.join(data_dir, str.split(case_split[0],"/")[-1], case_split[-1]))
        run_script, mod_dict = test_dict[case]
        
        run_args, ignore_args = run_parser.parse_known_args(str.split(run_script," ")+['--scale', '0', '-n', '1'])
        case_run_dir = case_run_directory(case_dir, run_args)

        mkdir_p(case_run_dir)
        if not args.force and os.listdir(case_run_dir):
            print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n") 
            print(  "WARNING: Test data already exists in:")
            print("\n"+case_run_dir+"\n")
            print(  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
            print(  "Some data may be overwritten. Either specify a different test directory, or ")
            print(  "re-run with the -f or --force option to ignore existing data.")
            print("\nExiting...")
            sys.exit(1)
        
        if args.performance and mod_dict:
            for mod in mod_dict:
                run_args, ignore_args = run_parser.parse_known_args(str.split(run_script," ")+mod_dict[mod].split())
                case_run_dir = case_run_directory(case_dir, run_args)

                mkdir_p(case_run_dir)
                if not args.force and os.listdir(case_run_dir):
                    print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n") 
                    print(  "WARNING: Test data already exists in:")
                    print("\n"+case_run_dir+"\n")
                    print(  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
                    print(  "Some data may be overwritten. Either specify a different test directory, or ")
                    print(  "re-run with the -f or --force option to ignore existing data.")
                    print("\nExiting...")
                    sys.exit(1)


    return data_dir


def cmake(args): 
    args.cmake_dir = []
    args.cmake_file = args.platform+'-cmake.sh'
    spec = args.platform.split("-")
    tried = []
    for i in reversed(range(len(spec)+1)):
        c_dir = os.path.join(args.cism_dir, 'builds', "-".join(spec[0:i]))
        if os.path.isdir(c_dir) and os.path.isfile(os.path.join(c_dir, args.cmake_file)):
            args.cmake_dir = c_dir
            break
        else:
            tried.append(c_dir)
            continue
	    
    if not args.cmake_dir:
        print("ERROR: cannot find your cmake file: ")
        print("    "+args.cmake_file)
        print("Tried looking in: ")
        for i in tried:
            print("    "+i)
        print("\nSee the cmake builds directory for supported platform-compiler combinations.")
        sys.exit(1)
    
    return args 

