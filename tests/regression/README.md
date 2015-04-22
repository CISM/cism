========================
Regression testing suite
========================

The regression testing suite primarily intended to allow users and developers to
quickly generate a set of regression tests for use with the Land Ice Validation
and Verification (LIVV) toolkit. 

Typically, this test suit will be used to either test a new installation of this
Community Ice Sheet Model (CISM) by users, or to test new features by
developers. In either case, the workflow is essentially the same. 

============================
A few notes on building CISM
============================

In order for this test suite to work, you need a copy of CISM and all of CISM's
dependencies installed on your system. Throughout this document, we will assume
the CISM version you want to test is located at `$CISM`.  This test suite then
lives at `$CISM/tests/regression/`. 

It is recommended that you follow the 
[installation instructions](http://oceans11.lanl.gov/cism/documentation.html) 
as they are written in the manual at least once before using this suite. The
CISM build scripts are located in the directory
`$CISM/builds/PLATFORM-COMPILER[-DYCORE]` where `PLATFORM` is the computing
platform you are building on, `COMPILER` is the compiler you are using to build
CISM, and `DYCORE` is an optional specifier for any additional dycores to build
CISM with. Example:

```
$CISM/builds/mac-gnu-biscicles/
```

would be the location of the build scripts used to build CISM on a Mac computer
using the `gnu` compiler and the biscicles dycore. Within that directory there
will likely be a number of build scripts that follow the forms

```
PLATFORM-COMPILER[-DYCORE][-DYCORE-SPEC]-cmake
PLATFORM-COMPILER[-DYCORE][-DYCORE-SPEC]-cmake.bash
PLATFORM-COMPILER[-DYCORE][-DYCORE-SPEC]-serial-cmake
```

where DYCORE-SPEC may be a dash separated list of optional specifiers, and the
`-serial` specifies a serial build. 

NOTE: Currently, additional dycores are NOT supported by this test suit, but is
planned for future development. Serial builds are NOT supported by this test
suite, and there is no plan to support them. 

This test suite only makes use of the `PLATFORM-COMPILER-cmake.bash` scripts. If
there isn't a `*-cmake.bash` script for your `PLATFORM-COMPILER` combination,
then your combination is not yet supported. If you would like your combination
supported, contact the development team, or open an issue on the CISM github
page. 

NOTE:
-----
The `*-cmake.bash` build scripts should only differ slightly from the `*-cmake`
build scripts (that are referenced in the user manual) and add support for out
of source builds. If you use bash as your terminal, they should essentially be
drop in replacements. 

If you needed to change any of the variables in the `*-cmake` build script when
doing an initial build following the manual, you will also have to mirror those
changes in the equivilent `*-cmake.bash` script!

==================================
Additional test suite dependencies
==================================

This test suite requires you to have python 2.7 or greater, numpy, and
netcdf4-python installed (or loaded). Installing matplotlib and scipy is also
recommended but not required. 


======================
Running the test suite
======================

The test suite is accessed by running the `build_and_test.py` script in the
`$CISM/tests/regression/` directory.

The default run
---------------

By default, if you run `build_and_test.py` without any options, it will assume
you are on linux machine, want to use the `gnu` compiler, `$CISM` is `../..`
from your current working directory, and will then look for a build script named
`linux-gnu-cmake.bash` in `$CISM/builds/linux-gnu/`. 

NOTE: It is typically best to run the test suite from its own directory, and
hereafter `./` will refere to `$CISM/tests/regression/`.

The test suite will then create the directory `./build` and build CISM into it.
The CISM executable will then be located at `./build/cism_driver/cism_driver`. 

Next, the test suite will create the directory
`./reg_test/linux-gnu/` and will run the default tests and copy
their data to a `TEST` subdirectory(ies). 

reg test will have the structure:

```
reg_test
└── linux-gnu
    └── higher-order
    ├── dome
    │   ├── dome.RESO.pPRC.*
    │   ├── dome.RESO.pPRC.*
    ├── ismip-hom
    │   ├── ismip-hom-a.RESO.pPRC.*
    │   ├── ismip-hom-c.RESO.pPRC.*
    │   ├── ismip-hom-f.0100.pPRC.*
    ├── shelf
    │   ├── shelf-circular.RESO.pPRC.*
    │   ├── shelf-confined.RESO.pPRC.*
    └── stream
        └── stream.RESO.pPRC.*
```

where `RESO` is a four-digit number indicating the model resolution (units are
test specific), and `pPRC` is an optional three-digit number, prefixed by a `p`,
indicating the number of processors used to run the model.

NOTE: `RESO` will always equal `0100` for ISMIP-HOM test F as it is always run
at a 100 km resolution.

Detailed options
----------------

This script takes a few options:

```

  -h, --help            show the help message and exit.

  -p PLATFORM, --platform PLATFORM
                        Your computer platform. 
                        (default: linux)
  
  -c COMPILER, --compiler COMPILER
                        Your compiler. 
                        (default: gnu)
  
  -i CISM_DIR, --cism-dir CISM_DIR
                        Location of the CISM source code. 
                        (default: ../..)
  
  -b BUILD_DIR, --build-dir BUILD_DIR
                        Location to build CISM. 
                        (default: ./build)

  -j J                  Number of processors to use when making CISM. 
                        (default: 8)
  
  -s, --skip-build      Skip build, and look for the cism driver:
                            BUILD_DIR/cism_driver/cism_driver
                        (default: False)

  -o OUT_DIR, --out-dir OUT_DIR
                        Location of the directoy to output the test data.
                        (default: reg_test)

  -f, --force           Supress any warning about possibly overwriting test data.
                        (default: False)
  
  --timing              Run the timing test. This is needed for creating a new
                        benchmark dataset. 
                        (default: False)
  
  --performance         Run the performance tests. 
                        (default: False -- unless HPC system is detected, then
                                           will be forced to true always)

  --sleep SLEEP         Number of seconds to sleep between checks on running
                        tests. Only used on personal machines, when HPC systems 
                        are detected,it does not run the tests and just creates 
                        the batch jobs.
                        (default: 30)
```

You can see these options anytime by running `./build_and_test.py -h`. 

Multiple version runs
---------------------

If you would like to run the regression tests for a variety of different CISM
commits or versions, you can use the `-o/--out-dir` option to specify a different
`reg_test` directory to output the test data to. For example:

```
./build_and_test.py --platform mac \
                    --cism-dir CISM \
                    --build-dir build_VERSION \
                    --test_dir reg_test_VERSION
```

Will build CISM on a Mac, using the CISM code found in `CISM`,
into a directory named `build_VERSION` and output the data to
`reg_test_VERSION`. 


Authors
=======
Joseph H. Kennedy, ORNL 
    github : jhkennedy
