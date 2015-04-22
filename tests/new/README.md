How to create a new CISM test
=============================

This directory contains a new test template called `test.py` which
should be used to create a new test in CISM. 

Steps to creating a new test
----------------------------

Say we want to create a higher-order test for Antarctica, which we will 
call `antarctica` The first thing we do is make the Antarctica test
directory, and then move into it:

```bash
$ mkdir $CISM/tests/higher-order/antarctica
$ cd $CISM/tests/higher-order/antarctica
```

Next, we copy the `runTest.py` python script template into your new test 
directory, as well as the netCDF tools module `netCDF.py` from the 
`$CISM/tests/new` directory. 

```bash
$ cp $CISM/tests/new/runTest.py runAntarctica.py
$ cp $CISM/test/new/netCDF.py ./
```

Next, we open up the `runAntarctica.py` script in our favorite editor 
and create our test. Most of the structure of the test has been created
from the template. To finish the test, we find the `#FIXME:` statements 
in the `runAntarctica.py`, read the associated comment text, and create
what we need for the test case. 




*Last updated: April 27, 2015 by Joseph H Kennedy at ORNL*
