#!/bin/csh -f
/bin/rm glLoc.txt
touch glLoc.txt
foreach i ( plot.stream-shelf.0*0.* )
echo $i
./glfaces2d.Linux.64.g++.gfortran.DEBUG.ex $i inputs.BISICLES.stream-shelf > temp1.txt
./centerlineGL.awk < temp1.txt > temp2.txt
cat glLoc.txt temp2.txt > temp1.txt
/bin/cp -f temp1.txt glLoc.txt
end
