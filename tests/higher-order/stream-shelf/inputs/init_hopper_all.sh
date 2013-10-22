#!/bin/ksh

#commonly changed parameters
NEND=800 # number of iterations we're doing
RUNHOME=/scratch2/scratchdirs/dmartin/Goldberg/bisicles-rerun/bisicles
#RUNHOME=.
#directory containing input nc files



INFILE_TEMPLATE=inputs.BISICLES.template
CONFIGFILE_TEMPLATE=stream-shelf.config.template
QSUBFILE_TEMPLATE=hopper.qsb.template


#function to get dt
getdt()
{
  case $STEPSPERYEAR in
    12) DT=0.08333333333333333;;
    24) DT=0.04166666666666667;;
    3)  DT=0.33333333333333333;;
    6)  DT=0.16666666666666667;;
    *) echo "Unanticipated STEPSPERYEAR val for dt";;
esac
}

for STEPSPERYEAR in 3 
do
echo "generating inputs for" $STEPSPERYEAR " steps per year"

POPRESULTSDIR=popResults
PERYEAR=_per_year
DIR=$RUNHOME/$STEPSPERYEAR$PERYEAR/
INFILE1_BASE=inputs.BISICLES.$STEPSPERYEAR$PERYEAR
CONFIGFILE_BASE=stream-shelf.$STEPSPERYEAR$PERYEAR
CONFIGSUFFIX=config
QSBSUFFIX=qsb

#INDIR=.
#INDIR="\/global\/project\/projectdirs\/iceocean\/bisicles-pop\/XylarToDan\/$STEPSPERYEAR$PERYEAR"
#INDIR="\/scratch\/scratchdirs\/xylar\/Goldberg\/coupling\/$STEPSPERYEAR$PERYEAR\/popResults"
INDIR="\/scratch2\/scratchdirs\/dmartin\/Goldberg\/bisicles-rerun\/xylar\/$STEPSPERYEAR$PERYEAR"
#directory containing output nc files
#OUTDIR=.
#OUTDIR="\/global\/project\/projectdirs\/iceocean\/bisicles-pop\/DanToXylar\/$STEPSPERYEAR$PERYEAR"
OUTDIR="\/global\/project\/projectdirs\/iceocean\/bisicles-pop\/DanToXylar\/$STEPSPERYEAR$PERYEAR"


echo $DIR
mkdir -p $DIR

 getdt   
 echo "dt = " $DT

 TSTART=0.0
 TEND=$DT
 
LASTSTEP=""
NSTEP=0
while [ $NSTEP -le $NEND ]
do 
#  echo $NSTEP
STEP=$(printf "%05i" $NSTEP)
#  echo $STEP

#    TSTART=$(($STEP * $DT))
#   echo $TSTART

    of1=$DIR$INFILE1_BASE.$STEP
    of2=$DIR$CONFIGFILE_BASE.$STEP.$CONFIGSUFFIX

#    echo $of1
#    echo $of2
#    echo $of3
#    echo $of4

    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/   $INFILE_TEMPLATE > $of1
    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/ -e s/@DT/$DT/ -e s/@TEND/$TEND/ -e s/@TSTART/$TSTART/ -e s/@INDIR/$INDIR/ -e s/@OUTDIR/$OUTDIR/  $CONFIGFILE_TEMPLATE > $of2

    TSTART=$TEND
    NSTEP=$((NSTEP + 1))
    TEND=$((TEND + DT))
    LASTSTEP=$STEP

#    echo $NSTEP
done

done 

#generate hopper qsub file

echo "generating qsub file for hopper"
QSUBFILE=$RUNHOME/BISICLES-all.$STEPSPERYEAR$PERYEAR.hopper

#clean up pre-existing files
rm -f $QSUBFILE

echo "#PBS -q regular" >> $QSUBFILE
echo "#PBS -l mppwidth=24" >> $QSUBFILE
echo "#PBS -l walltime=6:00:00" >> $QSUBFILE
echo "#PBS -N BISICLES-all-$STEPSPERYEAR$PERYEAR" >> $QSUBFILE
echo "#PBS -V " >> $QSUBFILE
echo "#PBS -m abe" >> $QSUBFILE
echo "  " >> $QSUBFILE
echo "cd \$PBS_O_WORKDIR/$STEPSPERYEAR$PERYEAR" >> $QSUBFILE

LASTSTEP=""
#NSTEP=2
NSTEP=98

while [ $NSTEP -le $NEND ]
do 
#  echo "step = "  $NSTEP
  STEP=$(printf "%05i" $NSTEP)

  echo "aprun -n 16 ../simple_bisicles stream-shelf.3_per_year.$STEP.config " >> $QSUBFILE

    NSTEP=$((NSTEP + 1))

done

exit 0


