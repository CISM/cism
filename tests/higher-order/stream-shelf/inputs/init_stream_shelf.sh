#!/bin/ksh

#commonly changed parameters
NEND=30 # number of iterations we're doing
#RUNHOME=/scratch/users/dmartin/newXylar
RUNHOME=.
#directory containing input nc files



INFILE_TEMPLATE=inputs.BISICLES.template
CONFIGFILE_TEMPLATE=stream-shelf.config.template
FIXNC_INFILE_TEMPLATE=inputs.fixNC.template
FIXNC_CONFIGFILE_TEMPLATE=fixNC.config.template
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

for STEPSPERYEAR in 12 24 6 3 
do
echo "generating inputs for" $STEPSPERYEAR " steps per year"

PERYEAR=_per_year
DIR=$RUNHOME/$STEPSPERYEAR$PERYEAR/
INFILE1_BASE=inputs.BISICLES.$STEPSPERYEAR$PERYEAR
CONFIGFILE_BASE=stream-shelf.$STEPSPERYEAR$PERYEAR
FIXNC_INFILE_BASE=inputs.fixNC.$STEPSPERYEAR$PERYEAR
FIXNC_CONFIGFILE_BASE=fixNC.$STEPSPERYEAR$PERYEAR
CONFIGSUFFIX=config
QSBSUFFIX=qsb

#INDIR=.
INDIR="\/global\/project\/projectdirs\/iceocean\/bisicles-pop\/XylarToDan\/$STEPSPERYEAR$PERYEAR"
#directory containing output nc files
#OUTDIR=.
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
STEP=$(printf "%05i" $NSTEP)
  echo $STEP
    
#    TSTART=$(($STEP * $DT))
#   echo $TSTART

    of1=$DIR$INFILE1_BASE.$STEP
    of2=$DIR$CONFIGFILE_BASE.$STEP.$CONFIGSUFFIX
    of3=$DIR$FIXNC_INFILE_BASE.$STEP
    of4=$DIR$FIXNC_CONFIGFILE_BASE.$STEP.$CONFIGSUFFIX

#    echo $of1
#    echo $of2
#    echo $of3
#    echo $of4

    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/   $INFILE_TEMPLATE > $of1
    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/ -e s/@DT/$DT/ -e s/@TEND/$TEND/ -e s/@TSTART/$TSTART/ -e s/@INDIR/$INDIR/ -e s/@OUTDIR/$OUTDIR/  $CONFIGFILE_TEMPLATE > $of2
    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/ -e s/@DT/$DT/ -e s/@TEND/$TEND/ -e s/@TSTART/$TSTART/  -e s/@INDIR/$INDIR/ -e s/@OUTDIR/$OUTDIR/  $FIXNC_INFILE_TEMPLATE > $of3
    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/ -e s/@DT/$DT/ -e s/@TEND/$TEND/ -e s/@TSTART/$TSTART/ -e s/@INDIR/$INDIR/ -e s/@OUTDIR/$OUTDIR/   $FIXNC_CONFIGFILE_TEMPLATE > $of4

    TSTART=$TEND
    NSTEP=$((NSTEP + 1))
    TEND=$((TEND + DT))
    LASTSTEP=$STEP

#    echo $NSTEP
done

done 

#generate hopper qsub files

echo "generating qsub files for hopper"
QSUBFILE_BASE=BISICLES.hopper
LASTSTEP=""
NSTEP=0
while [ $NSTEP -le $NEND ]
do 
#  echo $NSTEP
  if [ $NSTEP -le 9 ]
  then
    STEP="0000"$NSTEP
  elif [ $NSTEP -lt 99 ]
  then 
    STEP="000"$NSTEP
  elif [ $NSTEP -lt 999 ]
  then 
    STEP="00"$NSTEP
  elif [ $NSTEP -lt 9999 ]
  then 
    STEP="0"$NSTEP
  else
    STEP=$n
  fi

    qfile=$RUNHOME/$QSUBFILE_BASE.$STEP

    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@FIXSTEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/   $QSUBFILE_TEMPLATE > $qfile

    NSTEP=$((NSTEP + 1))

done

exit 0


