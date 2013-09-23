#!/bin/ksh
INFILE_TEMPLATE=inputs.BISICLES.template
CONFIGFILE_TEMPLATE=stream-shelf.config.template
FIXNC_INFILE_TEMPLATE=inputs.fixNC.template
FIXNC_CONFIGFILE_TEMPLATE=fixNC.config.template


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
echo "generating inputs for" $STEPSPERYEAR

PERYEAR=_per_year
DIR=$STEPSPERYEAR$PERYEAR/
INFILE1_BASE=inputs.BISICLES.$STEPSPERYEAR$PERYEAR
CONFIGFILE_BASE=stream-shelf.$STEPSPERYEAR$PERYEAR
FIXNC_INFILE_BASE=inputs.fixNC.$STEPSPERYEAR$PERYEAR
FIXNC_CONFIGFILE_BASE=fixNC.$STEPSPERYEAR$PERYEAR
CONFIGSUFFIX=config
mkdir -p $DIR

 getdt   
 echo "dt = " $DT

 TSTART=0.0
 TEND=$DT

LASTSTEP=""
NSTEP=0
for STEP in 00000 00001 00002 00003 00004 00005 00006 00007 00008 
do
    
#    TSTART=$(($STEP * $DT))
#   echo $TSTART

    of1=$DIR$INFILE1_BASE.$STEP
    of2=$DIR$CONFIGFILE_BASE.$STEP.$CONFIGSUFFIX
    of3=$DIR$FIXNC_INFILE_BASE.$STEP
    of4=$DIR$FIXNC_CONFIGFILE_BASE.$STEP.$CONFIGSUFFIX

    echo $of1
    echo $of2
    echo $of3
    echo $of4

    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/   $INFILE_TEMPLATE > $of1
    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/ -e s/@DT/$DT/ -e s/@TEND/$TEND/ -e s/@TSTART/$TSTART/   $CONFIGFILE_TEMPLATE > $of2
    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/ -e s/@DT/$DT/ -e s/@TEND/$TEND/ -e s/@TSTART/$TSTART/   $FIXNC_INFILE_TEMPLATE > $of3
    sed  -e s/@STEPSPERYEAR/$STEPSPERYEAR/ -e s/@STEP/$STEP/ -e s/@LASTSTEP/$LASTSTEP/ -e s/@DT/$DT/ -e s/@TEND/$TEND/ -e s/@TSTART/$TSTART/   $FIXNC_CONFIGFILE_TEMPLATE > $of4

    TSTART=$TEND
    NSTEP=$((NSTEP + 1))
    TEND=$((TEND + DT))
    LASTSTEP=$STEP

    echo $NSTEP
done

done 

exit 1
    outfile1="run.petsc-gamg.$RES"
    innerConvergename1="solverConverge/resid.petsc-gamg.$RES"
    outerConvergename1="solverConverge/resid.petsc-gamg.$RES.outer"
    poutname1="pout.petsc.l0.$RES.0"
    runcommand1="mpirun -np $NPROC $EXECFILE1 $of1 > $outfile1"
    echo "echo \"doing $RES run\" " >> $RUNFILE1
    echo $runcommand1 >> $RUNFILE1
    echo "$SCRIPTDIR/innerJFNK.awk < $poutname1 > $TEMPFILE1 " >> $RUNFILE1
    echo "$SCRIPTDIR/parseMg $TEMPFILE1  $innerConvergename1" >> $RUNFILE1
    echo "$SCRIPTDIR/jfnk.awk < $poutname1 > $TEMPFILE1" >> $RUNFILE1
    echo "$SCRIPTDIR/parseJFNK  $TEMPFILE1  $outerConvergename1" >> $RUNFILE1

    outfile2="run.MG-JFNK.$RES"
    innerConvergename2="solverConverge/resid.MG-JFNK.$RES"
    outerConvergename2="solverConverge/resid.MG-JFNK.$RES.outer"
    runcommand2="mpirun -np $NPROC $EXECFILE2 $of2 > $outfile2"
    echo "echo \"doing $RES run\" " >> $RUNFILE2
    poutname2="pout.MG-JFNK.l0.$RES.0"
    echo $runcommand2 >> $RUNFILE2
    echo "$SCRIPTDIR/innerJFNK.awk < $poutname2 > $TEMPFILE2 " >> $RUNFILE2
    echo "$SCRIPTDIR/parseMg $TEMPFILE2  $innerConvergename2" >> $RUNFILE2
    echo "$SCRIPTDIR/jfnk.awk < $poutname2 > $TEMPFILE2" >> $RUNFILE2
    echo "$SCRIPTDIR/parseJFNK  $TEMPFILE2  $outerConvergename2" >> $RUNFILE2



CRSERES=$RES
done 

done

exit 0


