#!/bin/bash

MYWORKDIR="/afs/cern.ch/work/j/jtsai/caculatePileUp/Run2012/CMSSW_7_4_2/src/OnLxplus"
MYJSON="Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt" # Must put under MYWORKDIR
#XINGMINLUMI=0.3
XINGMINLUMI=0.1
BATCHDIR=${PWD}

#########################################
export SCRAM_ARCH=slc6_amd64_gcc481
cd $MYWORKDIR 
eval `scram runtime -sh`

#########################################
cp -v "$MYWORKDIR/$MYJSON" $BATCHDIR

cd $BATCHDIR
echo "Running CMSSW job"
lumiCalc2.py lumibylsXing --xingMinLum $XINGMINLUMI -b stable -i $MYJSON -o lumi.csv
estimatePileup_makeJSON.py --csvInput lumi.csv  pileup_JSON.txt
exitcode=$?

echo "Copying file lumi.csv to $MYWORKDIR/lumi_XINGMINLUMI$XINGMINLUMI.csv"
cp -v lumi.csv "$MYWORKDIR/lumi_XINGMINLUMI$XINGMINLUMI.csv"
echo "Copying file pileup_JSON.txt to $MYWORKDIR/pileup_JSON_XINGMINLUMI$XINGMINLUMI.txt"
cp -v pileup_JSON.txt "$MYWORKDIR/pileup_JSON_XINGMINLUMI$XINGMINLUMI.txt"

exit $exitcode

