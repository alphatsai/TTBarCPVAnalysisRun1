#/bin/tcsh
if ( $3 == "" ) then
	echo "Please with following rule."	
	echo "	./reSubmitJob.csh [foderName] [samepleName] [numJob] <q=1nh>"	
	exit
endif
if ( !( -e $1 ) ) then
	echo "Here is no $1"
	exit	
endif
if ( !( -e $1/$2 ) ) then
	echo "Here is no $1/$2"
	exit	
endif
set q="1nh"
if ( $4 != "" ) then
	set q=$4 	
endif
cd $1
	set nowPath=`pwd`
	mv $2/output/job_$3.log $2
	rm -f $2/output/bprimeTobH_BkgABCDMethod_$3.root
	echo resubmit job_$3.sh...
	bsub -q $q -o $nowPath/$2/output/job_$3.log source $nowPath/$2/input/job_$3.sh
cd -
