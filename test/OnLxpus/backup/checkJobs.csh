#/bin/tcsh
echo "########################################################"
echo "###                                                  ###"
echo "### ./checkJob.csh [name] <rootname> <reSubmit> <q>  ###"
echo "###                                                  ###"
echo "########################################################"
if ( $1 == "" ) then
	echo "ERROR: Please input work folder name."
	echo "Ex: ./checkJobs.csh [name]"
	exit
endif
if ( ! ( -e $1 ) ) then
	echo "ERROR: Here is no work folder name $1 "
	exit
endif

cmsenv
#set start=`/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls eos/cms/store/user/jtsai/bpTobH/backgroundEstimationSkim | grep $1`
#if ( $start == "" ) then
#	echo "Nothing output..."
#	exit
#endif
if( $4 == "" ) then
	set sq="1nh" 
else
	set sq=$4
endif

if( $2 == "" ) then
    set rootname="SemiLeptanicAnalysis"
else
    set rootname=$2
endif

echo ">> Checking root file name $rootname..."
cd $1
	set nowPath=`pwd`
	rm -f tmp_.log tmp_check_.log
	set sampleName=`cat datasetList.txt | grep -v '#' | awk '{print $1}' | sed 's/^\///g' | sed 's/\//__/g'`
	set total=`echo $sampleName | wc -w`
	set doneS=0
	foreach sample($sampleName)
		touch tmp_.log tmp_check_.log
		set i=0
		set notDone=0
		echo "============================================================================================="
		echo "$sample"
		set jobNum=`ls -l $sample/input | grep '.sh' | wc -l`
		set lognum=`ls -l $sample/output | grep '.log' | wc -l`
		set killedJobs=`grep Killed $sample/output/*.log | grep -v 'cpu usage'| sed 's/.*job_\(.*\)\.sh.*/\1/g'`
		set kRoot=`grep 'Fatal Root Error' $sample/output/*.log | sed 's/.*job_\(.*\)\.log.*/\1/g'`
		set ksegJobs=`grep 'Segmentation' $sample/output/*.log | sed 's/.*job_\(.*\)\.sh.*/\1/g'`
		set kbadallocJobs_=`grep -r 'bad_alloc' $sample/output | sed 's/.*job_\(.*\)\.log.*/\1/g'`
		set kbadallocNum_=`echo $kbadallocJobs_ | wc -w` 
		set check=0
		if ( $kbadallocNum_ > 0  ) then
		while ( $check < $jobNum )
			set bad_=0
			foreach bad($kbadallocJobs_)
				if ( $bad == $check) then
					set bad_=1
				endif
			end
			if ( $bad_ == 1 ) then
				echo $check >> tmp_check_.log
			endif
			@ check++
		end
		endif
		set kbadallocJobs=`cat tmp_check_.log` 
		set kCPUJobs=`grep Killed $sample/output/*.log | grep 'cpu usage' | sed 's/.*job_\(.*\)\.log.*/\1/g'`
		set kCPUJobs2=`grep 'CPU time limit exceeded' $sample/output/*.log | grep 'sh:Exited' | sed 's/.*job_\(.*\)\.sh.*/\1/g'`
		set kCPUJobs3=`grep 'CPU time limit exceeded' $sample/output/*.log | grep 'sh: line' | sed 's/.*job_\(.*\)\.sh.*/\1/g'`
		set abJobs=`grep Aborted $sample/output/*.log | sed 's/.*job_\(.*\)\.sh.*/\1/g'`
		set kEOSRoot=`grep 'Error while doing the asyn writing' $sample/output/*.log | sed 's/.*job_\(.*\)\.log.*/\1/g'`
                #set failOpens=`grep 'Failed to open the file' $sample/output/*.log | sed 's/.*job_\(.*\)\.log.*/\1/g'`
                set failOpens=`grep -r 'Failed to open the file' $sample/output | sed 's/.*job_\(.*\)\.log.*/\1/g'`
                set eosCopys=`grep -r 'error: target file was not created!' $sample/output | sed 's/.*job_\(.*\)\.log.*/\1/g'`

                set failOpensNum=`echo $failOpens | wc -w`
                set eosCopysNum=`echo $eosCopys | wc -w`
		set kRootNum=`echo $kRoot | wc -w`
		set kEOSRootNum=`echo $kEOSRoot | wc -w`
		set ksegNum=`echo $ksegJobs | wc -w `
		set kbadallocNum=`echo $kbadallocJobs| wc -w `
		set kCPUNum=`echo $kCPUJobs | wc -w `
		set kCPUNum2=`echo $kCPUJobs2 | wc -w `
		set kCPUNum3=`echo $kCPUJobs3 | wc -w `
		set killedNum=`echo $killedJobs | wc -w `
		set abNum=`echo $abJobs | wc -w `
		set doneJobs=`ls -l $sample/output | grep root | awk '{print $9}' | sed "s/""$rootname""_\(.*\)\.root/\1/g"`
		set doneNum=`echo $doneJobs | wc -w`	
		set realdoneNum=`echo $doneNum'-'$killedNum'-'$abNum'-'$kCPUNum2'-'$kCPUNum3'-'$ksegNum'-'$kbadallocNum'-'$kRootNum'-'$failOpensNum'-'$eosCopysNum | bc`
		echo "Num Log Files $lognum/$jobNum"	
		echo "Status(root): $doneNum/$jobNum"
		echo "Status(real): $realdoneNum/$jobNum"
		if ( $kRootNum != 0 ) then
			echo "Input root broken: $kRoot"
		endif
		if ( $doneNum == 0 && $lognum == 0 ) then
			echo "Nothing output..."	
		else if ( $realdoneNum == $jobNum ) then
			@ doneS++
			echo "Done!"	
		else
			while ( $i < $jobNum )
				set done=0
				foreach job($doneJobs)	
					if ( $i == $job ) then
						set done=1
						#echo "error $i $done"	
					endif	
					#echo $i "'"$job"'" $done
				end
				if ( $done == 0 ) then
					#echo $i
					echo $i >> tmp_.log
					@ notDone++ 
				endif	
				@ i++
			end
		endif
		if ( $kCPUNum != 0 ) then
			echo "CPU Use Jobs: "$kCPUJobs 
		endif
		if ( $kCPUNum2 != 0 || $kCPUNum3 != 0 ) then
			echo "CPU Time Jobs: "$kCPUJobs2 $kCPUJobs3 
		endif
		if ( $killedNum != 0 ) then
			echo "Killed Jobs: "$killedJobs 
		endif
		if ( $abNum != 0 ) then
			echo "Aborted Jobs: "$abJobs 
		endif
		if ( $ksegNum != 0 ) then
			echo "Segmetation Jobs: "$ksegJobs 
		endif
		if ( $kbadallocNum != 0 ) then
			echo "bad_alloc: "$kbadallocJobs 
		endif
                if ( $kEOSRootNum != 0 ) then
                        echo "To eos fail: "$kEOSRoot
                endif
                if ( $failOpensNum != 0 ) then
                        echo "Failed open: "$failOpens
                endif
                if ( $eosCopysNum != 0 ) then
                        echo "Failed copy: "$eosCopys
                endif
		if ( $notDone != 0 && $rootname != 'Skim' ) then
			set notDonelist=`cat tmp_.log`	
			echo "No root Jobs: "$notDonelist 
		endif
		rm -f tmp_.log tmp_check_.log

		if ( $3 == 'reSubmit' && $notDone != 0 && $rootname != 'Skim' ) then
			foreach nn($notDonelist)
				mv $nowPath/$sample/output/job_$nn.log $nowPath/$sample
				echo resubmit job_$nn.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$nn.log source $nowPath/$sample/input/job_$nn.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $kRootNum != 0 ) then
			foreach kn($kRoot)
				mv $nowPath/$sample/output/job_$kn.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$kn.root
				echo resubmit job_$kn.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$kn.log source $nowPath/$sample/input/job_$kn.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $killedNum != 0 ) then
			foreach kn($killedJobs)
				mv $nowPath/$sample/output/job_$kn.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$kn.root
				echo resubmit job_$kn.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$kn.log source $nowPath/$sample/input/job_$kn.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $kCPUNum2 != 0 ) then
			foreach kcn($kCPUJobs2)
				mv $nowPath/$sample/output/job_$kcn.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$kcn.root
				echo resubmit job_$kcn.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$kcn.log source $nowPath/$sample/input/job_$kcn.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $kCPUNum3 != 0 ) then
			foreach kcn3($kCPUJobs3)
				mv $nowPath/$sample/output/job_$kcn3.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$kcn3.root
				echo resubmit job_$kcn3.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$kcn3.log source $nowPath/$sample/input/job_$kcn3.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $abNum != 0 ) then
			foreach an($abJobs)
				mv $nowPath/$sample/output/job_$an.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$an.root
				echo resubmit job_$an.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$an.log source $nowPath/$sample/input/job_$an.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $failOpensNum != 0 ) then
			foreach an($failOpens)
				mv $nowPath/$sample/output/job_$an.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$an.root
				echo resubmit job_$an.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$an.log source $nowPath/$sample/input/job_$an.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $eosCopysNum != 0 ) then
			foreach an($eosCopys)
				mv $nowPath/$sample/output/job_$an.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$an.root
				echo resubmit job_$an.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$an.log source $nowPath/$sample/input/job_$an.sh
			end	
		endif
		if ( $3 == 'reSubmit' && $kbadallocNum != 0 ) then
			foreach an($kbadallocJobs)
				mv $nowPath/$sample/output/job_$an.log $nowPath/$sample
				rm -f $sample/output/$rootname'_'$an.root
				echo resubmit job_$an.sh...
				bsub -q $sq -o $nowPath/$sample/output/job_$an.log source $nowPath/$sample/input/job_$an.sh
			end	
		endif
	end
	echo "============================================================================================="
        echo ">> Checking root file name $rootname.root"
	echo ">> Summerize: $doneS/$total"
	rm -f tmp_.log
cd -

