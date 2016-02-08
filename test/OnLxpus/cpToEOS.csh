#!/bin/tcsh
if ( $1 == "" ) then
    echo "ERROR: Please input work folder name."
    echo "Ex: ./cpToEOS.csh [name] [deleting]"
    echo "    [deleting] deleteRoot, deleteLog, deleteAll"
    exit
endif
if ( ! ( -e $1 ) ) then
	echo "ERROR: Here is no work folder name $1 "
	exit
endif

cmsenv
set eosDir="eos/cms/store/user/jtsai/TTBarCPV/results/Jan2016"
    eos mkdir -p $eosDir/$1
    set roots=`\ls -l "$1" | grep root | awk '{print $9}'`
    set nroots=`\ls -l "$1" | grep root | awk '{print $9}' | wc -l `
    echo "   Copy $nroots roots to $eosDir/$1..."
    foreach root($roots)
        xrdcp $1/$root xroot://eoscms.cern.ch//$eosDir/$1/$root
        if ( $2 == 'deleteRoot' || $2 == 'deleteAll' ) then
            rm -f $1/$root
        endif
    end

set subDirs=`\ls -l $1 | grep drwxr-xr-x | awk '{print $9}'` 
echo $subDirs
foreach subDir( $subDirs )
        if ( $2 == 'deleteRoot' || $2 == 'deleteAll' ) then
            echo "[INFO] Deleting $1/$subDir"
            rm -rf $1/$subDir
        endif
#
#    set myDir="$1/$subDir"
#    set toEOS="$eosDir/$1/$subDir"
#
#    echo ">> In $myDir"
#    eos mkdir -p $toEOS
#
#    set roots=`\ls -l "$myDir/output" | grep root | awk '{print $9}'`
#    set nroots=`\ls -l "$myDir/output" | grep root | awk '{print $9}' | wc -l `
#    echo "   Copy $nroots roots to $toEOS..."
#    foreach root($roots)
#        xrdcp $myDir/output/$root xroot://eoscms.cern.ch//$toEOS/$root
#        if ( $2 == 'deleteRoot' || $2 == 'deleteAll' ) then
#            rm -f $myDir/output/$root
#        endif
#    end
#    #set logs=`\ls -l "$myDir/output" | grep  'job_' | grep '.log' | awk '{print $9}'`
#    #set nlogs=`\ls -l "$myDir/output" | grep 'job_' | grep '.log' | awk '{print $9}' | wc -l `
#    #echo "   Copy $nlogs logs to $toEOS..."
#    #foreach log($logs)
#    #    #xrdcp $myDir/output/$log xroot://eoscms.cern.ch//$toEOS/$log
#    #    if ( $2 == 'deleteLog' || $2 == 'deleteAll' ) then
#    #        rm -f $myDir/output/$log
#    #    endif
#    #end
end
