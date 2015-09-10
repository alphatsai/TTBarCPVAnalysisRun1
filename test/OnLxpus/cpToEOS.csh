#!/bin/tcsh
if ( $1 == "" ) then
    echo "ERROR: Please input work folder name."
    echo "Ex: ./cpToEOS.csh [name]"
    exit
endif
if ( ! ( -e $1 ) ) then
	echo "ERROR: Here is no work folder name $1 "
	exit
endif

cmsenv
set eosDir="eos/cms/store/user/jtsai/TTBarCPV/results"
set subDirs=`\ls -l $1 | grep drwxr-xr-x | awk '{print $9}'` 
foreach subDir( $subDirs )

    set myDir="$1/$subDir"
    set toEOS="$eosDir/$1/$subDir"

    echo ">> In $myDir"
    eos mkdir -p $toEOS

    set roots=`\ls -l "$myDir/output" | grep root | awk '{print $9}'`
    set nroots=`\ls -l "$myDir/output" | grep root | awk '{print $9}' | wc -l `
    echo "   Copy $nroots files to $toEOS..."
    foreach root($roots)
        xrdcp $myDir/output/$root xroot://eoscms.cern.ch//$toEOS/$root
    end
end
