#/bin/tcsh
if ( $1 == "" ) then
    echo "ERROR: Please input work folder name."
    echo "Ex: ./combineDataRuns.csh [name]"
    exit
endif
if ( ! ( -e $1 ) ) then
    echo "ERROR: Here is no work folder name $1 "
    exit
endif

cmsenv
set dataNames='SingleElectron SingleMu'
set runName='Run2012'

cd $1
foreach data( $dataNames )

    echo ">*************************************************************"
    echo ">> Checking $data..."
    set rootfiles=`ls -l | grep $data | grep root | grep $runName | awk '{print $9}'`

    if (  `echo $rootfiles | wc -w` != 0 ) then
        if ( -e $data.root ) then
            rm -f $data.root
        endif
        echo ">> Combining..."
        hadd "$data.root" $rootfiles
        echo ">> [Done] $1/$data.root"
    else
        echo ">> [ERROR] No data $data!"
    endif

end
cd -
