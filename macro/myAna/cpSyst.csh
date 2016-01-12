#!/bin/tcsh
if ( $1 == "" ) then
    echo "[INFO] Please add the dir name"
    echo "       E.x: ./cpSyst.csh 11Jan_Syst"
    exit
endif

set systName="PU JER BTagSF TopPT elID muID muISO"
set tuneName="up down"
set pathName="/afs/cern.ch/work/j/jtsai/myAna/TTBarCPV/CMSSW_7_2_1_dev/src/TTBarCPV/TTBarCPVAnalysisRun1/test/OnLxpus"

foreach syst($systName)
    foreach tune($tuneName)
        set dir=$syst$tune
        echo "[INFO] Copying from $1_$dir..."
        if ( ! -e $dir ) then
            mkdir $dir
        endif 
        scp jtsai@lxplus.cern.ch:$pathName"/"$1"_"$dir"/"Final_histograms_SemiLeptanic.root $dir"/"
    end
end
