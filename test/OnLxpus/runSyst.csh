#!/bin/tcsh
if ( $1 == "" ) then
    echo "[INFO] Please input the name of work dir"
    echo "       ./runSyst.csh [myDir]"
    exit
endif
if ( -e $1 ) then
    echo "[ERROR] $1 has been exist!"
    exit
endif

set qName="cmscaf1nd"
set systName="JES JER BTagSF MuonIDSF MuonIsoSF ElectronIDSF PU TopPt"
set shifts="Up Down"
set dataset="datasetListLepJetSystUnc/datasetListLepJetSkimedSamples"
foreach syst($systName)
    foreach shift($shifts)
        echo "[INFO] Sumitting $syst $shift..."
        ./createAndSubmitJobs.py -w $1_"$syst""$shift" -d datasetListLepJetSystUnc/datasetListLepJetSkimedSamples"$syst""$shift".txt -c ../semiLeptanicAnalysis_cfg.py -q $qName
    end
end
