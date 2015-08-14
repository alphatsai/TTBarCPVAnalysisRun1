import FWCore.ParameterSet.Config as cms

infiniteMax =  1000000
infiniteMin = -1000000

defaultLeptonSelectionParameters = cms.PSet(
    # Lepton general selections
    lepType                       = cms.string('LooseLepton'),
    lepPtMin                      = cms.double(10),
    lepPtMax                      = cms.double(infiniteMax),
    lepEtMin                      = cms.double(20),
    lepEtMax                      = cms.double(infiniteMax),
    lepAbsEtaMin                  = cms.double(infiniteMin),
    lepAbsEtaMax                  = cms.double(2.5),
    lepAbsEtaExcludeMin           = cms.double(1.4442),
    lepAbsEtaExcludeMax           = cms.double(1.5660),
    lepRelIsoR03Min               = cms.double(infiniteMin),
    lepRelIsoR03Max               = cms.double(0.15),
    lepRelIsoR04Min               = cms.double(infiniteMin),
    lepRelIsoR04Max               = cms.double(0.2),

    # Electron special
    ElAbsTrackDxyPVMin            = cms.double(infiniteMin),
    ElAbsTrackDxyPVMax            = cms.double(infiniteMax),
    EgammaMVATrigMin              = cms.double(0),
    EgammaMVATrigMax              = cms.double(1),

    NumberOfExpectedInnerHitsMin  = cms.double(infiniteMin),
    NumberOfExpectedInnerHitsMax  = cms.double(infiniteMax),

    # Muon special
    MuAbsInnerTrackDxyPVMin       = cms.double(infiniteMin),
    MuAbsInnerTrackDxyPVMax       = cms.double(infiniteMax),
    MuGlobalNormalizedChi2Min     = cms.double(infiniteMin),
    MuGlobalNormalizedChi2Max     = cms.double(infiniteMax),

    MuNMuonhitsMin                = cms.double(infiniteMin),
    MuNMuonhitsMax                = cms.double(infiniteMax),
    MuNMatchedStationsMin         = cms.double(infiniteMin),
    MuNMatchedStationsMax         = cms.double(infiniteMax),
    MuNTrackLayersWMeasurementMin = cms.double(infiniteMin), 
    MuNTrackLayersWMeasurementMax = cms.double(infiniteMax),

    CheckGlobalMuon               = cms.bool(True),
    CheckTrackerMuon              = cms.bool(True)
)
