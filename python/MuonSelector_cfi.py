import FWCore.ParameterSet.Config as cms

from TTBarCPV.TTBarCPVAnalysisRun1.LeptonSelector_cfi import*

defaultMounSelectionParameters = defaultLeptonSelectionParameters.clone(
	# Lepton general selections
	lepType                       = cms.string('TightMuon'),
	lepPtMin                      = cms.double(26),
	lepAbsEtaMax                  = cms.double(2.1),
	lepRelIsoR04Max               = cms.double(0.12),

	# Muon special
	MuAbsInnerTrackDxyPVMax       = cms.double(0.2),
	MuGlobalNormalizedChi2Max     = cms.double(10),

	MuNMuonhitsMin                = cms.double(0),
	MuNMatchedStationsMin         = cms.double(1),
	MuNTrackLayersWMeasurementMin = cms.double(5), 

	CheckGlobalMuon               = cms.bool(True),
	CheckTrackerMuon              = cms.bool(False)
)
