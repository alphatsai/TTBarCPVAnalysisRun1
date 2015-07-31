import FWCore.ParameterSet.Config as cms

from TTBarCPV.TTBarCPVAnalysisRun1.LeptonSelector_cfi import*

defaultElectronSelectionParameters = defaultLeptonSelectionParameters.clone(
	# Lepton general selections
	lepType                       = cms.string('TightElectron'),
	lepEtMin                      = cms.double(30),
	lepAbsEtaMax                  = cms.double(2.1),
	lepRelIsoR03Max               = cms.double(0.1),

	# Electron special
	ElAbsTrackDxyPVMax            = cms.double(0.02),
	EgammaMVATrigMin              = cms.double(0.5),
	EgammaMVATrigMax              = cms.double(infiniteMax),
	NumberOfExpectedInnerHitsMax  = cms.double(1),

)
