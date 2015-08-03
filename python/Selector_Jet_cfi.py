import FWCore.ParameterSet.Config as cms

infiniteMax =  1000000
infiniteMin = -1000000

defaultJetSelectionParameters = cms.PSet(
	# jet general selections
	jetType                  = cms.string('TightJet'),
	jetPtMin                 = cms.double(30),
	jetPtMax                 = cms.double(infiniteMax),
	jetAbsEtaMin             = cms.double(infiniteMin),
	jetAbsEtaMax             = cms.double(2.4),

	# jet id
	jetNConstituentsMin      = cms.double(1),
	jetNConstituentsMax      = cms.double(infiniteMax),
	jetCHFMin                = cms.double(0),
	jetCHFMax                = cms.double(infiniteMax),
	jetNCHMin                = cms.double(0),
	jetNCHMax                = cms.double(infiniteMax),
	jetNEFMin                = cms.double(infiniteMin),
	jetNEFMax                = cms.double(0.99),
	jetNHFMin                = cms.double(infiniteMin),
	jetNHFMax                = cms.double(0.99),
	jetCEFMin                = cms.double(infiniteMin),
	jetCEFMax                = cms.double(0.99),

	# jet b tag
	jetCombinedSVBJetTagsMin = cms.double(infiniteMin),
	jetCombinedSVBJetTagsMax = cms.double(infiniteMax),
)
