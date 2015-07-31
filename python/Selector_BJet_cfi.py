import FWCore.ParameterSet.Config as cms

from TTBarCPV.TTBarCPVAnalysisRun1.Selector_Jet_cfi import*

defaultBJetSelectionParameters = defaultJetSelectionParameters.clone(
	# jet general selections
	jetType                  = cms.string('BJet'),

	# jet b tag
	jetCombinedSVBJetTagsMin = cms.double(0.679),
)
