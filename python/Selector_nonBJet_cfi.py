import FWCore.ParameterSet.Config as cms

from TTBarCPV.TTBarCPVAnalysisRun1.Selector_Jet_cfi import*

defaultNonBJetSelectionParameters = defaultJetSelectionParameters.clone(
	# jet general selections
	jetType                  = cms.string('NonBJet'),

	# jet b tag
	jetCombinedSVBJetTagsMax = cms.double(0.679),
)
