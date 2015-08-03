import FWCore.ParameterSet.Config as cms

infiniteMax =  1000000
infiniteMin = -1000000

defaultVertexSelectionParameters = cms.PSet(
	# Vertex general selections
	vtxType          = cms.string('GoodVertex'),
	vtxNdofMin       = cms.double(4),
	vtxNdofMax       = cms.double(infiniteMax),
	vtxAbsZMin       = cms.double(infiniteMin),
	vtxAbsZMax       = cms.double(24),
	vtxRhoMin        = cms.double(infiniteMin),
	vtxRhoMax        = cms.double(2),

	CheckIsOfflinePV = cms.bool(True),
	CheckIsFake      = cms.bool(True)
)
