import FWCore.ParameterSet.Config as cms

infiniteMax =  1000000
infiniteMin = -1000000

defaultJetSelectionParameters = cms.PSet(
    # https://twiki.cern.ch/twiki/bin/view/CMS/TopJMERun1#Jets
    # jet general selections
    jetType                  = cms.string('LooseJet'), # LooseJet, TightJet
    jetPtMin                 = cms.double(30),
    jetPtMax                 = cms.double(infiniteMax),
    jetAbsEtaMin             = cms.double(infiniteMin),
    jetAbsEtaMax             = cms.double(2.4),

    # jet id
    jetNConstituentsMin      = cms.double(1),
    jetNConstituentsMax      = cms.double(infiniteMax),
    jetNEFMin                = cms.double(infiniteMin),
    jetNEFMax                = cms.double(0.99),         # Loose <0.99, Tight <0.90
    jetNHFMin                = cms.double(infiniteMin),
    jetNHFMax                = cms.double(0.99),         # Loose <0.99, Tight <0.90
    # |eta| <= 2.4
    jetNCHMin                = cms.double(0),
    jetNCHMax                = cms.double(infiniteMax),
    jetCHFMin                = cms.double(0),
    jetCHFMax                = cms.double(infiniteMax),
    jetCEFMin                = cms.double(infiniteMin),
    jetCEFMax                = cms.double(0.99),        # Loose <0.99, Tight <0.90

    # jet b tag
    jetCombinedSVBJetTagsMin = cms.double(infiniteMin),
    jetCombinedSVBJetTagsMax = cms.double(infiniteMax),
)
