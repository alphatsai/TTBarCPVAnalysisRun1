import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
from inputFiles_cfi import * 
from inputJsons_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'Skim.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
    )
options.register('MaxEvents', -1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Run events max"
    )
options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
    )
options.register('ttreedir', 'bprimeKit',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of ROOT TTree dir: Either 'ntuple' or 'skim' or 'bVeto'"
    )
#options.register('NJet', 4,
options.register('NJet', 3,
    VarParsing.multiplicity.singleton, 
    VarParsing.varType.int,
    "Number of jets"
    )
options.register('NLep', 1,
    VarParsing.multiplicity.singleton, 
    VarParsing.varType.int,
    "Number of leptons"
    )
options.register('LepPt', 10,
    VarParsing.multiplicity.singleton, 
    VarParsing.varType.float,
    "Pt of lepton"
    )
options.register('JetPt', 20,
    VarParsing.multiplicity.singleton, 
    VarParsing.varType.float,
    "Pt of lepton"
    )
options.register('Debug', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out Debug info during run"
    )
options.parseArguments()

process = cms.Process("Skim")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) 

process.source = cms.Source("EmptySource")
process.TFileService = cms.Service("TFileService",
 fileName = cms.string(options.outFilename) 
 )

process.Skim = cms.EDAnalyzer('LeptonJetFilter',
    MaxEvents   = cms.int32(options.MaxEvents),
    ReportEvery = cms.int32(options.reportEvery),  
    InputTTree  = cms.string(options.ttreedir+'/root'),
    InputJsons  = cms.vstring(JsonNames), 
    InputFiles  = cms.vstring(FileNames), 
    #InputFiles  = cms.vstring(FileNames_BprimtKits_SemiLeptTest),
    #InputFiles  = cms.vstring(FileNames_BprimtKits_SemiLept),
    MuonHLT     = cms.vint32( 2868,3244,3542,4204,4205,4827,5106,5573  ), #HLT_IsoMu24_eta2p1_v*
    ElectronHLT = cms.vint32( 3155,3496,4002,4003,4004,5043 ), #HLT_Ele27_WP80_v 
    NJet  = cms.int32(options.NJet),
    NLep  = cms.int32(options.NLep),
    LepPt = cms.double(options.LepPt),
    JetPt = cms.double(options.JetPt),
    Debug = cms.bool(options.Debug) 
) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.Skim)

