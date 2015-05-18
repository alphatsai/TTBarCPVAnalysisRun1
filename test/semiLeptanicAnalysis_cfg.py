import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'SemiLeptanicAnalysis.root',
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
options.register('NJets', 3,
	VarParsing.multiplicity.singleton,	
	VarParsing.varType.int,
	"Number of jets"
	)
options.register('NonBjetCSVThr', 0.679,
	VarParsing.multiplicity.singleton,	
	VarParsing.varType.float,
	"Non B jet CSV threshold"
	)
options.register('IsoElePt', 45.,
	VarParsing.multiplicity.singleton,	
	VarParsing.varType.float,
	"Pt of isolated electron"
	)
options.register('IsoMuonPt', 35.,
	VarParsing.multiplicity.singleton,	
	VarParsing.varType.float,
	"Pt of isolated muon"
	)
options.register('Owrt', 'MT:3',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"Weight the obseverble in top mass"
	)
options.register('Debug', True,
	VarParsing.multiplicity.singleton,
	VarParsing.varType.bool,
	"Print out Debug info during run"
	)
options.parseArguments()

process = cms.Process("SemiLeptanic")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) ) 

process.source = cms.Source("EmptySource")
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(options.outFilename) 
	)

MT = 173.0
Oweight = 1.0
power = int(options.Owrt.split(':')[1])
i = 0
while ( i < power ):
	Oweight *= MT
	i+=1

process.SemiLeptanic = cms.EDAnalyzer('SemiLeptanicAnalysis',
	MaxEvents   = cms.int32(options.MaxEvents),
	ReportEvery = cms.int32(options.reportEvery),  
	InputTTree  = cms.string(options.ttreedir+'/root'),
	#InputFiles  = cms.vstring(FileNames), 
	InputFiles  = cms.vstring(FileNames_BprimtKits_SemiLeptTest),
	#InputFiles  = cms.vstring(FileNames_BprimtKits_SemiLept),
	MuonHLT     = cms.vint32( 2868,3244,3542,4204,4205,4827,5106,5573  ), #HLT_IsoMu24_eta2p1_v*
	ElectronHLT = cms.vint32( 3155,3496,4002,4003,4004,5043 ), #HLT_Ele27_WP80_v 
	#MuonHLT     = cms.vint32( 1169,1170,1812,1813,2283,2642,2909,3310,3587,4873,5195,5665 ), #HLT_Mu30_v1
	#ElectronHLT = cms.vint32( 724,725,973,728,975,1475,1476,2132,2502,2819,1478 ), #HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*
	#ElectronHLT = cms.vint32( 724,725,973,1475,1478 ),
	NJets = cms.double(options.NJets),
	IsoElePt = cms.double(options.IsoElePt),
	IsoMuonPt = cms.double(options.IsoMuonPt),
	NonBjetCSVThr = cms.double(options.NonBjetCSVThr), 
	Owrt    = cms.double(Oweight), 
	Debug    = cms.bool(options.Debug) 
	) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.SemiLeptanic)

