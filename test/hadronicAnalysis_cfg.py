import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from inputFiles_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'HadronicAnalysis.root',
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
options.register('ttreedir', 'ntuple',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"Name of ROOT TTree dir: Either 'ntuple' or 'skim' or 'bVeto'"
	)
#options.register('JetPtMin', 50.,
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.float,
#    "Minimum b jet Pt"
#    )
#options.register('doPUReweighting', True,
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.bool,
#    "Do pileup reweighting"
#)
options.parseArguments()

process = cms.Process("Hadronic")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.MaxEvents) ) 

process.source = cms.Source("EmptySource")
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(options.outFilename) 
	)

process.BprimebH = cms.EDAnalyzer('HadronicAnalysis',
	MaxEvents           = cms.int32(options.MaxEvents),
	ReportEvery         = cms.int32(options.reportEvery),  
	InputTTree          = cms.string(options.ttreedir+'/tree'),
	InputFiles          = cms.vstring(FileNames), 
	) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.Hadronic)

