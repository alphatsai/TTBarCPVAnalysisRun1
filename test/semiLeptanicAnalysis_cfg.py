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
options.register('Owrt', 'MT:3',
	VarParsing.multiplicity.singleton,
	VarParsing.varType.string,
	"Weight the obseverble in top mass"
	)
#options.register('JetPtMin', 50.,
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.float,
#    "Minimum b jet Pt"
#    )
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
	InputFiles  = cms.vstring(FileNames_BprimtKits_SemiLept),
	Owrt    = cms.double(Oweight), 
	Debug    = cms.bool(options.Debug) 
	) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.SemiLeptanic)

