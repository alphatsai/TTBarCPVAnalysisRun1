import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
from inputFiles_cfi import * 
from inputJsons_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'SemiLeptanicResultsCheck.root',
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
options.register('ttreedir', 'SemiLeptanic',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of ROOT TTree dir"
    )
options.register('MinChi2Cut',  0,
    VarParsing.multiplicity.singleton,    
    VarParsing.varType.float,
    "Event selection for hadronic top min of min_chi2"
    )
options.register('MaxChi2Cut', 40,
    VarParsing.multiplicity.singleton,    
    VarParsing.varType.float,
    "Event selection for hadronic top max of max_chi2"
    )
options.register('Owrt', 'MT:3',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Weight the obseverble in top mass"
    )
options.register('DoWrtEvt', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do event weight"
    ) 
options.register('Debug', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out Debug info during run"
    )
options.parseArguments()

process = cms.Process("CheckEventsOfLepJets")
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

isSkim = False
if options.ttreedir.lower() == 'skim':
    isSkim = True

process.CheckEventsOfLepJets = cms.EDAnalyzer('SemiLeptanicResultsCheckObs',
    MaxEvents             = cms.int32(options.MaxEvents),
    ReportEvery           = cms.int32(options.reportEvery),  
    InputTTree            = cms.string(options.ttreedir+'/root'),
    InputFiles            = cms.vstring(FileNames), 
    #InputFiles            = cms.vstring(FileNames_SemiLeptAnaTopResults), 
    MaxChi2Cut            = cms.double(options.MaxChi2Cut),
    MinChi2Cut            = cms.double(options.MinChi2Cut),
    Owrt                  = cms.double(Oweight), 
    DoWrtEvt              = cms.bool(options.DoWrtEvt),
    Debug                 = cms.bool(options.Debug),
) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.CheckEventsOfLepJets)

