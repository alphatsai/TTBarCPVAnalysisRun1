import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
from inputFiles_cfi import * 
from inputJsons_cfi import * 

options = VarParsing('python')

options.register('outFilename', 'PDFUncertaintyAna.root',
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
options.register('PdfName', 'CT10nnlo.LHgrid', # ls /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/lhapdf/5.9.1-cms4/share/lhapdf/PDFsets
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "PDF method name"
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

process.CheckEventsOfLepJets = cms.EDAnalyzer('PDFUncertaintyAna',
    MaxEvents   = cms.int32(options.MaxEvents),
    ReportEvery = cms.int32(options.reportEvery),  
    InputTTree  = cms.string(options.ttreedir+'/pdftree'),
    InputFiles  = cms.vstring(FileNames), 
    #InputFiles  = cms.vstring("root://eoscms//eos/cms/store/user/jtsai/TTBarCPV/results/08Feb_PDFTree/TTJets_SemiLeptMGDecays/SemiLeptanicAnalysis_35.root"), 
    MaxChi2Cut  = cms.double(options.MaxChi2Cut),
    MinChi2Cut  = cms.double(options.MinChi2Cut),
    PdfName     = cms.string(options.PdfName), 
    Debug       = cms.bool(options.Debug),
) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.CheckEventsOfLepJets)

