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
options.register('dRMatching', 1000,
    VarParsing.multiplicity.singleton,    
    VarParsing.varType.float,
    "Matching object and gen particle with deltaR( lepton, jet )"
    )
options.register('dRIsoLeptonFromJets', 0.5,
    VarParsing.multiplicity.singleton,    
    VarParsing.varType.float,
    "isolate lepton with deltaR( lepton, jet )"
    )
options.register('Owrt', 'MT:3',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Weight the obseverble in top mass"
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

from TTBarCPV.TTBarCPVAnalysisRun1.Selector_Vertex_cfi   import*
from TTBarCPV.TTBarCPVAnalysisRun1.Selector_Jet_cfi      import*
from TTBarCPV.TTBarCPVAnalysisRun1.Selector_BJet_cfi     import*
from TTBarCPV.TTBarCPVAnalysisRun1.Selector_nonBJet_cfi  import*
from TTBarCPV.TTBarCPVAnalysisRun1.Selector_Lepton_cfi   import*
from TTBarCPV.TTBarCPVAnalysisRun1.Selector_Muon_cfi     import*
from TTBarCPV.TTBarCPVAnalysisRun1.Selector_Electron_cfi import*

process.CheckEventsOfLepJets = cms.EDAnalyzer('SemiLeptanicResultsCheck',
    MaxEvents             = cms.int32(options.MaxEvents),
    ReportEvery           = cms.int32(options.reportEvery),  
    InputTTree            = cms.string(options.ttreedir+'/root'),
    #InputFiles            = cms.vstring(FileNames), 
    #InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/TTJets_SemiLeptMGDecays.root'), 
    InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/TTJets_FullLeptMGDecays.root'), 
    #InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/T_t-channel.root'), 
    #InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/Tbar_t-channel.root'), 
    #InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/Tbar_t-channel.root'), 
    #InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/WJetsToLNu.root'), 
    #InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/WZ.root'), 
    #InputFiles            = cms.vstring('file:OnLxpus/03Aug_MC_IsoLeptonFromJets_Tree/WW.root'), 
    InputJsons            = cms.vstring(JsonNames), 
    HLT_MuChannel         = cms.vint32( 2868,3244,3542,4204,4205,4827,5106,5573  ), # HLT_IsoMu24_eta2p1_v*
    HLT_ElChannel         = cms.vint32( 3155,3496,4002,4003,4004,5043 ),            # HLT_Ele27_WP80_v* 
    SelPars_Vertex        = defaultVertexSelectionParameters.clone(), 
    SelPars_Jet           = defaultJetSelectionParameters.clone(), 
    SelPars_BJet          = defaultBJetSelectionParameters.clone(), 
    SelPars_NonBJet       = defaultNonBJetSelectionParameters.clone(), 
    SelPars_LooseLepton   = defaultLeptonSelectionParameters.clone(), 
    SelPars_TightMuon     = defaultMounSelectionParameters.clone(), 
    SelPars_TightElectron = defaultElectronSelectionParameters.clone(),
    dR_Matching           = cms.double(options.dRMatching),
    dR_IsoLeptonFromJets  = cms.double(options.dRIsoLeptonFromJets),
    Owrt                  = cms.double(Oweight), 
    Debug                 = cms.bool(options.Debug),
) 

process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(process.CheckEventsOfLepJets)

