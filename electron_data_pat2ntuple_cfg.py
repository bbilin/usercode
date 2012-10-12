import FWCore.ParameterSet.Config as cms

process = cms.Process("CALIB")

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.pfTools import *
usePFIso( process )
#process.patElectrons.pfElectronSource = 'particleFlow'
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *



mytrigs = ['*']
# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')

inputJetCorrLabel =  ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']) 
inputJetCorrLabell = ('AK5Calo',['L2Relative', 'L3Absolute']) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'root://eoscms.cern.ch//eos/cms/store/data/Run2012A/DoubleElectron/RECO/30May2012-GR_P_V32_525p1-v1/0000/040D970A-B0AA-E111-A70B-0026189438F9.root'
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(  'GR_R_52_V9D::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

removeMCMatching(process, ['All'])
addPfMET(process, 'PF')

process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'



from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5CaloJets'),
		'AK5','Calo',
                 doJTA        = False,
                 doBTagging   = False,
                 jetCorrLabel = inputJetCorrLabell,
                 doType1MET   = True,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )


switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = inputJetCorrLabel,
                 doType1MET   = True,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )

process.patJets.addTagInfos = True
#process.patJets.tagInfoSources  = cms.VInputTag(
#    cms.InputTag("secondaryVertexTagInfosAOD"),
#    )


process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                   numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                   )
 #HB + HE noise filtering
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                         vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                          )


process.TFileService=cms.Service("TFileService",
    fileName=cms.string("elec_ntuple.root")
)



from HLTrigger.HLTfilters.hltHighLevel_cfi import *
if mytrigs is not None :
    process.hltSelection = hltHighLevel.clone(TriggerResultsTag = 'TriggerResults::HLT', HLTPaths = mytrigs)
    process.hltSelection.throw = False


#process.elec_HLT = cms.EDFilter("HLTHighLevel",
#     TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#     HLTPaths = cms.vstring('HLT_Photon10_L1R','HLT_Photon15_Cleaned_L1R','HLT_Ele15_SW_CaloEleId_L1R','HLT_Ele17_SW_CaloEleId_L1R','HLT_Ele17_SW_TightEleId_L1R','HLT_Ele17_SW_TighterEleIdIsol_L1R'),   ## provide list of HLT paths (or patterns) you want
#     eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
#     andOr = cms.bool(True),              # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
#     throw = cms.bool(False)    # throw exception on unknown path names
# ) 


process.patElectrons.electronIDSources = cms.PSet(
    #MVA
    mvaTrigV0 = cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
)

process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )

#add pat conversions
process.patConversions = cms.EDProducer("PATConversionProducer",
    # input collection
    #electronSource = cms.InputTag("gsfElectrons"),
    electronSource = cms.InputTag("selectedPatElectrons")  
    # this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer. ,
)


from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

process.demo = cms.EDAnalyzer("ntupleGenerator",

				PfJetAlg = cms.string("selectedPatJets"),
				CaloJetAlg =cms.string("selectedPatJetsAK5Calo"),
				elecTag  = cms.InputTag("selectedPatElectrons")
				
)
process.p = cms.Path(
#process.elec_HLT *
process.scrapingVeto*
    process.primaryVertexFilter*
    process.HBHENoiseFilter*
#process.hltSelection     * 
    process.mvaID + 
    process.patDefaultSequence+
    process.patConversions*
process.kt6PFJetsForIsolation*
 process.demo
)
process.out.outputCommands = cms.untracked.vstring('drop *')

#process.out.outputCommands +=[
#     'keep *_patConversions*_*_*'
#]
