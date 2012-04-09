import FWCore.ParameterSet.Config as cms

process = cms.Process("CALIB")

from PhysicsTools.PatAlgos.patTemplate_cfg import *
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
'file:/tmp/bbilin/MyOutputFile_1_1_0Jk.root'
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(  'GR_R_44_V15::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


#addJetID(process,cms.InputTag('ak5CaloJets'),'ak5')

removeSpecificPATObjects( process, ['Taus'] )
process.patDefaultSequence.remove( process.patTaus )

removeMCMatching(process, ['All'])
addPfMET(process, 'PF')


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
process.patJets.tagInfoSources  = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAOD"),
    )








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
process.patDefaultSequence *
 process.demo
)
process.out.outputCommands = cms.untracked.vstring('drop *')
