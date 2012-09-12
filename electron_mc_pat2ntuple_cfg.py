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

inputJetCorrLabel =  ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
inputJetCorrLabell = ('AK5Calo',['L2Relative', 'L3Absolute']) 

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'root://eoscms.cern.ch//eos/cms//store/user/bbilin/ntuples/mc/copy/MyOutputFile_1_1_aS5.root'
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(  'START52_V9::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


#addJetID(process,cms.InputTag('ak5CaloJets'),'ak5')

#removeSpecificPATObjects( process, ['Taus'] )
#process.patDefaultSequence.remove( process.patTaus )


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








process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')


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


process.demo = cms.EDAnalyzer("ntupleGenerator",

				PfJetAlg = cms.string("selectedPatJets"),
				CaloJetAlg =cms.string("selectedPatJetsAK5Calo"),
				elecTag  = cms.InputTag("selectedPatElectrons")
				
)

process.p = cms.Path(
#process.elec_HLT *
#process.hltSelection *    
    process.mvaID + 
    process.patDefaultSequence+
    process.patConversions*
 process.demo
)
process.out.outputCommands = cms.untracked.vstring('drop *')
