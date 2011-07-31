import FWCore.ParameterSet.Config as cms

process = cms.Process("CALIB")

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('PATSummaryTables')

inputJetCorrLabel = ('AK5PF', ['L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual'])

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/data/Run2010B/Electron/RECO/PromptReco-v2/000/149/294/88FA418D-00E5-DF11-B320-000423D94E70.root'
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(  'START41_V0::All')
process.load("Configuration.StandardSequences.MagneticField_cff")


#addJetID(process,cms.InputTag('ak5CaloJets'),'ak5')

removeMCMatching(process, ['All'])
addPfMET(process, 'PF')




from PhysicsTools.PatAlgos.tools.jetTools import *
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

removeSpecificPATObjects( process, ['Taus'] )
process.patDefaultSequence.remove( process.patTaus )



process.TFileService=cms.Service("TFileService",
    fileName=cms.string("2010_elec_ntuple.root")
)



process.demo = cms.EDAnalyzer("ntupleGenerator",
				PfJetAlg = cms.string("selectedPatJets"),
				elecTag  = cms.InputTag("selectedPatElectrons"),

)


process.p = cms.Path(
#process.elec_HLT *    
process.patDefaultSequence *
process.demo
)
