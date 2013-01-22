
// -*- C++ -*-
//
// Package:    Fwdcalib
// Class:      Fwdcalib
// 
/**\class Fwdcalib Fwdcalib.cc UserCode/Fwdcalib/src/Fwdcalib.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Bugra Bilin,8 R-004,+41227676479,
//         Created:  Tue May  3 16:39:40 CEST 2011
// $Id: FwdCalib.cc,v 1.11 2012/10/12 10:00:01 bbilin Exp $
//
//



//
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/EventSetupRecord.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "FWCore/Framework/interface/Run.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TROOT.h>
#include "TF1.h"
#include "TMath.h"
#include <TSystem.h>
#include "TFile.h"
#include <TCanvas.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <Math/VectorUtil.h>
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "PhysicsTools/PatAlgos/plugins/PATConversionProducer.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//elec

#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

//jets
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
// MET
#include "DataFormats/PatCandidates/interface/MET.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

//pf
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
//vtx
// not sure if all that includes are needed .....
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h" 
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TTree.h"


#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include <map>
#include "DataFormats/Common/interface/AssociationMap.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
class TH1F;
class TH2F;
class TStyle;
class TTree;

//
// class decleration
//

using namespace edm;
using namespace reco;
using namespace std;
using namespace ROOT::Math::VectorUtil;
using namespace HepMC;
using namespace isodeposit;


class ntupleGenerator : public edm::EDAnalyzer {
  public:
    explicit ntupleGenerator(const edm::ParameterSet&);
    ~ntupleGenerator();

    const edm::ValueMap<double>& getValueMap(const edm::Event& iEvent,edm::InputTag& inputTag);
      
    //void produce(const edm::Handle<SuperClusterCollection>& SuperClusters,
    //				      const HFEMClusterShapeAssociationCollection& AssocShapes,
    //	     RecoEcalCandidateCollection& RecoECand);


  private:
    //  virtual void beginJob(const edm::EventSetup&) ;
    virtual void beginJob();	
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------

    float DeltaPhi(float phi1, float phi2);
    float DeltaR(float eta1, float eta2, float phi1, float phi2);
    edm::InputTag elecTag_;
    string PfJetAlg;
    edm::Service<TFileService> fs;
    TTree * myTree;
    int event, run,lumi,bxnumber,realdata;
    int hlt_trigger_fired;
    float ElecPt[50],ElecE[50],ElecM[50],ElecPx[50], ElecPy[50],ElecPz[50],ElecEta[50],ElecPhi[50],ElecGsfTrk_d0[50],Elecdr03TkSumPt[50],Elecdr03EcalRecHitSumEt[50],Elecdr03HcalTowerSumEt[50],ElecscSigmaIEtaIEta[50],ElecdeltaPhiSuperClusterTrackAtVtx[50],ElecdeltaEtaSuperClusterTrackAtVtx[50],ElechadronicOverEm[50],ElecgsfTrack_numberOfLostHits[50],ElecDcotTheta[50],ElecMVATrigId[50],ElecMVANonTrigId[50],Elecd0vtx[50],Elecdzvtx[50], chIso03[50],nhIso03[50],phIso03[50],puChIso03[50],relIso[50],relIsodb[50],relIsorho[50],fMVAVar_fbrem[50],  fMVAVar_kfchi2[50], fMVAVar_kfhits[50],fMVAVar_gsfchi2[50],fMVAVar_detacalo[50],fMVAVar_see[50], fMVAVar_spp[50] , fMVAVar_etawidth[50],fMVAVar_phiwidth[50] ,fMVAVar_e1x5e5x5[50] , fMVAVar_R9[50] ,fMVAVar_EoP[50]  ,fMVAVar_IoEmIoP[50],fMVAVar_eleEoPout[50] ,fMVAVar_PreShowerOverRaw[50];
    int Elecindex, HFElecindex,ElecCharge[50],ElecIsEB[50], ElecIsEE[50], Elecooeoop[50]; 
    bool hasMatchedConversion[50]; 
    float hfelec_L9[50] ,hfelec_S9[50] ,hfelec_L25[50] ,hfelec_L1[50] ,hhfclustereta[50]   ,hhfclusterphi[50], hfelec_Pt[50] ,hfelec_Eta[50], hfelec_Px[50],hfelec_Py[50],hfelec_Pz[50],hfelec_Phi[50],hfelec_E[50],hfelec_M[50] ,hfelec_charge[50];
    //particle information
    int par_index, mom[500], daug[500];
    float ParticlePt[500], ParticleEta[500], ParticlePhi[500], ParticlePx[500], ParticlePy[500], ParticlePz[500], ParticleE[500], ParticleM[500];
    int ParticleId[500], ParticleStatus[500], ParticleMother[500][100], ParticleDaughter[500][100];
    //reco jets 
    int pfNjets,pfNtracks;
    int CaloNjets;
    float PFjet_mva[100], PFjetEta[100], PFjetPhi[100],PFjetPt[100],PFCorrjetPt[100],PFjetCEMF[100],PFjetNEMF[100],PFjetPx[100],PFjetPy[100],PFjetPz[100],PFjetE[100],PFjetM[100],PFHadEHF[100],PFEmEHF[100];
    float PFjetTrkVZ[100][100],PFjetTrkPT[100][100];
    float vtxZ[100],vtxZerr[100];
    float vtxY[100],vtxYerr[100];
    float vtxX[100],vtxXerr[100];
    int vtxisValid[100],vtxisFake[100];
    int nVertices,nGoodVertices; 
    int NumInteractions,BunchCrossing;
    float TrueNumInteractions,Tnpv_check;
    double MyWeight;
    int techTrigger[44];
    //met 
    float  pfMET, pfSET;
    float  pfMETX;
    float  pfMETY;
    float muCorrMET, muCorrSET;
    int tr_1, tr_2, tr_3;
    //  reco::helper::JetIDHelper *jetID;

};



ntupleGenerator::ntupleGenerator(const edm::ParameterSet& iConfig)
{   
  elecTag_ = iConfig.getParameter<edm::InputTag>("elecTag");
  PfJetAlg = iConfig.getParameter<string>("PfJetAlg");
}

ntupleGenerator::~ntupleGenerator()
{
}


// ------------ method called to for each event  ------------
void
ntupleGenerator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Handle<GenParticleCollection> genParticles_h;
  iEvent.getByLabel("genParticles", genParticles_h);
  const GenParticleCollection* genParticles  = genParticles_h.failedToGet () ? 0 : &*genParticles_h;


  Handle<reco::RecoEcalCandidateCollection> HFElectrons;
  iEvent.getByLabel("hfRecoEcalCandidate",HFElectrons);
  const RecoEcalCandidateCollection* hfelec = HFElectrons.failedToGet() ? 0 : &*HFElectrons;



  Handle<reco::HFEMClusterShapeAssociationCollection>  electronHFClusterAssociation;
  iEvent.getByLabel(edm::InputTag("hfEMClusters"),electronHFClusterAssociation);
 
  
  Handle< edm::View<pat::MET> > pfMEThandle;
  iEvent.getByLabel("patMETsPF", pfMEThandle); 

  Handle<ValueMap<float> > puJetIDMVA;
  iEvent.getByLabel(InputTag("puJetMva" ,"fullDiscriminant"),puJetIDMVA);


  edm::Handle<edm::View<pat::Conversion>> conversions;
  iEvent.getByLabel("patConversions",conversions);

  edm::Handle<double>  rho_;
//  iEvent.getByLabel(InputTag("kt6PFJetsForIsolation", "rho"), rho_);
iEvent.getByLabel(InputTag("kt6PFJets", "rho"), rho_);
  //double rhoIso = *(rho_.product());
  double rhoIso = *rho_;
  Handle<edm::View<pat::Jet> > PFjets;
  iEvent.getByLabel(PfJetAlg,PFjets);
  //const edm::View<pat::Jet> & calo_jet = *PFjets;

  Handle<edm::View<pat::Electron> > elecs_h;	
  iEvent.getByLabel(elecTag_,elecs_h);
  const edm::View<pat::Electron>* elec = elecs_h.failedToGet() ? 0 : &*elecs_h;

  // EcalClusterLazyTools const &myEcalCluster;

  Handle<BeamSpot> beamSpotHandle;
 
  if (!iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle)) {
    
	cout<< ">>> No beam spot found !!!" <<endl;
    return;
  }
  
  Handle<VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);  

  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.luminosityBlock();
  bxnumber = iEvent.bunchCrossing();
  realdata = iEvent.isRealData();


  if(!realdata){
    //cout<<"adasdasdasasasdas"<<endl;
    edm::LumiReWeighting LumiWeights_;

    LumiWeights_ = edm::LumiReWeighting("mc_pup.root", "data.root", "h1D", "pup");
 
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;

    NumInteractions=-99;
    TrueNumInteractions=-99;
    BunchCrossing=-99;
    float Tnpv = -1; 
    // (then, for example, you can do)

	for( PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

	  NumInteractions = PVI->getPU_NumInteractions();
	  BunchCrossing = PVI->getBunchCrossing();
	  TrueNumInteractions = PVI->getTrueNumInteractions();
    //	std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
      int BX = PVI->getBunchCrossing();
           
      if(BX == 0) { 
        Tnpv = PVI->getTrueNumInteractions();
        Tnpv_check = Tnpv;
        continue;
      }

    }

    MyWeight = LumiWeights_.weight(Tnpv);
    //std::cout <<MyWeight<<endl;
 
  }


  nVertices = pvHandle->size(); 

  const VertexCollection vertexColl = *(pvHandle.product());

  nGoodVertices = 0;
  
  for (VertexCollection::const_iterator vtx = vertexColl.begin(); vtx!= vertexColl.end(); ++vtx){
    if (vtx->isValid() && !vtx->isFake()) ++nGoodVertices;
  }

  //  if (nVertices != nGoodVertices) cout<<"@@@@@@@@@@@@@   ------> "<<nVertices<<"  "<<nGoodVertices<<endl;

  reco::VertexCollection vtxs = *pvHandle;
  for (int iv = 0;iv<nVertices;iv++){
    vtxX[iv] = vtxs[iv].x();
    vtxY[iv] = vtxs[iv].y();
    vtxZ[iv] = vtxs[iv].z();
    vtxXerr[iv] = vtxs[iv].xError();
    vtxYerr[iv] = vtxs[iv].yError();
    vtxZerr[iv] = vtxs[iv].zError();
    vtxisValid[iv] = vtxs[iv].isValid();
    vtxisFake[iv] = vtxs[iv].isFake();
  }

  //trigger primitive version//
  tr_1=0;
  tr_2=0;
  tr_3=0;

  //if(realdata){

 
  int ntrigs;
  vector<string> trigname;
  vector<bool> trigaccept;
  edm::Handle< edm::TriggerResults > HLTResHandle;
  edm::InputTag HLTTag = edm::InputTag( "TriggerResults", "", "HLT");
  iEvent.getByLabel(HLTTag, HLTResHandle);


  if ( HLTResHandle.isValid() && !HLTResHandle.failedToGet() ) {
    edm::RefProd<edm::TriggerNames> trigNames( &(iEvent.triggerNames( *HLTResHandle )) );
    ntrigs = (int)trigNames->size();

    for (int i = 0; i < ntrigs; i++) {
      trigname.push_back(trigNames->triggerName(i));
      trigaccept.push_back(HLTResHandle->accept(i));
      //cout<<trigname[i]<<endl;
      if (trigaccept[i]){

        if(std::string(trigname[i]).find("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL")!=std::string::npos)tr_1=1;

        if(std::string(trigname[i]).find("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL")!=std::string::npos)tr_2=1;

        if(std::string(trigname[i]).find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL")!=std::string::npos)tr_3=1;

      }
    }
  }

//}
//if(tr_3!=0)cout<<tr_3<<endl;
 
  //----------------------------------------------------------------------------


  //--------------------met-------------------------------------------------
 
  pfMET = (pfMEThandle->front()).et();
  pfSET = (pfMEThandle->front()).sumEt(); 


  pfMETX = (pfMEThandle->front()).px();
  pfMETY = (pfMEThandle->front()).py();

 // muCorrMET      = (muCorrMEThandle->front()).et();
 // muCorrSET      = (muCorrMEThandle->front()).sumEt(); 
  //------------------------------------------------------------------------

  par_index = 0;
  for (int i=0;i<500;i++){
    ParticlePt[i] = -999; 
    ParticleEta[i] = -999; 
    ParticlePhi[i] = -999; 
    ParticlePx[i] = -999; 
    ParticlePy[i] = -999;
    ParticlePz[i] = -999; 
    ParticleE[i] = -999;
    ParticleM[i] = -999;
    ParticleId[i] = -999; 
    ParticleStatus[i] = -999; 
    for (int j=0;j<100;j++){
      ParticleMother[i][j] = -999; 
      ParticleDaughter[i][j] = -999;
    }
  }

  //--------------------particles-------------------------------------------
  
  if (!realdata && genParticles){	
    par_index = 0;
    for (size_t i=0; i<genParticles->size(); ++i){
      //const GenParticle & p = (*genParticles)[i];
	  const Candidate & p = (*genParticles)[i];
	  int id = p.pdgId();
	  int st = p.status();
	  if (st==3 || (st==1 && (abs(id)==13 || abs(id)==11))){
		ParticlePt[par_index] = p.pt();
		ParticleEta[par_index] = p.eta();
		ParticlePhi[par_index] = p.phi();
		ParticlePx[par_index] = p.px();
		ParticlePy[par_index] = p.py();
		ParticlePz[par_index] = p.pz();
		ParticleE[par_index] =p.energy();
		ParticleM[par_index] =p.mass();
		ParticleId[par_index] = id;
		ParticleStatus[par_index] = st;
		++par_index;
	  }
	}   
  }



  for (int jj = 0; jj<50; ++jj){
    Elecd0vtx[jj] = -99.;
    Elecdzvtx[jj] = -99.;
    ElecPt[jj] = -99.;
    ElecPx[jj] = -99.;
    ElecPy[jj] = -99.;
    ElecPz[jj] = -99.;
    ElecEta[jj] = -99.;
    ElecPhi[jj] = -99.;
    ElecE[jj] = -99.;
    ElecM[jj] = -99.;
    ElecCharge[jj] = -99.;
    ElecGsfTrk_d0[jj] = -99.;
    ElecIsEB[jj] = -99.;
    ElecIsEE[jj] = -99.;
    Elecdr03TkSumPt[jj] = -99.;
    Elecdr03EcalRecHitSumEt[jj] = -99.;
    Elecdr03HcalTowerSumEt[jj] = -99.;
    ElecscSigmaIEtaIEta[jj] = -99.;
    ElecdeltaPhiSuperClusterTrackAtVtx[jj] = -99.;
    ElecdeltaEtaSuperClusterTrackAtVtx[jj] = -99.;
    ElechadronicOverEm[jj] = -99.;
    ElecgsfTrack_numberOfLostHits[jj] = -99.;
    ElecDcotTheta[jj] = -99.;
    ElecMVATrigId[jj]=-99.;
    ElecMVANonTrigId[jj]=-99.;
    Elecd0vtx[jj]=-99.;
    Elecdzvtx[jj]=-99.; 
    chIso03[jj]=-99.;
    nhIso03[jj]=-99.;
    phIso03[jj]=-99.;
    puChIso03[jj]=-99.;
    relIso[jj]=-99.;
    relIsodb[jj]=-99.;
    relIsorho[jj]=-99.;
    fMVAVar_fbrem[jj]=-99;
    fMVAVar_kfchi2[jj]=-99; 
    fMVAVar_kfhits[jj]=-99;
    fMVAVar_gsfchi2[jj]=-99;
    fMVAVar_detacalo[jj]=-99;
    fMVAVar_see[jj]=-99; 
    fMVAVar_spp[jj] =-99; 
    fMVAVar_etawidth[jj]=-99;
    fMVAVar_phiwidth[jj] =-99;
    fMVAVar_e1x5e5x5[jj] =-99; 
    fMVAVar_R9[jj] =-99;
    fMVAVar_EoP[jj]  =-99;
    fMVAVar_IoEmIoP[jj]=-99;
    fMVAVar_eleEoPout[jj] =-99;
    fMVAVar_PreShowerOverRaw[jj]=-99;
  }
  Elecindex = 0;

  double probMin =1e-6; 
  double lxyMin=2.0;  
  double nHitsBeforeVtxMax=0;

  //for (unsigned int i=0; i < elec->size()-1;++i){
  // pat::Electron el = elec->at(i);

  for(edm::View<pat::Electron>::const_iterator el = elec->begin(); el != elec->end();  ++el) {
  
    // apply to neutrals

     Elecooeoop[Elecindex] =1.0/ el->ecalEnergy() - (el->eSuperClusterOverP()/el->ecalEnergy());

     if (pvHandle->size() > 0) {
        reco::VertexRef vtx(pvHandle, 0);    
        Elecd0vtx[Elecindex] = el->gsfTrack()->dxy(vtx->position());
        Elecdzvtx[Elecindex] = el->gsfTrack()->dz(vtx->position());
    }
	 else {
        Elecd0vtx[Elecindex] = el->gsfTrack()->dxy();
        Elecdzvtx[Elecindex] = el->gsfTrack()->dz();
    }

	hasMatchedConversion[Elecindex] = false;
    for (edm::View<pat::Conversion>::const_iterator iconv=conversions->begin() ; iconv != conversions->end(); ++iconv) {

      if (iconv->index() != Elecindex) continue;
      if (!(iconv->vtxProb()>probMin &&  iconv->lxy()>lxyMin && iconv->nHitsMax()<=nHitsBeforeVtxMax)) continue;
      hasMatchedConversion[Elecindex] = true;
    }


    const string mvaTrigV0 = "mvaTrigV0";
    const string mvaNonTrigV0 = "mvaNonTrigV0";
    //cout<<el->r9()<<endl;
    //cout<<el->electronID(mvaTrigV0)<<endl;
	ElecMVATrigId[Elecindex]= el->electronID(mvaTrigV0);
	ElecMVANonTrigId[Elecindex]= el->electronID(mvaNonTrigV0);
    //cout<<el->electronID(mvaTrigV0)<<"  "<<el->electronID(mvaNonTrigV0)<<endl;
	ElecPt[Elecindex] = el->pt();
//	cout<<ElecPt[Elecindex]<<"  "<<endl;
	ElecPx[Elecindex] = el->px();
	ElecPy[Elecindex] = el->py();
	ElecPz[Elecindex] = el->pz();
	ElecEta[Elecindex] = el->superCluster()->eta();
	ElecPhi[Elecindex] = el->phi();	
	ElecE[Elecindex] = el->energy();
    //cout<<run<<"  "<<el->superCluster()->eta()<<" "<<el->r9()<<"  "<<el->energy()<<"  "<<el->ecalEnergy()<<"  "<<el->ecalEnergy()/el->energy()<<endl;
	ElecM[Elecindex] = el->mass();
	ElecCharge[Elecindex] = el->charge();
    //	ElecGsfTrk_d0[Elecindex] = el->gsfTrack()->d0();
	ElecGsfTrk_d0[Elecindex] =  el->convDist(); //new (28 02 2012)
	ElecDcotTheta[Elecindex] = el->convDcot();//new (28 02 2012)
	ElecIsEB[Elecindex] = el->isEB();
	ElecIsEE[Elecindex] = el->isEE();
	Elecdr03TkSumPt[Elecindex] = el->dr03TkSumPt();
	Elecdr03EcalRecHitSumEt[Elecindex] = el->dr03EcalRecHitSumEt();
	Elecdr03HcalTowerSumEt[Elecindex] = el->dr03HcalTowerSumEt();
	ElecscSigmaIEtaIEta[Elecindex] = el->scSigmaIEtaIEta();
	ElecdeltaPhiSuperClusterTrackAtVtx[Elecindex] = el->deltaPhiSuperClusterTrackAtVtx();
	ElecdeltaEtaSuperClusterTrackAtVtx[Elecindex] = el->deltaEtaSuperClusterTrackAtVtx();
	ElechadronicOverEm[Elecindex] = el->hadronicOverEm();
    //	ElecgsfTrack_numberOfLostHits[Elecindex] = el->gsfTrack()->numberOfLostHits();
    ElecgsfTrack_numberOfLostHits[Elecindex] = el->gsfTrack()->trackerExpectedHitsInner().numberOfHits();//new (28 02 2012)
    double   AEff;
//    if(realdata){
      AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, el->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
//    }
//	else{
//      AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, el->superCluster()->eta(), ElectronEffectiveArea::kEleEAFall11MC);
//    }

    //cout<<rhoIso<<endl; 
    double rhoPrime = std::max(rhoIso, 0.0);

    chIso03[Elecindex] = el->chargedHadronIso();
    nhIso03[Elecindex] = el->neutralHadronIso();
    phIso03[Elecindex] = el->photonIso();
    puChIso03[Elecindex] = el->puChargedHadronIso();

    relIso[Elecindex] = ( chIso03[Elecindex] + nhIso03[Elecindex] + phIso03[Elecindex] ) / el->pt() ;
    relIsodb[Elecindex] = ( chIso03[Elecindex] + max(0.0, nhIso03[Elecindex] + phIso03[Elecindex] - 0.5*puChIso03[Elecindex]) )/ el->pt();
    relIsorho[Elecindex] = ( chIso03[Elecindex] + max(0.0, nhIso03[Elecindex] + phIso03[Elecindex] - rhoPrime*AEff) )/ el->pt();
//    cout<<relIso[Elecindex]<<"  "<<relIsorho[Elecindex]<<"   "<<rhoPrime<<"   "<<AEff<<endl;

    bool validKF= false; 
    reco::TrackRef myTrackRef =el->closestCtfTrackRef();
    validKF = (myTrackRef.isAvailable());
    validKF = (myTrackRef.isNonnull());  


    fMVAVar_fbrem[Elecindex]           =  el->fbrem();
    fMVAVar_kfchi2[Elecindex]          =  (validKF) ? myTrackRef->normalizedChi2() : 0 ;
    fMVAVar_kfhits[Elecindex]          =  (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; 
    fMVAVar_gsfchi2[Elecindex]         =  el->gsfTrack()->normalizedChi2();  
    fMVAVar_detacalo[Elecindex]        =  el->deltaEtaSeedClusterTrackAtCalo();
    fMVAVar_see[Elecindex]             =  el->sigmaIetaIeta();    //EleSigmaIEtaIEta
    //std::vector<float> vCov = EcalClusterLazyTools localCovariances(*(el->superCluster()->seed())) ;
    //if (!isnan(vCov[2])) fMVAVar_spp[Elecindex] = sqrt (vCov[2]);   //EleSigmaIPhiIPhi
    //else fMVAVar_spp[Elecindex] = 0.;    
    fMVAVar_etawidth[Elecindex]        =  el->superCluster()->etaWidth();
    fMVAVar_phiwidth[Elecindex]        =  el->superCluster()->phiWidth();
    fMVAVar_e1x5e5x5[Elecindex]        =  (el->e5x5()) !=0. ? 1.-(el->e1x5()/el->e5x5()) : -1. ;

    //fMVAVar_R9[Elecindex]              =  EcalClusterLazyTools::e3x3(*(el->superCluster()->seed())) / el->superCluster()->rawEnergy();
    fMVAVar_EoP[Elecindex]             =  el->eSuperClusterOverP();
    fMVAVar_IoEmIoP[Elecindex]         =  (1.0/el->ecalEnergy()) - (1.0 / el->p());  
    fMVAVar_eleEoPout[Elecindex]       =  el->eEleClusterOverPout();
    fMVAVar_PreShowerOverRaw[Elecindex]=  el->superCluster()->preshowerEnergy() / el->superCluster()->rawEnergy();
	++Elecindex; 
  }

  HFElecindex = 0;
  for(reco::RecoEcalCandidateCollection::const_iterator hfel = hfelec->begin(); hfel != hfelec->end();  ++hfel) {

	const reco::RecoEcalCandidate& HFcan = (*hfel);
	reco::SuperClusterRef theClusRef=HFcan.superCluster();
 	const reco::SuperCluster& hfECALSuperCluster=*theClusRef;
	const reco::HFEMClusterShapeRef clusShapeRef=(*electronHFClusterAssociation).find(theClusRef)->val;
	const reco::HFEMClusterShape& clusShape=*clusShapeRef;

	hfelec_Eta[HFElecindex] = hfel->eta();
	hfelec_Pt[HFElecindex] = hfel->pt();
	hfelec_Px[HFElecindex] = hfel->px();
	hfelec_Py[HFElecindex] = hfel->py();
	hfelec_Pz[HFElecindex] = hfel->pz();
	hfelec_Phi[HFElecindex] = hfel->phi();
	hfelec_E[HFElecindex] = hfel->energy();
	hfelec_M[HFElecindex] = hfel->mass();
	hfelec_charge[HFElecindex] = hfel->charge();
	hfelec_L9[HFElecindex] = clusShape.eLong3x3();
	hfelec_S9[HFElecindex] = clusShape.eShort3x3();
	hfelec_L25[HFElecindex] = clusShape.eLong5x5();
	hfelec_L1[HFElecindex] = clusShape.eLong1x1();
	hhfclustereta[HFElecindex] = hfECALSuperCluster.eta();
	hhfclusterphi[HFElecindex] = hfECALSuperCluster.phi();

    ++HFElecindex; 
  }

  
  //  int NJets = 50;

  for(int i=0; i<100; ++i){
    PFjet_mva[i]=-99;
    PFjetEta[i]=-99.;
    PFjetPhi[i]=-99.;
    PFjetPt[i]=-99.;
    PFjetCEMF[i]=-99.;
    PFjetNEMF[i]=-99.; 
    PFCorrjetPt[i]=-99;
    PFjetPx[i]=-99;
    PFjetPy[i]=-99;
    PFjetPz[i]=-99;
    PFjetE[i]=-99;
    PFjetM[i]=-99;
    PFHadEHF[i]=-99;
    PFEmEHF[i]=-99;
  }
  //----------------------------PF jets-------------------------------------------------------------------------------------
  for (int kkk = 0;kkk<100;kkk++){
    for (int jjj = 0;jjj<100;jjj++){
	  PFjetTrkVZ[kkk][jjj] = -99;
	  PFjetTrkPT[kkk][jjj] = -99;
    }
  }
  pfNjets=0;

  for ( unsigned int i=0; i<PFjets->size(); ++i ) {

    // for(edm::View<pat::Jet>::const_iterator i_jet = PFjets->begin(); i_jet != PFjets->end() && pfNjets < NJets; i_jet++){
    //    const math::XYZTLorentzVector theJet = i_jet->p4();

    //float mva   = (*puJetIdMVA)[jets->refAt(i)];
    ////      int    idflag = (*puJetIdFlag)[jets->refAt(i)];

    //float mva = i_jet-> puJetIdMVA();
 
	const pat::Jet & i_jet = PFjets->at(i);

    if(i_jet.pt()<15.)continue;

    PFjet_mva[pfNjets]   = (*puJetIDMVA)[PFjets->refAt(i)];

    //if(i_jet.eta()>2.9)cout<<PFjet_mva[pfNjets]<<endl;

    PFjetEta[pfNjets]=i_jet.eta();
    PFjetPhi[pfNjets]=i_jet.phi();
    PFjetPt[pfNjets]=i_jet.correctedJet(0).pt();
	PFjetPx[pfNjets]=i_jet.correctedJet(0).px();
    PFjetPy[pfNjets]=i_jet.correctedJet(0).py();
    PFjetPz[pfNjets]=i_jet.correctedJet(0).pz();
    PFjetE[pfNjets]=i_jet.correctedJet(0).energy();
    PFjetM[pfNjets]=i_jet.correctedJet(0).mass();
    PFCorrjetPt[pfNjets] = i_jet.pt();
    PFjetCEMF[pfNjets]=i_jet.chargedEmEnergyFraction();
    PFjetNEMF[pfNjets]=i_jet.neutralEmEnergyFraction();
    //cout<<pfNjets<<"   "<<PFCorrjetPt[pfNjets]/PFjetPt[pfNjets]<<endl;
    PFHadEHF[pfNjets]=i_jet.HFHadronEnergy();
    PFEmEHF[pfNjets]=i_jet.HFEMEnergy();    
    PFHadEHF[pfNjets]=i_jet.HFHadronEnergy();
    PFEmEHF[pfNjets]=i_jet.HFEMEnergy();    
    pfNtracks=0;

    const reco::TrackRefVector &tracks = i_jet.associatedTracks();

	for (reco::TrackRefVector::const_iterator iTrack = tracks.begin();iTrack != tracks.end(); ++iTrack) {
      PFjetTrkVZ[pfNjets][pfNtracks] = (**iTrack).vz();
	  PFjetTrkPT[pfNjets][pfNtracks] = (**iTrack).pt();
	  pfNtracks++;
    }
    pfNjets++;

  }

  //-----------------------------end of pf jets--------------------------------------------------------------------------
  
  myTree->Fill();//!!!!!!

}


// ------------ method called once each job just before starting event loop  ------------
//void ntupleGenerator::beginJob(const edm::EventSetup&)
void ntupleGenerator::beginJob()
{
  TFileDirectory TestDir = fs->mkdir("test");
  myTree = new TTree("Tree","Tree");
  myTree->Branch("event", &event, "event/I");
  myTree->Branch("run", &run, "run/I");
  myTree->Branch("lumi", &lumi, "lumi/I");
  myTree->Branch("bxnumber", &bxnumber, "bxnumber/I");
  myTree->Branch("realdata", &realdata, "realdata/I");
  myTree->Branch("Elecindex",&Elecindex,"Elecindex/I");
  myTree->Branch("hasMatchedConversion",hasMatchedConversion,"hasMatchedConversion[Elecindex]/B");
  myTree->Branch("Elecooeoop",Elecooeoop,"Elecooeoop[Elecindex]/F");
  myTree->Branch("ElecPt",ElecPt,"ElecPt[Elecindex]/F");
  myTree->Branch("ElecPx",ElecPx,"ElecPx[Elecindex]/F");
  myTree->Branch("ElecPy",ElecPy,"ElecPy[Elecindex]/F");
  myTree->Branch("ElecPz",ElecPz,"ElecPz[Elecindex]/F");
  myTree->Branch("ElecEta",ElecEta,"ElecEta[Elecindex]/F");
  myTree->Branch("ElecPhi",ElecPhi,"ElecPhi[Elecindex]/F");
  myTree->Branch("ElecE",ElecE,"ElecE[Elecindex]/F");  
  myTree->Branch("ElecM",ElecM,"ElecM[Elecindex]/F");
  myTree->Branch("ElecCharge",ElecCharge,"ElecCharge[Elecindex]/I");
  myTree->Branch("ElecGsfTrk_d0",ElecGsfTrk_d0,"ElecGsfTrk_d0[Elecindex]/F");
  myTree->Branch("ElecIsEB",ElecIsEB,"ElecIsEB[Elecindex]/I");
  myTree->Branch("ElecIsEE",ElecIsEE,"ElecIsEE[Elecindex]/I");
  myTree->Branch("Elecdr03TkSumPt",Elecdr03TkSumPt,"Elecdr03TkSumPt[Elecindex]/F");
  myTree->Branch("Elecdr03EcalRecHitSumEt",Elecdr03EcalRecHitSumEt,"Elecdr03EcalRecHitSumEt[Elecindex]/F");
  myTree->Branch("Elecdr03HcalTowerSumEt",Elecdr03HcalTowerSumEt,"Elecdr03HcalTowerSumEt[Elecindex]/F");
  myTree->Branch("ElecscSigmaIEtaIEta",ElecscSigmaIEtaIEta,"ElecscSigmaIEtaIEta[Elecindex]/F");
  myTree->Branch("ElecdeltaPhiSuperClusterTrackAtVtx",ElecdeltaPhiSuperClusterTrackAtVtx,"ElecdeltaPhiSuperClusterTrackAtVtx[Elecindex]/F");
  myTree->Branch("ElecdeltaEtaSuperClusterTrackAtVtx",ElecdeltaEtaSuperClusterTrackAtVtx,"ElecdeltaEtaSuperClusterTrackAtVtx[Elecindex]/F");
  myTree->Branch("ElechadronicOverEm",ElechadronicOverEm,"ElechadronicOverEm[Elecindex]/F");
  myTree->Branch("ElecgsfTrack_numberOfLostHits",ElecgsfTrack_numberOfLostHits,"ElecgsfTrack_numberOfLostHits[Elecindex]/F");
  myTree->Branch("ElecDcotTheta",ElecDcotTheta,"ElecDcotTheta[Elecindex]/F");
  myTree->Branch("ElecMVATrigId",ElecMVATrigId,"ElecMVATrigId[Elecindex]/F");
  myTree->Branch("ElecMVANonTrigId",ElecMVANonTrigId,"ElecMVANonTrigId[Elecindex]/F");
  myTree->Branch("Elecd0vtx",Elecd0vtx,"Elecd0vtx[Elecindex]/F");
  myTree->Branch("Elecdzvtx",Elecdzvtx,"Elecdzvtx[Elecindex]/F");
  myTree->Branch("relIso",relIso,"relIso[Elecindex]/F");
  myTree->Branch("relIsodb",relIsodb,"relIsodb[Elecindex]/F");
  myTree->Branch("relIsorho",relIsorho,"relIsorho[Elecindex]/F");
  myTree->Branch("fMVAVar_fbrem",fMVAVar_fbrem,"fMVAVar_fbrem[Elecindex]/F");
  myTree->Branch("fMVAVar_kfchi2",fMVAVar_kfchi2,"fMVAVar_kfchi2[Elecindex]/F");
  myTree->Branch("fMVAVar_kfhits",fMVAVar_kfhits,"fMVAVar_kfhits[Elecindex]/F");
  myTree->Branch("fMVAVar_gsfchi2",fMVAVar_gsfchi2,"fMVAVar_gsfchi2[Elecindex]/F");
  myTree->Branch("fMVAVar_detacalo",fMVAVar_detacalo,"fMVAVar_detacalo[Elecindex]/F");
  myTree->Branch("fMVAVar_see",fMVAVar_see,"fMVAVar_see[Elecindex]/F");
  myTree->Branch("fMVAVar_spp",fMVAVar_spp,"fMVAVar_spp[Elecindex]/F");
  myTree->Branch("fMVAVar_etawidth",fMVAVar_etawidth,"fMVAVar_etawidth[Elecindex]/F");
  myTree->Branch("fMVAVar_phiwidth",fMVAVar_phiwidth,"fMVAVar_phiwidth[Elecindex]/F");
  myTree->Branch("fMVAVar_e1x5e5x5",fMVAVar_e1x5e5x5,"fMVAVar_e1x5e5x5[Elecindex]/F");
  myTree->Branch("fMVAVar_R9",fMVAVar_R9,"fMVAVar_R9[Elecindex]/F");
  myTree->Branch("fMVAVar_EoP",fMVAVar_EoP,"fMVAVar_EoP[Elecindex]/F");
  myTree->Branch("fMVAVar_IoEmIoP",fMVAVar_IoEmIoP,"fMVAVar_IoEmIoP[Elecindex]/F");
  myTree->Branch("fMVAVar_eleEoPout",fMVAVar_eleEoPout,"fMVAVar_eleEoPout[Elecindex]/F");
  myTree->Branch("fMVAVar_PreShowerOverRaw",fMVAVar_PreShowerOverRaw,"fMVAVar_PreShowerOverRaw[Elecindex]/F");
  myTree->Branch("HFElecindex",&HFElecindex,"HFElecindex/I");
  myTree->Branch("hfelec_Pt",hfelec_Pt,"hfelec_Pt[HFElecindex]/F");
  myTree->Branch("hfelec_Eta",hfelec_Eta,"hfelec_Eta[HFElecindex]/F");
  myTree->Branch("hfelec_L9",hfelec_L9,"hfelec_L9[HFElecindex]/F");
  myTree->Branch("hfelec_S9",hfelec_S9,"hfelec_S9[HFElecindex]/F");
  myTree->Branch("hfelec_L25",hfelec_L25,"hfelec_L25[HFElecindex]/F");
  myTree->Branch("hfelec_L1",hfelec_L1,"hfelec_L1[HFElecindex]/F");
  myTree->Branch("hhfclustereta",hhfclustereta,"hhfclustereta[HFElecindex]/F");
  myTree->Branch("hhfclusterphi",hhfclusterphi,"hhfclusterphi[HFElecindex]/F");
  myTree->Branch("hfelec_Px",hfelec_Px,"hfelec_Px[HFElecindex]/F");
  myTree->Branch("hfelec_Py",hfelec_Py,"hfelec_Py[HFElecindex]/F");
  myTree->Branch("hfelec_Pz",hfelec_Pz,"hfelec_Pz[HFElecindex]/F");
  myTree->Branch("hfelec_Phi",hfelec_Phi,"hfelec_Phi[HFElecindex]/F");
  myTree->Branch("hfelec_E",hfelec_E,"hfelec_E[HFElecindex]/F");
  myTree->Branch("hfelec_M",hfelec_M,"hfelec_M[HFElecindex]/F");
  myTree->Branch("hfelec_charge",hfelec_charge,"hfelec_charge[HFElecindex]/F");
  myTree->Branch("techTrigger",techTrigger, "techTrigger[44]/I");
  //particles  
  myTree->Branch("par_index", &par_index, "par_index/I");
  myTree->Branch("ParticlePt", ParticlePt, "ParticlePt[par_index]/F");
  myTree->Branch("ParticleEta", ParticleEta, "ParticleEta[par_index]/F");
  myTree->Branch("ParticlePhi", ParticlePhi, "ParticlePhi[par_index]/F");
  myTree->Branch("ParticlePx", ParticlePx, "ParticlePx[par_index]/F");
  myTree->Branch("ParticlePy", ParticlePy, "ParticlePy[par_index]/F");
  myTree->Branch("ParticlePz", ParticlePz, "ParticlePz[par_index]/F");
  myTree->Branch("ParticleE", ParticleE, "ParticleE[par_index]/F");
  myTree->Branch("ParticleM", ParticleM, "ParticleM[par_index]/F");  
  myTree->Branch("ParticleId", ParticleId, "ParticleId[par_index]/I");
  myTree->Branch("mom", mom, "mom[par_index]/I");
  myTree->Branch("daug", daug, "daug[par_index]/I");
  myTree->Branch("ParticleStatus", ParticleStatus, "ParticleStatus[par_index]/I");
  myTree->Branch("ParticleMother", ParticleMother, "ParticleMother[par_index][5]/I");
  myTree->Branch("ParticleDaughter", ParticleDaughter, "ParticleDaughter[par_index][5]/I");
  myTree->Branch("pfMET", &pfMET, "pfMET/F");
  myTree->Branch("pfSET", &pfSET, "pfSET/F"); 
  myTree->Branch("pfMETX", &pfMETX, "pfMETX/F");
  myTree->Branch("pfMETY", &pfMETY, "pfMETY/F");
  myTree->Branch("muCorrMET", &muCorrMET, "muCorrMET/F");
  myTree->Branch("muCorrSET", &muCorrSET, "muCorrSET/F");	
  myTree->Branch("pfNjets",&pfNjets,"pfNjets/I");
  myTree->Branch("PFjet_mva",PFjet_mva,"PFjet_mva[pfNjets]/F");
  myTree->Branch("PFjetEta",PFjetEta,"PFjetEta[pfNjets]/F");
  myTree->Branch("PFjetPhi",PFjetPhi,"PFjetPhi[pfNjets]/F");
  myTree->Branch("PFjetPt",PFjetPt,"PFjetPt[pfNjets]/F");
  myTree->Branch("PFjetPx",PFjetPx,"PFjetPx[pfNjets]/F");
  myTree->Branch("PFjetPy",PFjetPy,"PFjetPy[pfNjets]/F"); 
  myTree->Branch("PFjetPz",PFjetPz,"PFjetPz[pfNjets]/F");
  myTree->Branch("PFjetE",PFjetE,"PFjetE[pfNjets]/F");
  myTree->Branch("PFjetM",PFjetM,"PFjetM[pfNjets]/F");
  myTree->Branch("PFCorrjetPt",PFCorrjetPt,"PFCorrjetPt[pfNjets]/F");
  myTree->Branch("PFjetCEMF",PFjetCEMF,"PFjetCEMF[pfNjets]/F");
  myTree->Branch("PFjetNEMF",PFjetNEMF,"PFjetNEMF[pfNjets]/F");
  myTree->Branch("PFHadEHF",PFHadEHF,"PFHadEHF[pfNjets]/F");
  myTree->Branch("PFEmEHF",PFEmEHF,"PFEmEHF[pfNjets]/F");
  myTree->Branch("pfNtracks",&pfNtracks,"pfNtracks/I");
  myTree->Branch("PFjetTrkVZ",PFjetTrkVZ,"PFjetTrkVZ[pfNjets][30]/F");
  myTree->Branch("PFjetTrkPT",PFjetTrkPT,"PFjetTrkPT[pfNjets][30]/F");
  myTree->Branch("nVertices",&nVertices,"nVertices/I");
  myTree->Branch("nGoodVertices",&nGoodVertices,"nGoodVertices/I");
  myTree->Branch("NumInteractions", &NumInteractions, "NumInteractions/I");
  myTree->Branch("TrueNumInteractions", &TrueNumInteractions, "TrueNumInteractions/F");
  myTree->Branch("Tnpv_check", &Tnpv_check, "Tnpv_check/F");
  myTree->Branch("BunchCrossing", &BunchCrossing,"BunchCrossing/I");
  myTree->Branch("MyWeight",&MyWeight,"MyWeight/D");
  myTree->Branch("vtxX",vtxX,"vtxX[nVertices]/F");
  myTree->Branch("vtxY",vtxY,"vtxY[nVertices]/F");
  myTree->Branch("vtxZ",vtxZ,"vtxZ[nVertices]/F");
  myTree->Branch("vtxXerr",vtxXerr,"vtxXerr[nVertices]/F");
  myTree->Branch("vtxYerr",vtxYerr,"vtxYerr[nVertices]/F");
  myTree->Branch("vtxZerr",vtxZerr,"vtxZerr[nVertices]/F");
  myTree->Branch("vtxisValid",vtxisValid,"vtxisValid[nVertices]/I");
  myTree->Branch("vtxisFake",vtxisFake,"vtxisFake[nVertices]/I");
  myTree->Branch("tr_1",&tr_1,"tr_1/I");
  myTree->Branch("tr_2",&tr_2,"tr_2/I");
  myTree->Branch("tr_3",&tr_3,"tr_3/I");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntupleGenerator::endJob() {
  myTree->Print();//!!!
  //  f->Write();
  //  f->Close();
}




float ntupleGenerator::DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}

float ntupleGenerator::DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}


//
//define this as a plug-in
DEFINE_FWK_MODULE(ntupleGenerator);

