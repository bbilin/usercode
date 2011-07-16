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
// $Id: FwdCalib.cc,v 1.1 2011/06/09 08:31:56 bbilin Exp $
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

//elec

#include "DataFormats/PatCandidates/interface/Electron.h"

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
  string CaloJetAlg;
 
  edm::Service<TFileService> fs;
  TTree * myTree;
  int event, run,lumi,bxnumber,realdata;
  int hlt_trigger_fired;

float ElecPt[50],ElecE[50],ElecM[50],ElecPx[50], ElecPy[50],ElecPz[50],ElecEta[50],ElecPhi[50],ElecGsfTrk_d0[50],Elecdr03TkSumPt[50],Elecdr03EcalRecHitSumEt[50],Elecdr03HcalTowerSumEt[50],ElecscSigmaIEtaIEta[50],ElecdeltaPhiSuperClusterTrackAtVtx[50],ElecdeltaEtaSuperClusterTrackAtVtx[50],ElechadronicOverEm[50],ElecgsfTrack_numberOfLostHits[50];
  int Elecindex,ElecCharge[50],ElecIsEB[50], ElecIsEE[50]; 
 

  //particle information
  int par_index, mom[500], daug[500];
  float ParticlePt[500], ParticleEta[500], ParticlePhi[500], ParticlePx[500], ParticlePy[500], ParticlePz[500], ParticleE[500], ParticleM[500];
  int ParticleId[500], ParticleStatus[500], ParticleMother[500][100], ParticleDaughter[500][100];



  //reco jets 
  int pfNjets,pfNtracks;
  
  int CaloNjets;

  float PFjetEta[100], PFjetPhi[100],PFjetPt[100],PFCorrjetPt[100],PFjetCEMF[100],PFjetNEMF[100],PFjetPx[100],PFjetPy[100],PFjetPz[100],PFjetE[100],PFjetM[100],PFHadEHF[100],PFEmEHF[100];

 float CalojetEta[100], CalojetPhi[100],CalojetPt[100],CaloCorrjetPt[100],CalojetCEMF[100],CalojetNEMF[100],CalojetPx[100],CalojetPy[100],CalojetPz[100],CalojetE[100],CalojetM[100],CaloHadEHF[100],CaloEmEHF[100];


  float PFjetTrkVZ[100][100],PFjetTrkPT[100][100];
  float vtxZ[100],vtxZerr[100];
  float vtxY[100],vtxYerr[100];
  float vtxX[100],vtxXerr[100];
  int vtxisValid[100],vtxisFake[100];
  int nVertices,nGoodVertices; 
  int techTrigger[44];
  //met 
  float caloMET, caloSET, pfMET, pfSET;
  float caloMETX, pfMETX;
  float caloMETY, pfMETY;
  float muCorrMET, muCorrSET;



//  reco::helper::JetIDHelper *jetID;

};



ntupleGenerator::ntupleGenerator(const edm::ParameterSet& iConfig)
{   

  elecTag_ = iConfig.getParameter<edm::InputTag>("elecTag");
  PfJetAlg = iConfig.getParameter<string>("PfJetAlg");
  CaloJetAlg = iConfig.getParameter<string>("CaloJetAlg");

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




 
 
  Handle< edm::View<pat::MET> > caloMEThandle;
  iEvent.getByLabel("patMETs", caloMEThandle); 

  Handle< edm::View<pat::MET> > pfMEThandle;
  iEvent.getByLabel("patMETsPF", pfMEThandle); 


  Handle<edm::View<pat::Jet> > Calojets;
  iEvent.getByLabel(CaloJetAlg,Calojets);
//  edm::View<pat::Jet>::const_iterator i_jet;
//const edm::View<pat::Jet> & i_jet = *Calojets;
 

  Handle<edm::View<pat::Jet> > PFjets;
  iEvent.getByLabel(PfJetAlg,PFjets);
//const edm::View<pat::Jet> & calo_jet = *PFjets;

 Handle<edm::View<pat::Electron> > elecs_h;	
  iEvent.getByLabel(elecTag_,elecs_h);
  const edm::View<pat::Electron>* elec = elecs_h.failedToGet() ? 0 : &*elecs_h;



  Handle<BeamSpot> beamSpotHandle;
  if (!iEvent.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle)) {
    cout<< ">>> No beam spot found !!!" <<endl;
    return;
  }
  
  Handle<VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);  

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

  




  event = iEvent.id().event();
  run = iEvent.id().run();
  lumi = iEvent.luminosityBlock();
  bxnumber = iEvent.bunchCrossing();
  realdata = iEvent.isRealData();


 
  //----------------------------------------------------------------------------


  //--------------------met-------------------------------------------------
  caloMET = (caloMEThandle->front()).et();
  caloSET = (caloMEThandle->front()).sumEt();
  pfMET = (pfMEThandle->front()).et();
  pfSET = (pfMEThandle->front()).sumEt(); 

  caloMETX = (caloMEThandle->front()).px();
  caloMETY = (caloMEThandle->front()).py();

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

  }
  Elecindex = 0;
  for(edm::View<pat::Electron>::const_iterator el = elec->begin(); el != elec->end();  ++el) {
	ElecPt[Elecindex] = el->pt();
        ElecPx[Elecindex] = el->px();
	ElecPy[Elecindex] = el->py();
	ElecPz[Elecindex] = el->pz();
	ElecEta[Elecindex] = el->eta();
	ElecPhi[Elecindex] = el->phi();	
        ElecE[Elecindex] = el->energy();
	ElecM[Elecindex] = el->mass();
	ElecCharge[Elecindex] = el->charge();
	ElecGsfTrk_d0[Elecindex] = el->gsfTrack()->d0();
	ElecIsEB[Elecindex] = el->isEB();
	ElecIsEE[Elecindex] = el->isEE();
	Elecdr03TkSumPt[Elecindex] = el->dr03TkSumPt();
	Elecdr03EcalRecHitSumEt[Elecindex] = el->dr03EcalRecHitSumEt();
	Elecdr03HcalTowerSumEt[Elecindex] = el->dr03HcalTowerSumEt();
	ElecscSigmaIEtaIEta[Elecindex] = el->scSigmaIEtaIEta();
	ElecdeltaPhiSuperClusterTrackAtVtx[Elecindex] = el->deltaPhiSuperClusterTrackAtVtx();
	ElecdeltaEtaSuperClusterTrackAtVtx[Elecindex] = el->deltaEtaSuperClusterTrackAtVtx();
	ElechadronicOverEm[Elecindex] = el->hadronicOverEm();
	ElecgsfTrack_numberOfLostHits[Elecindex] = el->gsfTrack()->numberOfLostHits();
        ++Elecindex; 
  }

  
  int NJets = 20;

for(int i=0; i<100; ++i){
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

    CalojetEta[i]=-99.;
    CalojetPhi[i]=-99.;
    CalojetPt[i]=-99.;
    CalojetCEMF[i]=-99.;
    CalojetNEMF[i]=-99.; 
    CaloCorrjetPt[i]=-99;
    CalojetPx[i]=-99;
    CalojetPy[i]=-99;
    CalojetPz[i]=-99;
    CalojetE[i]=-99;
    CalojetM[i]=-99;
    CaloHadEHF[i]=-99;
    CaloEmEHF[i]=-99;
}
  //----------------------------PF jets-------------------------------------------------------------------------------------
  for (int kkk = 0;kkk<100;kkk++){
    for (int jjj = 0;jjj<100;jjj++){
	PFjetTrkVZ[kkk][jjj] = -99;
	PFjetTrkPT[kkk][jjj] = -99;
    }
  }
  pfNjets=0;
  CaloNjets=0;
  for(edm::View<pat::Jet>::const_iterator i_jet = PFjets->begin(); i_jet != PFjets->end(); i_jet++){
    const math::XYZTLorentzVector theJet = i_jet->p4();
    PFjetEta[pfNjets]=i_jet->eta();
    PFjetPhi[pfNjets]=i_jet->phi();
    PFjetPt[pfNjets]=i_jet->correctedJet(0).pt();
     PFjetPx[pfNjets]=i_jet->correctedJet(0).px();
   PFjetPy[pfNjets]=i_jet->correctedJet(0).py();
   PFjetPz[pfNjets]=i_jet->correctedJet(0).pz();
   PFjetE[pfNjets]=i_jet->correctedJet(0).energy();
   PFjetM[pfNjets]=i_jet->correctedJet(0).mass();
    PFCorrjetPt[pfNjets] = i_jet->pt();
    PFjetCEMF[pfNjets]=i_jet->chargedEmEnergyFraction();
    PFjetNEMF[pfNjets]=i_jet->neutralEmEnergyFraction();

    PFHadEHF[pfNjets]=i_jet->HFHadronEnergy();
    PFEmEHF[pfNjets]=i_jet->HFEMEnergy();    

          cout<<PFjetEta[pfNjets]<<"  "<<PFjetPt[pfNjets]<<"  "<<PFCorrjetPt[pfNjets]<<endl;
    // ---- This accesses to the vertex Z of the track-base constituents of pfjets
    pfNtracks=0;



    const reco::TrackRefVector &tracks = i_jet->associatedTracks();
    for (reco::TrackRefVector::const_iterator iTrack = tracks.begin();iTrack != tracks.end(); ++iTrack) {
          PFjetTrkVZ[pfNjets][pfNtracks] = (**iTrack).vz();
          PFjetTrkPT[pfNjets][pfNtracks] = (**iTrack).pt();
          pfNtracks++;
    }
    pfNjets++;
//
  }

  //-----------------------------end of pf jets--------------------------------------------------------------------------

 for(edm::View<pat::Jet>::const_iterator i_jet = Calojets->begin(); i_jet != Calojets->end(); i_jet++){
    const math::XYZTLorentzVector theJet = i_jet->p4();
   CalojetEta[CaloNjets]=i_jet->eta();
    CalojetPhi[CaloNjets]=i_jet->phi();
    CalojetPt[CaloNjets]=i_jet->correctedJet(0).pt();
     CalojetPx[CaloNjets]=i_jet->correctedJet(0).px();
   CalojetPy[CaloNjets]=i_jet->correctedJet(0).py();
   CalojetPz[CaloNjets]=i_jet->correctedJet(0).pz();
   CalojetE[CaloNjets]=i_jet->correctedJet(0).energy();
   CalojetM[CaloNjets]=i_jet->correctedJet(0).mass();
    CaloHadEHF[CaloNjets]=i_jet->hadEnergyInHF();
    CaloEmEHF[CaloNjets]=i_jet->emEnergyInHF();    
 //   CaloCorrjetPt[CaloNjets] = i_jet->pt();
//    CalojetCEMF[CaloNjets]=i_jet->chargedEmEnergyFraction();
//    CalojetNEMF[CaloNjets]=i_jet->neutralEmEnergyFraction();
 //           cout<<CalojetEta[CaloNjets]<<"  "<<CalojetPt[CaloNjets]<<"  "<<CaloCorrjetPt[CaloNjets]<<endl;
    // ---- This accesses to the vertex Z of the track-base constituents of Calojets

    
    CaloNjets++;
//
  }



  
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
  

  myTree->Branch("caloMET", &caloMET, "caloMET/F");
  myTree->Branch("caloSET", &caloSET, "caloSET/F");
  myTree->Branch("pfMET", &pfMET, "pfMET/F");
  myTree->Branch("pfSET", &pfSET, "pfSET/F"); 
  myTree->Branch("caloMETX", &caloMETX, "caloMETX/F");
  myTree->Branch("pfMETX", &pfMETX, "pfMETX/F");
  myTree->Branch("caloMETY", &caloMETY, "caloMETY/F");
  myTree->Branch("pfMETY", &pfMETY, "pfMETY/F");
  myTree->Branch("muCorrMET", &muCorrMET, "muCorrMET/F");
  myTree->Branch("muCorrSET", &muCorrSET, "muCorrSET/F");	


  myTree->Branch("CaloNjets",&CaloNjets,"CaloNjets/I");
  myTree->Branch("CalojetEta",CalojetEta,"CalojetEta[CaloNjets]/F");
  myTree->Branch("CalojetPhi",CalojetPhi,"CalojetPhi[CaloNjets]/F");
  myTree->Branch("CalojetPt",CalojetPt,"CalojetPt[CaloNjets]/F");
  myTree->Branch("CalojetPx",CalojetPx,"CalojetPx[CaloNjets]/F");
  myTree->Branch("CalojetPy",CalojetPy,"CalojetPy[CaloNjets]/F"); 
  myTree->Branch("CalojetPz",CalojetPz,"CalojetPz[CaloNjets]/F");
  myTree->Branch("CalojetE",CalojetE,"CalojetE[CaloNjets]/F");
  myTree->Branch("CalojetM",CalojetM,"CalojetM[CaloNjets]/F");
  myTree->Branch("CaloHadEHF",CaloHadEHF,"CaloHadEHF[CaloNjets]/F");
  myTree->Branch("CaloEmEHF",CaloEmEHF,"CaloEmEHF[CaloNjets]/F");

//  myTree->Branch("CaloCorrjetPt",CaloCorrjetPt,"CaloCorrjetPt[CaloNjets]/F");
//  myTree->Branch("CalojetCEMF",CalojetCEMF,"CalojetCEMF[CaloNjets]/F");
//  myTree->Branch("CalojetNEMF",CalojetNEMF,"CalojetNEMF[CaloNjets]/F");

  myTree->Branch("pfNjets",&pfNjets,"pfNjets/I");
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
  myTree->Branch("PFHadEHF",PFHadEHF,"PFHadEHF[CaloNjets]/F");
  myTree->Branch("PFEmEHF",PFEmEHF,"PFEmEHF[CaloNjets]/F");
  myTree->Branch("pfNtracks",&pfNtracks,"pfNtracks/I");
  myTree->Branch("PFjetTrkVZ",PFjetTrkVZ,"PFjetTrkVZ[pfNjets][30]/F");
  myTree->Branch("PFjetTrkPT",PFjetTrkPT,"PFjetTrkPT[pfNjets][30]/F");
	

  myTree->Branch("nVertices",&nVertices,"nVertices/I");
  myTree->Branch("nGoodVertices",&nGoodVertices,"nGoodVertices/I");
  myTree->Branch("vtxX",vtxX,"vtxX[nVertices]/F");
  myTree->Branch("vtxY",vtxY,"vtxY[nVertices]/F");
  myTree->Branch("vtxZ",vtxZ,"vtxZ[nVertices]/F");
  myTree->Branch("vtxXerr",vtxXerr,"vtxXerr[nVertices]/F");
  myTree->Branch("vtxYerr",vtxYerr,"vtxYerr[nVertices]/F");
  myTree->Branch("vtxZerr",vtxZerr,"vtxZerr[nVertices]/F");
  myTree->Branch("vtxisValid",vtxisValid,"vtxisValid[nVertices]/I");
  myTree->Branch("vtxisFake",vtxisFake,"vtxisFake[nVertices]/I");

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
