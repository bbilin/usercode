#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "math.h"
#include <fstream>
#include <string>
#include <iostream>
#include <TStyle.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TF2.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TPostScript.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDBase.h"
#include "TLatex.h"

#include <TMath.h>
//#include <TFitter.h>
#include <TVector.h>
#include <TF1.h>
#include <TGraph.h>
//#include "TUnfold.h"
#include "TSystem.h"


using namespace std;
void ZtoEEanalyzer()
{


  gStyle->SetPalette(1);
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  float DeltaPhi(float phi1, float phi2);

  TChain myTree("demo/Tree");
int ismc =0, istau=0;
TString dataname;
//2012 DATA 

myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/2012A.root");
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/2012B_1.root");
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/2012B_2.root");
TFile *theFile = new TFile("rootfiles_new/02_10_data.root","RECREATE");
 dataname = "DATA";


//MC
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_tarball_1.root");
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_tarball_2.root");
 ismc =1;
TFile *theFile = new TFile("rootfiles_new/02_10_mc.root","RECREATE");
 dataname = "MC";
*/
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_10_50.root");
 ismc =1;
TFile *theFile = new TFile("rootfiles_new/02_10_mc_low.root","RECREATE");
 dataname = "MC_LOW";
*/

//TTbar
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_tt.root");
TFile *theFile = new TFile("rootfiles_new/02_10_bg_tt.root","RECREATE");
 dataname = "TTBAR";
*/

//WJets
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_WJets.root");
TFile *theFile = new TFile("rootfiles_new/02_10_bg_wjets.root","RECREATE");
 dataname = "WJets";
*/

//ZTauTau
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_tarball_1.root");
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_tarball_2.root");
 istau =1;
TFile *theFile = new TFile("rootfiles_new/02_10_bg_ztautau.root","RECREATE");
 dataname = "ZTauTau";
*/

//WW
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_ww.root");
TFile *theFile = new TFile("rootfiles_new/02_10_ww.root","RECREATE");
 dataname = "WW";
*/

//WZ
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_wz.root");
TFile *theFile = new TFile("rootfiles_new/02_10_wz.root","RECREATE");
 dataname = "WZ";
*/

//ZZ
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_zz.root");
TFile *theFile = new TFile("rootfiles_new/02_10_zz.root","RECREATE");
 dataname = "ZZ";
*/


TH1::AddDirectory(true);

  int event,run,lumi,bxnumber,realdata;
//  int hlt_trigger_fired;
  int Elecindex;
  float ElecPt[50], ElecEta[50], ElecPhi[50],ElecPx[50], ElecPy[50], ElecPz[50], ElecM[50], ElecE[50],  ElecscSigmaIEtaIEta[50], ElecdeltaPhiSuperClusterTrackAtVtx[50], ElecdeltaEtaSuperClusterTrackAtVtx[50], ElechadronicOverEm[50],ElecgsfTrack_numberOfLostHits[50],ElecMVATrigId[50],ElecMVANonTrigId[50],Elecd0vtx[50],Elecdzvtx[50],relIso[50],relIsodb[50],relIsorho[50],fMVAVar_fbrem[50],  fMVAVar_kfchi2[50], fMVAVar_kfhits[50],fMVAVar_gsfchi2[50],fMVAVar_detacalo[50],fMVAVar_see[50], fMVAVar_spp[50] , fMVAVar_etawidth[50],fMVAVar_phiwidth[50] ,fMVAVar_e1x5e5x5[50] , fMVAVar_R9[50] ,fMVAVar_EoP[50]  ,fMVAVar_IoEmIoP[50],fMVAVar_eleEoPout[50] ,fMVAVar_PreShowerOverRaw[50];

double MyWeight;

  //particle information
  int par_index, mom[50], daug[50];
 
  int ParticleId[50], ParticleStatus[50], ParticleMother[50][10], ParticleDaughter[50][10];
  //  int id_elec[20];
  int ElecCharge[50];
  //reco jets 
  int pfNjets;

  int techTrigger[50],nGoodVertices;

 bool hasMatchedConversion[50];

  float PFjetEta[50], PFjetPhi[50],PFCorrjetPt[50],PFHadEHF[50],PFEmEHF[50];
 //
  //met 


  int nVertices, tr_1, tr_2, tr_3;

  myTree.SetBranchAddress("event",  &event);
  myTree.SetBranchAddress("run", &run);
  myTree.SetBranchAddress("lumi", &lumi);

  myTree.SetBranchAddress("tr_1",&tr_1);
  myTree.SetBranchAddress("tr_2",&tr_2);
  myTree.SetBranchAddress("tr_3",&tr_3);
  myTree.SetBranchAddress("bxnumber", &bxnumber);
  myTree.SetBranchAddress("realdata",&realdata);
  myTree.SetBranchAddress("Elecindex",&Elecindex);
  myTree.SetBranchAddress("hasMatchedConversion",hasMatchedConversion);
  myTree.SetBranchAddress("ElecPt",ElecPt);
  myTree.SetBranchAddress("ElecEta",ElecEta);
  myTree.SetBranchAddress("ElecPhi",ElecPhi);
  myTree.SetBranchAddress("ElecPx",ElecPx);
  myTree.SetBranchAddress("ElecPy",ElecPy);
  myTree.SetBranchAddress("ElecPz",ElecPz);
  myTree.SetBranchAddress("ElecE",ElecE);
  myTree.SetBranchAddress("ElecM",ElecM);
  myTree.SetBranchAddress("ElecCharge",ElecCharge);
  myTree.SetBranchAddress("ElecMVATrigId",ElecMVATrigId);
  myTree.SetBranchAddress("ElecMVANonTrigId",ElecMVANonTrigId);
  myTree.SetBranchAddress("ElecscSigmaIEtaIEta", ElecscSigmaIEtaIEta);
  myTree.SetBranchAddress("ElecdeltaPhiSuperClusterTrackAtVtx", ElecdeltaPhiSuperClusterTrackAtVtx);
  myTree.SetBranchAddress("ElecdeltaEtaSuperClusterTrackAtVtx", ElecdeltaEtaSuperClusterTrackAtVtx);
  myTree.SetBranchAddress("ElechadronicOverEm", ElechadronicOverEm);
  myTree.SetBranchAddress("ElecgsfTrack_numberOfLostHits", ElecgsfTrack_numberOfLostHits);
myTree.SetBranchAddress("Elecd0vtx",Elecd0vtx);
myTree.SetBranchAddress("Elecdzvtx",Elecdzvtx);
//myTree.SetBranchAddress("chIso03",chIso03);
//myTree.SetBranchAddress("nhIso03",nhIso03);
//myTree.SetBranchAddress("phIso03",phIso03);
//myTree.SetBranchAddress("puChIso03",puChIso03);
myTree.SetBranchAddress("relIso",relIso);
myTree.SetBranchAddress("relIsodb",relIsodb);
myTree.SetBranchAddress("relIsorho",relIsorho);
myTree.SetBranchAddress("fMVAVar_fbrem",fMVAVar_fbrem);
myTree.SetBranchAddress("fMVAVar_kfchi2",fMVAVar_kfchi2);
myTree.SetBranchAddress("fMVAVar_kfhits",fMVAVar_kfhits);
myTree.SetBranchAddress("fMVAVar_gsfchi2",fMVAVar_gsfchi2);
myTree.SetBranchAddress("fMVAVar_detacalo",fMVAVar_detacalo);
myTree.SetBranchAddress("fMVAVar_see",fMVAVar_see);
myTree.SetBranchAddress("fMVAVar_spp",fMVAVar_spp);
myTree.SetBranchAddress("fMVAVar_etawidth",fMVAVar_etawidth);
myTree.SetBranchAddress("fMVAVar_phiwidth",fMVAVar_phiwidth);
myTree.SetBranchAddress("fMVAVar_e1x5e5x5",fMVAVar_e1x5e5x5);
myTree.SetBranchAddress("fMVAVar_R9",fMVAVar_R9);
myTree.SetBranchAddress("fMVAVar_EoP",fMVAVar_EoP);
myTree.SetBranchAddress("fMVAVar_IoEmIoP",fMVAVar_IoEmIoP);
myTree.SetBranchAddress("fMVAVar_eleEoPout",fMVAVar_eleEoPout);
myTree.SetBranchAddress("fMVAVar_PreShowerOverRaw",fMVAVar_PreShowerOverRaw);


  myTree.SetBranchAddress("techTrigger",techTrigger);

  
  myTree.SetBranchAddress("par_index", &par_index);
  myTree.SetBranchAddress("ParticleId", ParticleId);
  myTree.SetBranchAddress("mom", mom);
  myTree.SetBranchAddress("daug", daug);
  myTree.SetBranchAddress("ParticleStatus", ParticleStatus);
  myTree.SetBranchAddress("ParticleMother", ParticleMother);
  myTree.SetBranchAddress("ParticleDaughter", ParticleDaughter);

myTree.SetBranchAddress("nGoodVertices", &nGoodVertices);
 
myTree.SetBranchAddress("MyWeight", &MyWeight );

  myTree.SetBranchAddress("pfNjets",&pfNjets);
  myTree.SetBranchAddress("PFjetEta",PFjetEta);
  myTree.SetBranchAddress("PFjetPhi",PFjetPhi);
  myTree.SetBranchAddress("PFCorrjetPt",PFCorrjetPt);
  myTree.SetBranchAddress("PFHadEHF",PFHadEHF);
  myTree.SetBranchAddress("PFEmEHF",PFEmEHF);

 

  myTree.SetBranchAddress("nVertices",&nVertices);

  theFile->cd();
TH1F * h_MZe_l=  new TH1F ("h_MZe_l","h_MZe_l",80,0,200);

TH1F * h_MZe_l_mva0_do02_pfiso01_mhits0=  new TH1F ("h_MZe_l_mva0_do02_pfiso01_mhits0","h_MZe_l_mva0_do02_pfiso01_mhits0",80,0,200);
TH1F * h_eEta_mva0_do02_pfiso01_mhits0=  new TH1F ("h_eEta_mva0_do02_pfiso01_mhits0","h_eEta_mva0_do02_pfiso01_mhits0",120,-5.2,5.2);
TH1F * h_ePhi_mva0_do02_pfiso01_mhits0=  new TH1F ("h_ePhi_mva0_do02_pfiso01_mhits0","h_ePhi_mva0_do02_pfiso01_mhits0",120,-6.,6.);
TH1F * h_ePt_mva0_do02_pfiso01_mhits0=  new TH1F ("h_ePt_mva0_do02_pfiso01_mhits0","h_ePt_mva0_do02_pfiso01_mhits0",120,0.,500.);
TH1F * h_ePx_mva0_do02_pfiso01_mhits0=  new TH1F ("h_ePx_mva0_do02_pfiso01_mhits0","h_ePx_mva0_do02_pfiso01_mhits0",120,0.,500.);
TH1F * h_ePy_mva0_do02_pfiso01_mhits0=  new TH1F ("h_ePy_mva0_do02_pfiso01_mhits0","h_ePy_mva0_do02_pfiso01_mhits0",120,0.,500.);
TH1F * h_ePz_mva0_do02_pfiso01_mhits0=  new TH1F ("h_ePz_mva0_do02_pfiso01_mhits0","h_ePz_mva0_do02_pfiso01_mhits0",120,0.,500.);
TH1F * h_eM_mva0_do02_pfiso01_mhits0=  new TH1F ("h_eM_mva0_do02_pfiso01_mhits0","h_eM_mva0_do02_pfiso01_mhits0",120,0.,1000.);
TH1F * h_eE_mva0_do02_pfiso01_mhits0=  new TH1F ("h_eE_mva0_do02_pfiso01_mhits0","h_eE_mva0_do02_pfiso01_mhits0",120,0.,500.);
TH1F * h_eMVATrigId_mva0_do02_pfiso01_mhits0=  new TH1F ("h_eMVATrigId_mva0_do02_pfiso01_mhits0","h_eMVATrigId_mva0_do02_pfiso01_mhits0",120,-1.,1.);
TH1F * h_escSigmaIEtaIEta_mva0_do02_pfiso01_mhits0=  new TH1F ("h_escSigmaIEtaIEta_mva0_do02_pfiso01_mhits0","h_escSigmaIEtaIEta_mva0_do02_pfiso01_mhits0",120,0,0.1);
TH1F * h_edeltaPhiSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0=  new TH1F ("h_edeltaPhiSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0","h_edeltaPhiSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0",120,-8.,8.);
TH1F * h_edeltaEtaSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0 =  new TH1F ("h_edeltaEtaSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0","h_edeltaEtaSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0 ",120,-8.,8);
TH1F * h_ehadronicOverEm_mva0_do02_pfiso01_mhits0=  new TH1F ("h_ehadronicOverEm_mva0_do02_pfiso01_mhits0","h_ehadronicOverEm_mva0_do02_pfiso01_mhits0",120,0.,200.);
TH1F * h_egsfTrack_numberOfLostHits_mva0_do02_pfiso01_mhits0=  new TH1F ("h_egsfTrack_numberOfLostHits_mva0_do02_pfiso01_mhits0","h_egsfTrack_numberOfLostHits_mva0_do02_pfiso01_mhits0",8,-0.5,7.5);
TH1F * h_ed0vtx_mva0_do02_pfiso01_mhits0=  new TH1F ("h_ed0vtx_mva0_do02_pfiso01_mhits0","h_ed0vtx_mva0_do02_pfiso01_mhits0",120,-80,80);
TH1F * h_edzvtx_mva0_do02_pfiso01_mhits0=  new TH1F ("h_edzvtx_mva0_do02_pfiso01_mhits0","h_edzvtx_mva0_do02_pfiso01_mhits0",600,-400.,1000.);
//TH1F * h_echIso03=  new TH1F ("h_echIso03","h_echIso03",300,0.,1000.);
//TH1F * h_enhIso03=  new TH1F ("h_enhIso03","h_enhIso03",300,0,1000);
//TH1F * h_ephIso03=  new TH1F ("h_ephIso03","h_ephIso03",300,0,1000);
//TH1F * h_epuChIso03=  new TH1F ("h_epuChIso03","h_epuChIso03",300,0,1000);
TH1F * h_erelIso_mva0_do02_pfiso01_mhits0=  new TH1F ("h_erelIso_mva0_do02_pfiso01_mhits0","h_erelIso_mva0_do02_pfiso01_mhits0",300,0,1000);
TH1F * h_erelIsodb_mva0_do02_pfiso01_mhits0 =  new TH1F ("h_erelIsodb_mva0_do02_pfiso01_mhits0","h_erelIsodb_mva0_do02_pfiso01_mhits0",300,0,1000);
TH1F * h_erelIsorho_mva0_do02_pfiso01_mhits0=  new TH1F ("h_erelIsorho_mva0_do02_pfiso01_mhits0","h_erelIsorho_mva0_do02_pfiso01_mhits0_mva0_do02_pfiso01_mhits0",300,0,1000);
TH1F * h_efMVAVar_fbrem_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_fbrem_mva0_do02_pfiso01_mhits0","h_efMVAVar_fbrem_mva0_do02_pfiso01_mhits0",600,-1400.,0.);
TH1F * h_efMVAVar_kfchi2_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_kfchi2_mva0_do02_pfiso01_mhits0","h_efMVAVar_kfchi2_mva0_do02_pfiso01_mhits0",120,0,40.);
TH1F * h_efMVAVar_kfhits_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_kfhits_mva0_do02_pfiso01_mhits0","h_efMVAVar_kfhits_mva0_do02_pfiso01_mhits0",12,-2.5,9.5);
TH1F * h_efMVAVar_gsfchi2_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_gsfchi2_mva0_do02_pfiso01_mhits0","h_efMVAVar_gsfchi2_mva0_do02_pfiso01_mhits0",100000,0,900000.);
TH1F * h_efMVAVar_detacalo_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_detacalo_mva0_do02_pfiso01_mhits0","h_efMVAVar_detacalo_mva0_do02_pfiso01_mhits0",120,-8,8);
TH1F * h_efMVAVar_see_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_see_mva0_do02_pfiso01_mhits0","h_efMVAVar_see_mva0_do02_pfiso01_mhits0",480,-8,8);
TH1F * h_efMVAVar_spp_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_spp_mva0_do02_pfiso01_mhits0","h_efMVAVar_spp_mva0_do02_pfiso01_mhits0",120,-8,8);
TH1F * h_efMVAVar_etawidth_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_etawidth_mva0_do02_pfiso01_mhits0","h_efMVAVar_etawidth_mva0_do02_pfiso01_mhits0",120,0,12);
TH1F * h_efMVAVar_phiwidth_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_phiwidth_mva0_do02_pfiso01_mhits0","h_efMVAVar_phiwidth_mva0_do02_pfiso01_mhits0",120,0,12);
TH1F * h_efMVAVar_e1x5e5x5_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_e1x5e5x5_mva0_do02_pfiso01_mhits0","h_efMVAVar_e1x5e5x5_mva0_do02_pfiso01_mhits0",120,-50.,50.);
TH1F * h_efMVAVar_R9_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_R9_mva0_do02_pfiso01_mhits0","h_efMVAVar_R9_mva0_do02_pfiso01_mhits0",120,0,1000);
TH1F * h_efMVAVar_EoP_mva0_do02_pfiso01_mhits0 =  new TH1F ("h_efMVAVar_EoP_mva0_do02_pfiso01_mhits0","h_efMVAVar_EoP_mva0_do02_pfiso01_mhits0",800,0,2500);
TH1F * h_efMVAVar_IoEmIoP_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_IoEmIoP_mva0_do02_pfiso01_mhits0","h_efMVAVar_IoEmIoP_mva0_do02_pfiso01_mhits0",120,-4,4);
TH1F * h_efMVAVar_eleEoPout_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_eleEoPout_mva0_do02_pfiso01_mhits0","h_efMVAVar_eleEoPout_mva0_do02_pfiso01_mhits0",800,0,2000);
TH1F * h_efMVAVar_PreShowerOverRaw_mva0_do02_pfiso01_mhits0=  new TH1F ("h_efMVAVar_PreShowerOverRaw_mva0_do02_pfiso01_mhits0","h_efMVAVar_PreShowerOverRaw_mva0_do02_pfiso01_mhits0",120,0,20);

TH1F* h_nelec_fire_loose = new TH1F ("h_nelec_fire_loose","h_nelec_fire_loose",50,-0.5,49.5);
TH1F* h_no_good_vtx = new TH1F ("h_no_good_vtx","h_no_good_vtx",60,-0.5,59.5);
TH1F* h_no_good_vtx_noWeight = new TH1F ("h_no_good_vtx_noWeight","h_no_good_vtx_noWeight",60,-0.5,59.5);


TH1F*h_no_good_vtx_zevents[2][2]; TH1F*h_no_good_vtx_noWeight_zevents[2][2]; TH1F *h_deltaphi_z_jet[2][2];TH1F *h_mz[2][2]; TH1F *h_mz_same_sign[2][2];TH1F *h_mz_same_sign_anti_iso[2][2];TH1F *h_mz_same_sign_anti_iso_oppsign[2][2];TH1F *h_dielecphi[2][2];TH1F *h_dielec_PT[2][2];TH1F *h_dielec_rapidity[2][2];TH1F *h_mz_1j_all[2][2];TH1F *h_dielecphi_1j_all[2][2];TH1F *h_dielec_PT_1j_all[2][2];TH1F *h_dielec_rapidity_1j_all[2][2];TH1F *h_mz_1j_c[2][2];TH1F *h_dielecphi_1j_c[2][2];TH1F *h_dielec_PT_1j_c[2][2];TH1F *h_dielec_rapidity_1j_c[2][2];TH1F *h_mz_1j_hf[2][2];TH1F *h_dielecphi_1j_hf[2][2];TH1F *h_dielec_PT_1j_hf[2][2];TH1F *h_dielec_rapidity_1j_hf[2][2];TH1F *h_deltar_elec_PFjet[2][2];TH1F *h_numberofPFjets[2][2];TH1F *h_numberofPFjets_hbhe[2][2];TH1F *h_numberofPFjets_hf[2][2];TH1F *h_leading_jet_pt[2][2];TH1F *h_leading_jet_eta[2][2];TH1F *h_jet_pt_all[2][2];TH1F *h_jet_eta_all[2][2];TH1F *h_jet_pt_hbhe[2][2];TH1F *h_jet_eta_hbhe[2][2];TH1F *h_jet_pt_hf[2][2];TH1F *h_jet_eta_hf[2][2];
TH1F *h_mz_2j_all[2][2];TH1F *h_dielecphi_2j_all[2][2];TH1F *h_dielec_PT_2j_all[2][2];TH1F *h_dielec_rapidity_2j_all[2][2];TH1F *h_mz_2j_c[2][2];TH1F *h_dielecphi_2j_c[2][2];TH1F *h_dielec_PT_2j_c[2][2];TH1F *h_dielec_rapidity_2j_c[2][2];TH1F *h_mz_2j_hf[2][2];TH1F *h_dielecphi_2j_hf[2][2];TH1F *h_dielec_PT_2j_hf[2][2];TH1F *h_dielec_rapidity_2j_hf[2][2];

TH1F *h_mz_3j_all[2][2];TH1F *h_dielecphi_3j_all[2][2];TH1F *h_dielec_PT_3j_all[2][2];TH1F *h_dielec_rapidity_3j_all[2][2];TH1F *h_mz_3j_c[2][2];TH1F *h_dielecphi_3j_c[2][2];TH1F *h_dielec_PT_3j_c[2][2];TH1F *h_dielec_rapidity_3j_c[2][2];TH1F *h_mz_3j_hf[2][2];TH1F *h_dielecphi_3j_hf[2][2];TH1F *h_dielec_PT_3j_hf[2][2];TH1F *h_dielec_rapidity_3j_hf[2][2];

TH1F *h_mz_4j_all[2][2];TH1F *h_dielecphi_4j_all[2][2];TH1F *h_dielec_PT_4j_all[2][2];TH1F *h_dielec_rapidity_4j_all[2][2];TH1F *h_mz_4j_c[2][2];TH1F *h_dielecphi_4j_c[2][2];TH1F *h_dielec_PT_4j_c[2][2];TH1F *h_dielec_rapidity_4j_c[2][2];TH1F *h_mz_4j_hf[2][2];TH1F *h_dielecphi_4j_hf[2][2];TH1F *h_dielec_PT_4j_hf[2][2];TH1F *h_dielec_rapidity_4j_hf[2][2];

TH1F *h_mz_5j_all[2][2];TH1F *h_dielecphi_5j_all[2][2];TH1F *h_dielec_PT_5j_all[2][2];TH1F *h_dielec_rapidity_5j_all[2][2];TH1F *h_mz_5j_c[2][2];TH1F *h_dielecphi_5j_c[2][2];TH1F *h_dielec_PT_5j_c[2][2];TH1F *h_dielec_rapidity_5j_c[2][2];TH1F *h_mz_5j_hf[2][2];TH1F *h_dielecphi_5j_hf[2][2];TH1F *h_dielec_PT_5j_hf[2][2];TH1F *h_dielec_rapidity_5j_hf[2][2];

TH1F *h_jet_pt_all_2j[2][2];TH1F *h_jet_eta_all_2j[2][2];TH1F *h_jet_pt_hbhe_2j[2][2];TH1F *h_jet_eta_hbhe_2j[2][2];TH1F *h_jet_pt_hf_2j[2][2];TH1F *h_jet_eta_hf_2j[2][2];

TH1F *h_jet_pt_all_3j[2][2];TH1F *h_jet_eta_all_3j[2][2];TH1F *h_jet_pt_hbhe_3j[2][2];TH1F *h_jet_eta_hbhe_3j[2][2];TH1F *h_jet_pt_hf_3j[2][2];TH1F *h_jet_eta_hf_3j[2][2];

TH1F *h_jet_pt_all_4j[2][2];TH1F *h_jet_eta_all_4j[2][2];TH1F *h_jet_pt_hbhe_4j[2][2];TH1F *h_jet_eta_hbhe_4j[2][2];TH1F *h_jet_pt_hf_4j[2][2];TH1F *h_jet_eta_hf_4j[2][2];

TH1F *h_jet_pt_all_5j[2][2];TH1F *h_jet_eta_all_5j[2][2];TH1F *h_jet_pt_hbhe_5j[2][2];TH1F *h_jet_eta_hbhe_5j[2][2];TH1F *h_jet_pt_hf_5j[2][2];TH1F *h_jet_eta_hf_5j[2][2];TH1F *h_dielec_PT_same_sign[2][2];


//char name[100][100];
//char name2[100][100][100];
TString name2[200];

//sprintf(name2[0],"h_no_good_vtx_zevents_%i",i);
for(int i=0;i<2;i++){
for(int j=0;j<2;j++){
name2[0] = "h_no_good_vtx_zevents_";
name2[1] ="h_no_good_vtx_noWeight_zevents_";
name2[2] ="h_deltaphi_z_jet_";
name2[3] ="h_mz_same_sign_anti_iso_";
name2[4] ="h_mz_same_sign_anti_iso_oppsign_";
name2[5] ="h_mz_";
name2[6]="h_mz_same_sign_";
name2[7]="h_dielecphi_";
name2[8]="h_dielec_PT_";
name2[9]="h_dielec_rapidity_";
name2[10]="h_mz_1j_all_";
name2[11] ="h_dielecphi_1j_all_";
name2[12] ="h_dielec_PT_1j_all_";
name2[13] ="h_dielec_rapidity_1j_all_";
name2[14] ="h_mz_1j_c_";
name2[15] ="h_dielecphi_1j_c_";
name2[16] ="h_dielec_PT_1j_c_";
name2[17] ="h_dielec_rapidity_1j_c_";
name2[18] ="h_mz_1j_hf_";
name2[19] ="h_dielecphi_1j_hf_";
name2[20] ="h_dielec_PT_1j_hf_";
name2[21] ="h_dielec_rapidity_1j_hf_";
name2[22] ="h_deltar_elec_PFjet_";
name2[23] ="h_deltar_elec_PFjet_hf_";
name2[24] ="h_numberofPFjets_";
name2[25] ="h_numberofPFjets_hbhe_";
name2[26] ="h_numberofPFjets_hf_";
name2[27] ="h_leading_jet_pt_";
name2[28] ="h_leading_jet_eta_";
name2[29] ="h_jet_pt_all_";
name2[30] ="h_jet_eta_all_";
name2[31] ="h_jet_pt_hbhe_";
name2[32] ="h_jet_eta_hbhe_";
name2[33] ="h_jet_pt_hf_";
name2[34] ="h_jet_eta_hf_";
name2[35]="h_mz_2j_all_";
name2[36] ="h_dielecphi_2j_all_";
name2[37] ="h_dielec_PT_2j_all_";
name2[38] ="h_dielec_rapidity_2j_all_";
name2[39] ="h_mz_2j_c_";
name2[40] ="h_dielecphi_2j_c_";
name2[41] ="h_dielec_PT_2j_c_";
name2[42] ="h_dielec_rapidity_2j_c_";
name2[43] ="h_mz_2j_hf_";
name2[44] ="h_dielecphi_2j_hf_";
name2[45] ="h_dielec_PT_2j_hf_";
name2[46] ="h_dielec_rapidity_2j_hf_";
name2[47]="h_mz_3j_all_";
name2[48] ="h_dielecphi_3j_all_";
name2[49] ="h_dielec_PT_3j_all_";
name2[50] ="h_dielec_rapidity_3j_all_";
name2[51] ="h_mz_3j_c_";
name2[52] ="h_dielecphi_3j_c_";
name2[53] ="h_dielec_PT_3j_c_";
name2[54] ="h_dielec_rapidity_3j_c_";
name2[55] ="h_mz_3j_hf_";
name2[56] ="h_dielecphi_3j_hf_";
name2[57] ="h_dielec_PT_3j_hf_";
name2[58] ="h_dielec_rapidity_3j_hf_";
name2[59]="h_mz_4j_all_";
name2[60] ="h_dielecphi_4j_all_";
name2[61] ="h_dielec_PT_4j_all_";
name2[62] ="h_dielec_rapidity_4j_all_";
name2[63] ="h_mz_4j_c_";
name2[64] ="h_dielecphi_4j_c_";
name2[65] ="h_dielec_PT_4j_c_";
name2[66] ="h_dielec_rapidity_4j_c_";
name2[67] ="h_mz_4j_hf_";
name2[68] ="h_dielecphi_4j_hf_";
name2[69] ="h_dielec_PT_4j_hf_";
name2[70] ="h_dielec_rapidity_4j_hf_";
name2[71]="h_mz_5j_all_";
name2[72] ="h_dielecphi_5j_all_";
name2[73] ="h_dielec_PT_5j_all_";
name2[74] ="h_dielec_rapidity_5j_all_";
name2[75] ="h_mz_5j_c_";
name2[76] ="h_dielecphi_5j_c_";
name2[77] ="h_dielec_PT_5j_c_";
name2[78] ="h_dielec_rapidity_5j_c_";
name2[79] ="h_mz_5j_hf_";
name2[80] ="h_dielecphi_5j_hf_";
name2[81] ="h_dielec_PT_5j_hf_";
name2[82] ="h_dielec_rapidity_5j_hf_";
name2[83] ="h_jet_pt_all_2j_";
name2[84] ="h_jet_eta_all_2j_";
name2[85] ="h_jet_pt_hbhe_2j_";
name2[86] ="h_jet_eta_hbhe_2j_";
name2[87] ="h_jet_pt_hf_2j_";
name2[88] ="h_jet_eta_hf_2j_";
name2[89] ="h_jet_pt_all_3j_";
name2[90] ="h_jet_eta_all_3j_";
name2[91] ="h_jet_pt_hbhe_3j_";
name2[92] ="h_jet_eta_hbhe_3j_";
name2[93] ="h_jet_pt_hf_3j_";
name2[94] ="h_jet_eta_hf_3j_";
name2[95] ="h_jet_pt_all_4j_";
name2[96] ="h_jet_eta_all_4j_";
name2[97] ="h_jet_pt_hbhe_4j_";
name2[98] ="h_jet_eta_hbhe_4j_";
name2[99] ="h_jet_pt_hf_4j_";
name2[100] ="h_jet_eta_hf_4j_";
name2[101] ="h_jet_pt_all_5j_";
name2[102] ="h_jet_eta_all_5j_";
name2[103] ="h_jet_pt_hbhe_5j_";
name2[104] ="h_jet_eta_hbhe_5j_";
name2[105] ="h_jet_pt_hf_5j_";
name2[106] ="h_jet_eta_hf_5j_";
name2[107] ="h_dielec_PT_same_sign_";
for(int k=0;k<108;k++){
name2[k]+=(j);
name2[k]+="_";
name2[k]+=(i);
}
h_no_good_vtx_zevents[j][i]=new TH1F (name2[0],name2[0],60,-0.5,59.5);
h_no_good_vtx_noWeight_zevents[j][i]=new TH1F (name2[1],name2[1],60,-0.5,59.5);
h_deltaphi_z_jet[j][i]=  new TH1F (name2[2],name2[2],35,0,3.5);
h_mz_same_sign_anti_iso[j][i]=  new TH1F (name2[3],name2[3],70,0,1000);
h_mz_same_sign_anti_iso_oppsign[j][i]=  new TH1F (name2[4],name2[4],70,0,1000);
//h_elec_eta[j][i]=  new TH1F (name2[4],name2[4],80,-5,5);
h_mz[j][i]=  new TH1F (name2[5],name2[5],80,0,200);
h_mz_same_sign[j][i]=  new TH1F (name2[6],name2[6],80,0,200);
h_dielecphi[j][i]=  new TH1F ( name2[7], name2[7],20,-3,3);
h_dielec_PT[j][i]=  new TH1F ( name2[8], name2[8],100,0,200);
h_dielec_rapidity[j][i]=  new TH1F ( name2[9], name2[9],50,-5,5);
h_mz_1j_all[j][i]=  new TH1F ( name2[10], name2[10],80,0,200);
h_dielecphi_1j_all[j][i]=  new TH1F ( name2[11], name2[11],20,-3,3);
h_dielec_PT_1j_all[j][i]=  new TH1F ( name2[12], name2[12],100,0,200);
h_dielec_rapidity_1j_all[j][i]=  new TH1F ( name2[13], name2[13],50,-5,5);
h_mz_1j_c[j][i]=  new TH1F ( name2[14], name2[14],80,0,200);
h_dielecphi_1j_c[j][i]=  new TH1F ( name2[15], name2[15],20,-3,3);
h_dielec_PT_1j_c[j][i]=  new TH1F ( name2[16], name2[16],100,0,200);
h_dielec_rapidity_1j_c[j][i]=  new TH1F ( name2[17], name2[17],50,-5,5);
h_mz_1j_hf[j][i]=  new TH1F ( name2[18], name2[18],80,0,200);
h_dielecphi_1j_hf[j][i]=  new TH1F ( name2[19], name2[19],20,-3,3);
h_dielec_PT_1j_hf[j][i]=  new TH1F ( name2[20], name2[20],100,0,200);
h_dielec_rapidity_1j_hf[j][i]=  new TH1F ( name2[21], name2[21],50,-5,5);
h_deltar_elec_PFjet[j][i]=  new TH1F (name2[22], name2[22],50,0,20);
//h_deltar_elec_PFjet_hf[j][i]=  new TH1F (name2[23], name2[23],50,0,20);
h_numberofPFjets[j][i] = new TH1F (name2[24], name2[24],14,-0.5,13.5);
h_numberofPFjets_hbhe[j][i] = new TH1F (name2[25], name2[25],14,-0.5,13.5);
h_numberofPFjets_hf[j][i] = new TH1F (name2[26], name2[26],14,-0.5,13.5);
h_leading_jet_pt[j][i]=  new TH1F (name2[27], name2[27],80,0,200);
h_leading_jet_eta[j][i]=  new TH1F (name2[28], name2[28],40,-6,6);
h_jet_pt_all[j][i]=  new TH1F (name2[29], name2[29],80,0,200);
h_jet_eta_all[j][i]=  new TH1F (name2[30], name2[30],40,-6,6);
h_jet_pt_hbhe[j][i]=  new TH1F (name2[31], name2[31],80,0,200);
h_jet_eta_hbhe[j][i]=  new TH1F (name2[32], name2[32],40,-6,6);
h_jet_pt_hf[j][i]=  new TH1F (name2[33], name2[33],80,0,200);
h_jet_eta_hf[j][i]=  new TH1F (name2[34], name2[34],40,-6,6);
h_mz_2j_all[j][i]=  new TH1F ( name2[35], name2[35],80,0,200);
h_dielecphi_2j_all[j][i]=  new TH1F ( name2[36], name2[36],20,-3,3);
h_dielec_PT_2j_all[j][i]=  new TH1F ( name2[37], name2[37],100,0,200);
h_dielec_rapidity_2j_all[j][i]=  new TH1F ( name2[38], name2[38],50,-5,5);
h_mz_2j_c[j][i]=  new TH1F ( name2[39], name2[39],80,0,200);
h_dielecphi_2j_c[j][i]=  new TH1F ( name2[40], name2[40],20,-3,3);
h_dielec_PT_2j_c[j][i]=  new TH1F ( name2[41], name2[41],100,0,200);
h_dielec_rapidity_2j_c[j][i]=  new TH1F ( name2[42], name2[42],50,-5,5);
h_mz_2j_hf[j][i]=  new TH1F ( name2[43], name2[43],80,0,200);
h_dielecphi_2j_hf[j][i]=  new TH1F ( name2[44], name2[44],20,-3,3);
h_dielec_PT_2j_hf[j][i]=  new TH1F ( name2[45], name2[45],100,0,200);
h_dielec_rapidity_2j_hf[j][i]=  new TH1F ( name2[46], name2[46],50,-5,5);
h_mz_3j_all[j][i]=  new TH1F ( name2[47], name2[47],80,0,200);
h_dielecphi_3j_all[j][i]=  new TH1F ( name2[48], name2[48],20,-3,3);
h_dielec_PT_3j_all[j][i]=  new TH1F ( name2[49], name2[49],100,0,200);
h_dielec_rapidity_3j_all[j][i]=  new TH1F ( name2[50], name2[50],50,-5,5);
h_mz_3j_c[j][i]=  new TH1F ( name2[51], name2[51],80,0,200);
h_dielecphi_3j_c[j][i]=  new TH1F ( name2[52], name2[52],20,-3,3);
h_dielec_PT_3j_c[j][i]=  new TH1F ( name2[53], name2[53],100,0,200);
h_dielec_rapidity_3j_c[j][i]=  new TH1F ( name2[54], name2[54],50,-5,5);
h_mz_3j_hf[j][i]=  new TH1F ( name2[55], name2[55],80,0,200);
h_dielecphi_3j_hf[j][i]=  new TH1F ( name2[56], name2[56],20,-3,3);
h_dielec_PT_3j_hf[j][i]=  new TH1F ( name2[57], name2[57],100,0,200);
h_dielec_rapidity_3j_hf[j][i]=  new TH1F ( name2[58], name2[58],50,-5,5);
h_mz_4j_all[j][i]=  new TH1F ( name2[59], name2[59],80,0,200);
h_dielecphi_4j_all[j][i]=  new TH1F ( name2[60], name2[60],20,-3,3);
h_dielec_PT_4j_all[j][i]=  new TH1F ( name2[61], name2[61],100,0,200);
h_dielec_rapidity_4j_all[j][i]=  new TH1F ( name2[62], name2[62],50,-5,5);
h_mz_4j_c[j][i]=  new TH1F ( name2[63], name2[63],80,0,200);
h_dielecphi_4j_c[j][i]=  new TH1F ( name2[64], name2[64],20,-3,3);
h_dielec_PT_4j_c[j][i]=  new TH1F ( name2[65], name2[65],100,0,200);
h_dielec_rapidity_4j_c[j][i]=  new TH1F ( name2[66], name2[66],50,-5,5);
h_mz_4j_hf[j][i]=  new TH1F ( name2[67], name2[67],80,0,200);
h_dielecphi_4j_hf[j][i]=  new TH1F ( name2[68], name2[68],20,-3,3);
h_dielec_PT_4j_hf[j][i]=  new TH1F ( name2[69], name2[69],100,0,200);
h_dielec_rapidity_4j_hf[j][i]=  new TH1F ( name2[70], name2[70],50,-5,5);
h_mz_5j_all[j][i]=  new TH1F ( name2[71], name2[71],80,0,200);
h_dielecphi_5j_all[j][i]=  new TH1F ( name2[72], name2[72],20,-3,3);
h_dielec_PT_5j_all[j][i]=  new TH1F ( name2[73], name2[73],100,0,200);
h_dielec_rapidity_5j_all[j][i]=  new TH1F ( name2[74], name2[74],50,-5,5);
h_mz_5j_c[j][i]=  new TH1F ( name2[75], name2[75],80,0,200);
h_dielecphi_5j_c[j][i]=  new TH1F ( name2[76], name2[76],20,-3,3);
h_dielec_PT_5j_c[j][i]=  new TH1F ( name2[77], name2[77],100,0,200);
h_dielec_rapidity_5j_c[j][i]=  new TH1F ( name2[78], name2[78],50,-5,5);
h_mz_5j_hf[j][i]=  new TH1F ( name2[79], name2[79],80,0,200);
h_dielecphi_5j_hf[j][i]=  new TH1F ( name2[80], name2[80],20,-3,3);
h_dielec_PT_5j_hf[j][i]=  new TH1F ( name2[81], name2[81],100,0,200);
h_dielec_rapidity_5j_hf[j][i]=  new TH1F ( name2[82], name2[82],50,-5,5);
h_jet_pt_all_2j[j][i]=  new TH1F (name2[83], name2[83],80,0,200);
h_jet_eta_all_2j[j][i]=  new TH1F (name2[84], name2[84],40,-6,6);
h_jet_pt_hbhe_2j[j][i]=  new TH1F (name2[85], name2[85],80,0,200);
h_jet_eta_hbhe_2j[j][i]=  new TH1F (name2[86], name2[86],40,-6,6);
h_jet_pt_hf_2j[j][i]=  new TH1F (name2[87], name2[87],80,0,200);
h_jet_eta_hf_2j[j][i]=  new TH1F (name2[88], name2[88],40,-6,6);
h_jet_pt_all_3j[j][i]=  new TH1F (name2[89], name2[89],80,0,200);
h_jet_eta_all_3j[j][i]=  new TH1F (name2[90], name2[90],40,-6,6);
h_jet_pt_hbhe_3j[j][i]=  new TH1F (name2[91], name2[91],80,0,200);
h_jet_eta_hbhe_3j[j][i]=  new TH1F (name2[92], name2[92],40,-6,6);
h_jet_pt_hf_3j[j][i]=  new TH1F (name2[93], name2[93],80,0,200);
h_jet_eta_hf_3j[j][i]=  new TH1F (name2[94], name2[94],40,-6,6);
h_jet_pt_all_4j[j][i]=  new TH1F (name2[95], name2[95],80,0,200);
h_jet_eta_all_4j[j][i]=  new TH1F (name2[96], name2[96],40,-6,6);
h_jet_pt_hbhe_4j[j][i]=  new TH1F (name2[97], name2[97],80,0,200);
h_jet_eta_hbhe_4j[j][i]=  new TH1F (name2[98], name2[98],40,-6,6);
h_jet_pt_hf_4j[j][i]=  new TH1F (name2[99], name2[99],80,0,200);
h_jet_eta_hf_4j[j][i]=  new TH1F (name2[100], name2[100],40,-6,6);
h_jet_pt_all_5j[j][i]=  new TH1F (name2[101], name2[101],80,0,200);
h_jet_eta_all_5j[j][i]=  new TH1F (name2[102], name2[102],40,-6,6);
h_jet_pt_hbhe_5j[j][i]=  new TH1F (name2[103], name2[103],80,0,200);
h_jet_eta_hbhe_5j[j][i]=  new TH1F (name2[104], name2[104],40,-6,6);
h_jet_pt_hf_5j[j][i]=  new TH1F (name2[105], name2[105],80,0,200);
h_jet_eta_hf_5j[j][i]=  new TH1F (name2[106], name2[106],40,-6,6);
h_dielec_PT_same_sign[j][i]=  new TH1F ( name2[107], name2[107],100,0,200);
}
}
 TH1F *h_short=  new TH1F ( "h_short", "h_short",100,0,200);
 TH1F *h_long=  new TH1F ( "h_long", "h_long",100,0,200);







Int_t nevent = myTree.GetEntries();
// Int_t nevent =1000000;  

for (Int_t iev=0;iev<nevent;iev++) {
//if(iev!=21994 &&iev!=21983 )continue;
 	  if (iev%100000 == 0) cout<<dataname<<" ===>>> "<<iev<<"/"<<nevent<<endl;

    	myTree.GetEntry(iev); 

if(realdata && tr_3!=1) continue;

if (realdata) MyWeight=1; 

h_no_good_vtx-> Fill(nGoodVertices,MyWeight);
h_no_good_vtx_noWeight-> Fill(nGoodVertices);


int elec_ind_loose[50],elec_ind_anti_iso[50];

for(int i=0; i<50; i++){elec_ind_loose[i]=0;elec_ind_anti_iso[i]=0;}
bool MatchedConversion[50];
int eCharge[50]={},nelec_fire_loose=0, n_anti_iso=0;
float ePt[50] = {}, eEta[50] = {}, ePhi[50]= {}, ePx[50]={}, ePy[50]={}, ePz[50] = {}, eM[50]= {}, eE[50]={}, eMVATrigId[50]={}, escSigmaIEtaIEta[50]={}, edeltaPhiSuperClusterTrackAtVtx[50]={}, edeltaEtaSuperClusterTrackAtVtx[50]= {}, ehadronicOverEm[50]={} , egsfTrack_numberOfLostHits[50]={},ed0vtx[50]={},edzvtx[50]={}, /*echIso03[50]={},enhIso03[50]={},ephIso03[50]={},epuChIso03[50]={},*/erelIso[50]={},erelIsodb[50]={},erelIsorho[50]={},efMVAVar_fbrem[50]={},  efMVAVar_kfchi2[50]={}, efMVAVar_kfhits[50]={},efMVAVar_gsfchi2[50]={},efMVAVar_detacalo[50]={},efMVAVar_see[50]={}, efMVAVar_spp[50]={} , efMVAVar_etawidth[50]={},efMVAVar_phiwidth[50]={} ,efMVAVar_e1x5e5x5[50]={} , efMVAVar_R9[50]={} ,efMVAVar_EoP[50]={}  ,efMVAVar_IoEmIoP[50]={},efMVAVar_eleEoPout[50]={} ,efMVAVar_PreShowerOverRaw[50]={} ;  /*,eMVANonTrigId[50]={}*/; 
float jetEta[50]={}, jetPhi[50]={}, jetPt[50]={},jPt[50]={}, HadEHF[200]={}, EmEHF[200]={};
int sorted_jet_index[50]={};

for(int ii=0;ii<pfNjets;ii++){
jPt[ii]=0; jetPt[ii]=0; jetEta[ii]=0; jetPhi[ii]=0; jetPt[ii]=0,HadEHF[ii]=0, EmEHF[ii]=0; 
}



    int particle_id[50] = {};

    int nparticle_gen =0;


int mcIsElec=0;
int mcIsTau=0;

    for (int ppo = 0;ppo<50;++ppo) particle_id[ppo] = 9999; 
    for (int j = 0;j<par_index;++j){
      if (ParticleStatus[j] != 3) continue;
      particle_id[nparticle_gen] = ParticleId[j];
	if (particle_id[nparticle_gen]==11) mcIsElec=1;
	if (particle_id[nparticle_gen]==15) mcIsTau=1;

nparticle_gen++ ;//9 or 10 only
}


if(!realdata && ((ismc ==1 && mcIsElec==0)||(istau==1 && mcIsTau==0))) continue;

for(int ii=0;ii<pfNjets;ii++){
jPt[ii]=fabs(PFCorrjetPt[ii]);
}

int sort_jet = 0;
int ind=0;
TMath::Sort(pfNjets,jPt,sorted_jet_index);

for(int ii = 0 ; ii < pfNjets ; ii++){

 sort_jet = sorted_jet_index[ind];

jetEta[ind]=PFjetEta[sort_jet]  ;
jetPhi[ind]=PFjetPhi[sort_jet] ;
jetPt[ind]= jPt[sort_jet];
HadEHF[ind]=PFHadEHF[sort_jet] ;
EmEHF[ind]=PFEmEHF[sort_jet] ;
ind++;
}

for(int ii = 0 ; ii < Elecindex ; ii++){
MatchedConversion[ii]=hasMatchedConversion[ii];
eCharge[ii]=ElecCharge[ii];
eEta[ii] = ElecEta[ii];
ePhi[ii]= ElecPhi[ii];
ePt[ii] = ElecPt[ii];
ePx[ii]= ElecPx[ii];
ePy[ii]=  ElecPy[ii];
ePz[ii] =ElecPz[ii];
eM[ii]= ElecM[ii]; 
eE[ii]= ElecE[ii]; 
eMVATrigId[ii]=ElecMVATrigId[ii];
escSigmaIEtaIEta[ii]= ElecscSigmaIEtaIEta[ii]; 
edeltaPhiSuperClusterTrackAtVtx[ii]=ElecdeltaPhiSuperClusterTrackAtVtx[ii];
edeltaEtaSuperClusterTrackAtVtx[ii]= ElecdeltaEtaSuperClusterTrackAtVtx[ii]; 
ehadronicOverEm[ii]= ElechadronicOverEm[ii];
egsfTrack_numberOfLostHits[ii]= ElecgsfTrack_numberOfLostHits[ii];
ed0vtx[ii]=Elecd0vtx[ii];
edzvtx[ii]= Elecdzvtx[ii]; 
erelIso[ii]=relIso[ii];
erelIsodb[ii]=relIsodb[ii];
erelIsorho[ii]=relIsorho[ii];
efMVAVar_fbrem[ii]=  fMVAVar_fbrem[ii];  
efMVAVar_kfchi2[ii]= fMVAVar_kfchi2[ii]; 
efMVAVar_kfhits[ii]=fMVAVar_kfhits[ii];
efMVAVar_gsfchi2[ii]=fMVAVar_gsfchi2[ii];
efMVAVar_detacalo[ii]=fMVAVar_detacalo[ii];
efMVAVar_see[ii]= fMVAVar_see[ii] ;
efMVAVar_spp[ii]=  fMVAVar_spp[ii] ; 
efMVAVar_etawidth[ii]=fMVAVar_etawidth[ii];
efMVAVar_phiwidth[ii]= fMVAVar_phiwidth[ii]; 
efMVAVar_e1x5e5x5[ii]=  fMVAVar_e1x5e5x5[ii];  
efMVAVar_R9[ii]= fMVAVar_R9[ii] ;
efMVAVar_EoP[ii]=  fMVAVar_EoP[ii];  
efMVAVar_IoEmIoP[ii]=fMVAVar_IoEmIoP[ii];
efMVAVar_eleEoPout[ii]= fMVAVar_eleEoPout[ii]; 
efMVAVar_PreShowerOverRaw[ii]=fMVAVar_PreShowerOverRaw[ii];

if(eMVATrigId[ii]>0.&& ePt[ii] >30. 
&&(fabs(eEta[ii]) <1.4442 || fabs(eEta[ii])>1.566) 
&&fabs(eEta[ii])<2.5
&& fabs(ed0vtx[ii])<0.02 
&& egsfTrack_numberOfLostHits[ii]<=0 
&& !(MatchedConversion[ii])){
//
if(erelIsorho[ii]>0.1){
elec_ind_anti_iso[n_anti_iso]=ii;
n_anti_iso++;
			}
if(erelIsorho[ii]<0.1){
elec_ind_loose[nelec_fire_loose]=ii;
//cout<<ii<<"  "<<nelec_fire_loose<<"  "<<elec_ind_loose[nelec_fire_loose]<<endl;
nelec_fire_loose++;
			}
}
}
//cout<<"newwww"<<endl;
//if(nelec_fire_loose>2)cout<<"2den fazla electron var"<<endl;
h_nelec_fire_loose->Fill(nelec_fire_loose,MyWeight);
int a=-99;
int b=-99;
int aa=-99;
int bb=-99;
double    p1dotp2_l=-99.;
double MZe_l = -99.;
double    p1dotp2_ll=-99.;
double MZe_ll = -99.;
double dielecpt =-99.;
double dielecrapidity=-99.;
double dielecphi=-99.;
int selectt=0;
int select_ss=0;
if(n_anti_iso==2){
 aa= elec_ind_anti_iso[0], bb = elec_ind_anti_iso[1];
	       p1dotp2_ll = ePx[aa]*ePx[bb]+ePy[aa]*ePy[bb]+ePz[aa]*ePz[bb];
               MZe_ll = fabs(sqrt(pow(eM[aa],2)+pow(eM[bb],2)+2*(eE[aa]*eE[bb]-p1dotp2_ll)));
		if (eCharge[aa]*eCharge[bb] != -1)h_mz_same_sign_anti_iso[1][0]->Fill(MZe_ll,MyWeight);
		if (eCharge[aa]*eCharge[bb] == -1)h_mz_same_sign_anti_iso_oppsign[1][0]->Fill(MZe_ll,MyWeight);
}


if(nelec_fire_loose==2){
 a= elec_ind_loose[0], b = elec_ind_loose[1];
//           cout<<a<<"  "<<b<<endl;
	       p1dotp2_l = ePx[a]*ePx[b]+ePy[a]*ePy[b]+ePz[a]*ePz[b];
               MZe_l = fabs(sqrt(pow(eM[a],2)+pow(eM[b],2)+2*(eE[a]*eE[b]-p1dotp2_l)));


h_no_good_vtx_zevents[1][0]-> Fill(nGoodVertices,MyWeight); 
h_no_good_vtx_noWeight_zevents[1][0]-> Fill(nGoodVertices); 

                  dielecpt= sqrt(pow((ePx[a]+ePx[b]),2) +pow((ePy[a]+ePy[b]),2));
                  dielecrapidity = 0.5*log((eE[a]+eE[b]+ePz[a]+ePz[b])/(eE[a]+eE[b]-ePz[a]-ePz[b]));
                  dielecphi = atan2(ePy[a]+ePy[b],ePx[a]+ePx[b]);

if(eCharge[a]*eCharge[b]==-1)h_MZe_l->Fill(MZe_l,MyWeight);
if (eCharge[a]*eCharge[b] != -1){
      h_mz_same_sign[1][1]->Fill(MZe_l,MyWeight); 
      h_dielec_PT_same_sign[1][1]->Fill(dielecpt,MyWeight); 
      h_dielecphi[1][1]->Fill(dielecphi,MyWeight); 
      h_dielec_rapidity[1][1]->Fill(dielecrapidity,MyWeight); 

if(MZe_l<120 && MZe_l>60) select_ss=1;
}

if(eCharge[a]*eCharge[b]==-1 ){

h_no_good_vtx_zevents[1][0]-> Fill(nGoodVertices,MyWeight); 
h_no_good_vtx_noWeight_zevents[1][0]-> Fill(nGoodVertices); 

      h_mz[1][0]->Fill(MZe_l,MyWeight); 

      h_dielecphi[1][0]->Fill(dielecphi,MyWeight); 
      h_dielec_PT[1][0]->Fill(dielecpt,MyWeight); 
      h_dielec_rapidity[1][0]->Fill(dielecrapidity,MyWeight); 


h_MZe_l_mva0_do02_pfiso01_mhits0->Fill(MZe_l);
if(MZe_l<120. && MZe_l>60.){
selectt=1;
h_eEta_mva0_do02_pfiso01_mhits0->Fill( eEta[a],MyWeight);
h_ePhi_mva0_do02_pfiso01_mhits0->Fill( ePhi[a],MyWeight);
h_ePt_mva0_do02_pfiso01_mhits0->Fill( ePt[a],MyWeight);
h_ePx_mva0_do02_pfiso01_mhits0->Fill( ePx[a],MyWeight);
h_ePy_mva0_do02_pfiso01_mhits0->Fill(  ePy[a],MyWeight);
h_ePz_mva0_do02_pfiso01_mhits0->Fill(ePz[a],MyWeight);
h_eM_mva0_do02_pfiso01_mhits0->Fill( eM[a],MyWeight); 
h_eE_mva0_do02_pfiso01_mhits0->Fill( eE[a],MyWeight); 
h_eMVATrigId_mva0_do02_pfiso01_mhits0->Fill(eMVATrigId[a],MyWeight);
h_escSigmaIEtaIEta_mva0_do02_pfiso01_mhits0->Fill( escSigmaIEtaIEta[a],MyWeight); 
h_edeltaPhiSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0->Fill(edeltaPhiSuperClusterTrackAtVtx[a],MyWeight);
h_edeltaEtaSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0->Fill( edeltaEtaSuperClusterTrackAtVtx[a],MyWeight); 
h_ehadronicOverEm_mva0_do02_pfiso01_mhits0->Fill( ehadronicOverEm[a],MyWeight);
h_egsfTrack_numberOfLostHits_mva0_do02_pfiso01_mhits0->Fill( egsfTrack_numberOfLostHits[a],MyWeight);
h_ed0vtx_mva0_do02_pfiso01_mhits0->Fill(ed0vtx[a],MyWeight);
h_edzvtx_mva0_do02_pfiso01_mhits0->Fill( edzvtx[a],MyWeight); 
h_erelIso_mva0_do02_pfiso01_mhits0->Fill(erelIso[a],MyWeight);
h_erelIsodb_mva0_do02_pfiso01_mhits0->Fill(erelIsodb[a],MyWeight);
h_erelIsorho_mva0_do02_pfiso01_mhits0->Fill(erelIsorho[a],MyWeight);
h_efMVAVar_fbrem_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_fbrem[a],MyWeight);  
h_efMVAVar_kfchi2_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_kfchi2[a],MyWeight); 
h_efMVAVar_kfhits_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_kfhits[a],MyWeight);
h_efMVAVar_gsfchi2_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_gsfchi2[a],MyWeight);
h_efMVAVar_detacalo_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_detacalo[a],MyWeight);
h_efMVAVar_see_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_see[a],MyWeight) ;
h_efMVAVar_spp_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_spp[a],MyWeight) ; 
h_efMVAVar_etawidth_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_etawidth[a],MyWeight);
h_efMVAVar_phiwidth_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_phiwidth[a],MyWeight); 
h_efMVAVar_e1x5e5x5_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_e1x5e5x5[a],MyWeight);  
h_efMVAVar_R9_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_R9[a],MyWeight) ;
h_efMVAVar_EoP_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_EoP[a],MyWeight);  
h_efMVAVar_IoEmIoP_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_IoEmIoP[a],MyWeight);
h_efMVAVar_eleEoPout_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_eleEoPout[a],MyWeight); 
h_efMVAVar_PreShowerOverRaw_mva0_do02_pfiso01_mhits0->Fill(efMVAVar_PreShowerOverRaw[a],MyWeight);
}
}
}

if(select_ss==1){

      double deltaphi_ss = 999.0;
      int j_in_ss=0;
      int markPFjet_ss;

           int numberofPFjets_ss = 0;
           int jetindex_ss[30];
           int numberofPFjets_hbhe_ss = 0;
           int numberofPFjets_hf_ss = 0;

//cout<<pfNjets<<endl;
      for (int PF=0;PF < pfNjets ;PF++){
	  float deltar_ss = 99.;
	  float deltar2_ss = 99.;
	  markPFjet_ss = 0;
	   deltar_ss = fabs(DeltaR(jetEta[PF],eEta[a],jetPhi[PF],ePhi[a]));
	   deltar2_ss = fabs(DeltaR(jetEta[PF],eEta[b],jetPhi[PF],ePhi[b]));

	h_deltar_elec_PFjet[1][1] -> Fill(deltar_ss,MyWeight); 
	h_deltar_elec_PFjet[1][1] -> Fill(deltar2_ss,MyWeight); 

	  if (deltar_ss<0.5) markPFjet_ss = 1; 
	  if (deltar2_ss<0.5) markPFjet_ss = 1; 

 if (markPFjet_ss == 1)continue;


	if (jetPt[PF] < 30.) continue;


	  numberofPFjets_ss += 1;
	  jetindex_ss[j_in_ss]=PF;

//cout<<jetPt[PF]<<endl;

                  if (fabs(jetEta[PF])<2.4){
	  			numberofPFjets_hbhe_ss += 1;
         				     }

	                if(fabs(jetEta[PF])>2.4){
			        numberofPFjets_hf_ss += 1;

						  }
j_in_ss++;
}//jet for loop

h_numberofPFjets[1][1]->Fill(numberofPFjets_ss,MyWeight); 
h_numberofPFjets_hf[1][1]->Fill(numberofPFjets_hf_ss,MyWeight); 
h_numberofPFjets_hbhe[1][1]->Fill(numberofPFjets_hbhe_ss,MyWeight);

if(numberofPFjets_ss>0){

int jetindex1_ss= jetindex_ss[0];
//cout<<jetPt[jetindex1_ss]<<endl;
//cout<<numberofPFjets_ss<<endl;
//cout<<jetindex_ss[0]<<"   "<<jetindex1_ss<<endl;
deltaphi_ss = fabs(DeltaPhi(dielecphi,jetPhi[jetindex1_ss]));

if(numberofPFjets_ss==1)h_deltaphi_z_jet[1][1]->Fill(deltaphi_ss,MyWeight); 


	h_leading_jet_pt[1][1]->Fill(jetPt[jetindex1_ss],MyWeight); 
	h_leading_jet_eta[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);


if(numberofPFjets_ss==1){
//cout<<MyWeight<<endl;
   h_jet_pt_all[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_all[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
   h_mz_1j_all[1][1]->Fill(MZe_l,MyWeight); 
   h_dielecphi_1j_all[1][1]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_1j_all[1][1]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_1j_all[1][1]->Fill(dielecrapidity,MyWeight); 

if (fabs(jetEta[jetindex1_ss])<2.4){
if(numberofPFjets_hbhe_ss!=1) cout<<"DANGER THERE IS AN ERROR 1"<<endl;
  h_jet_pt_hbhe[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_hbhe[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_1j_c[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_1j_c[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_1j_c[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_1j_c[1][1]->Fill(dielecrapidity,MyWeight); 

}

if (fabs(jetEta[jetindex1_ss])>2.4){
if(numberofPFjets_hf_ss!=1) cout<<"DANGER THERE IS AN ERROR 2"<<endl;
  h_jet_pt_hf[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
  h_jet_eta_hf[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_1j_hf[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_1j_hf[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_1j_hf[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_1j_hf[1][1]->Fill(dielecrapidity,MyWeight); 
}

}

if(numberofPFjets_ss==2){
//cout<<MyWeight<<endl;
   h_jet_pt_all_2j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_all_2j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
   h_mz_2j_all[1][1]->Fill(MZe_l,MyWeight); 
   h_dielecphi_2j_all[1][1]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_2j_all[1][1]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_2j_all[1][1]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe_ss==2){
  h_jet_pt_hbhe_2j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_hbhe_2j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_2j_c[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_2j_c[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_2j_c[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_2j_c[1][1]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf_ss==2){
  h_jet_pt_hf_2j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
  h_jet_eta_hf_2j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_2j_hf[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_2j_hf[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_2j_hf[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_2j_hf[1][1]->Fill(dielecrapidity,MyWeight); 
}

}

if(numberofPFjets_ss==3){
//cout<<MyWeight<<endl;
   h_jet_pt_all_3j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_all_3j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
   h_mz_3j_all[1][1]->Fill(MZe_l,MyWeight); 
   h_dielecphi_3j_all[1][1]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_3j_all[1][1]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_3j_all[1][1]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe_ss==3){
  h_jet_pt_hbhe_3j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_hbhe_3j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_3j_c[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_3j_c[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_3j_c[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_3j_c[1][1]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf_ss==3){
  h_jet_pt_hf_3j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
  h_jet_eta_hf_3j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_3j_hf[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_3j_hf[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_3j_hf[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_3j_hf[1][1]->Fill(dielecrapidity,MyWeight); 
}

}

if(numberofPFjets_ss==4){
//cout<<MyWeight<<endl;
   h_jet_pt_all_4j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_all_4j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
   h_mz_4j_all[1][1]->Fill(MZe_l,MyWeight); 
   h_dielecphi_4j_all[1][1]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_4j_all[1][1]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_4j_all[1][1]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe_ss==4){
  h_jet_pt_hbhe_4j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_hbhe_4j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_4j_c[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_4j_c[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_4j_c[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_4j_c[1][1]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf_ss==4){
  h_jet_pt_hf_4j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
  h_jet_eta_hf_4j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_4j_hf[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_4j_hf[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_4j_hf[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_4j_hf[1][1]->Fill(dielecrapidity,MyWeight); 
}

}


if(numberofPFjets_ss==5){
//cout<<MyWeight<<endl;
   h_jet_pt_all_5j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_all_5j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
   h_mz_5j_all[1][1]->Fill(MZe_l,MyWeight); 
   h_dielecphi_5j_all[1][1]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_5j_all[1][1]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_5j_all[1][1]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe_ss==5){
  h_jet_pt_hbhe_5j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
   h_jet_eta_hbhe_5j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_5j_c[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_5j_c[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_5j_c[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_5j_c[1][1]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf_ss==5){
  h_jet_pt_hf_5j[1][1]->Fill(jetPt[jetindex1_ss],MyWeight);   
  h_jet_eta_hf_5j[1][1]->Fill(jetEta[jetindex1_ss],MyWeight);   
  h_mz_5j_hf[1][1]->Fill(MZe_l,MyWeight); 
  h_dielecphi_5j_hf[1][1]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_5j_hf[1][1]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_5j_hf[1][1]->Fill(dielecrapidity,MyWeight); 
}
}



}//if number of jets > 0

}//if samesign pair in mass window

//if(selectt==0&&(a==-99 || b==-99))cout<<a<<" "<<b<<"  "<<selectt<<endl; 
if(selectt==1){


      double deltaphi = 999.0;
      int j_in=0;
      int markPFjet;

           int numberofPFjets = 0;
           int jetindex[30];
	   double Long[50]={};
	   double Short[50]={};
           int numberofPFjets_hbhe = 0;
           int numberofPFjets_hf = 0;

//cout<<pfNjets<<endl;
      for (int PF=0;PF < pfNjets ;PF++){
	  float deltar = 99.;
	  float deltar2 = 99.;
	  markPFjet = 0;


	   deltar = fabs(DeltaR(jetEta[PF],eEta[a],jetPhi[PF],ePhi[a]));
	   deltar2 = fabs(DeltaR(jetEta[PF],eEta[b],jetPhi[PF],ePhi[b]));

	h_deltar_elec_PFjet[1][0] -> Fill(deltar,MyWeight); 
	h_deltar_elec_PFjet[1][0] -> Fill(deltar2,MyWeight); 

	  if (deltar<0.5) markPFjet = 1; 
	  if (deltar2<0.5) markPFjet = 1; 

 //if (deltar<0.5||deltar2<0.5)cout<<deltar<<"   "<<deltar2<<"  "<<PF<<"   "<<jetPt[PF]<<"  "<<markPFjet<<"  "<<j_in<<"  "<<iev<<endl;
 if (markPFjet == 1)continue;


	if (jetPt[PF] < 30.) continue;


	  numberofPFjets += 1;
	  jetindex[j_in]=PF;

//cout<<jetPt[PF]<<endl;

                  if (fabs(jetEta[PF])<2.4){
	  			numberofPFjets_hbhe += 1;
         				     }

	                if(fabs(jetEta[PF])>2.4){
			        numberofPFjets_hf += 1;
				Long[j_in]= 0.5 * HadEHF[PF] + EmEHF[PF];
				Short[j_in]= 0.5 * HadEHF[PF];
				h_short->Fill(Short[j_in],MyWeight);
				h_long->Fill(Long[j_in],MyWeight);
						  }
j_in++;
}

h_numberofPFjets[1][0]->Fill(numberofPFjets,MyWeight); 
h_numberofPFjets_hf[1][0]->Fill(numberofPFjets_hf,MyWeight); 
h_numberofPFjets_hbhe[1][0]->Fill(numberofPFjets_hbhe,MyWeight);

if(numberofPFjets>0){

int jetindex1= jetindex[0];
//cout<<jetPt[jetindex1]<<endl;
//cout<<numberofPFjets<<endl;
//cout<<jetindex[0]<<"   "<<jetindex1<<endl;
deltaphi = fabs(DeltaPhi(dielecphi,jetPhi[jetindex1]));

if(numberofPFjets==1)h_deltaphi_z_jet[1][0]->Fill(deltaphi,MyWeight); 


	h_leading_jet_pt[1][0]->Fill(jetPt[jetindex1],MyWeight); 
	h_leading_jet_eta[1][0]->Fill(jetEta[jetindex1],MyWeight);


if(numberofPFjets==1){
//cout<<MyWeight<<endl;
   h_jet_pt_all[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_all[1][0]->Fill(jetEta[jetindex1],MyWeight);   
   h_mz_1j_all[1][0]->Fill(MZe_l,MyWeight); 
   h_dielecphi_1j_all[1][0]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_1j_all[1][0]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_1j_all[1][0]->Fill(dielecrapidity,MyWeight); 

if (fabs(jetEta[jetindex1])<2.4){
if(numberofPFjets_hbhe!=1) cout<<"DANGER THERE IS AN ERROR 1"<<endl;
  h_jet_pt_hbhe[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_hbhe[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_1j_c[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_1j_c[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_1j_c[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_1j_c[1][0]->Fill(dielecrapidity,MyWeight); 

}

if (fabs(jetEta[jetindex1])>2.4){
if(numberofPFjets_hf!=1) cout<<"DANGER THERE IS AN ERROR 2"<<endl;
  h_jet_pt_hf[1][0]->Fill(jetPt[jetindex1],MyWeight);   
  h_jet_eta_hf[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_1j_hf[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_1j_hf[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_1j_hf[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_1j_hf[1][0]->Fill(dielecrapidity,MyWeight); 
}

}

if(numberofPFjets==2){
//cout<<MyWeight<<endl;
   h_jet_pt_all_2j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_all_2j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
   h_mz_2j_all[1][0]->Fill(MZe_l,MyWeight); 
   h_dielecphi_2j_all[1][0]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_2j_all[1][0]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_2j_all[1][0]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe==2){
  h_jet_pt_hbhe_2j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_hbhe_2j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_2j_c[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_2j_c[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_2j_c[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_2j_c[1][0]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf==2){
  h_jet_pt_hf_2j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
  h_jet_eta_hf_2j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_2j_hf[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_2j_hf[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_2j_hf[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_2j_hf[1][0]->Fill(dielecrapidity,MyWeight); 
}

}

if(numberofPFjets==3){
//cout<<MyWeight<<endl;
   h_jet_pt_all_3j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_all_3j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
   h_mz_3j_all[1][0]->Fill(MZe_l,MyWeight); 
   h_dielecphi_3j_all[1][0]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_3j_all[1][0]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_3j_all[1][0]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe==3){
  h_jet_pt_hbhe_3j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_hbhe_3j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_3j_c[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_3j_c[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_3j_c[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_3j_c[1][0]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf==3){
  h_jet_pt_hf_3j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
  h_jet_eta_hf_3j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_3j_hf[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_3j_hf[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_3j_hf[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_3j_hf[1][0]->Fill(dielecrapidity,MyWeight); 
}

}

if(numberofPFjets==4){
//cout<<MyWeight<<endl;
   h_jet_pt_all_4j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_all_4j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
   h_mz_4j_all[1][0]->Fill(MZe_l,MyWeight); 
   h_dielecphi_4j_all[1][0]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_4j_all[1][0]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_4j_all[1][0]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe==4){
  h_jet_pt_hbhe_4j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_hbhe_4j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_4j_c[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_4j_c[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_4j_c[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_4j_c[1][0]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf==4){
  h_jet_pt_hf_4j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
  h_jet_eta_hf_4j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_4j_hf[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_4j_hf[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_4j_hf[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_4j_hf[1][0]->Fill(dielecrapidity,MyWeight); 
}

}


if(numberofPFjets==5){
//cout<<MyWeight<<endl;
   h_jet_pt_all_5j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_all_5j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
   h_mz_5j_all[1][0]->Fill(MZe_l,MyWeight); 
   h_dielecphi_5j_all[1][0]->Fill(dielecphi,MyWeight); 
   h_dielec_PT_5j_all[1][0]->Fill(dielecpt,MyWeight); 
   h_dielec_rapidity_5j_all[1][0]->Fill(dielecrapidity,MyWeight); 

if (numberofPFjets_hbhe==5){
  h_jet_pt_hbhe_5j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
   h_jet_eta_hbhe_5j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_5j_c[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_5j_c[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_5j_c[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_5j_c[1][0]->Fill(dielecrapidity,MyWeight); 

}

if (numberofPFjets_hf==5){
  h_jet_pt_hf_5j[1][0]->Fill(jetPt[jetindex1],MyWeight);   
  h_jet_eta_hf_5j[1][0]->Fill(jetEta[jetindex1],MyWeight);   
  h_mz_5j_hf[1][0]->Fill(MZe_l,MyWeight); 
  h_dielecphi_5j_hf[1][0]->Fill(dielecphi,MyWeight); 
  h_dielec_PT_5j_hf[1][0]->Fill(dielecpt,MyWeight); 
  h_dielec_rapidity_5j_hf[1][0]->Fill(dielecrapidity,MyWeight); 
}

}

}//if number of jets > 0

}//if opposite sign elec pair in mass window

}//end of event loop


  theFile->Write();
  theFile->Close();


}//end void  

int main(int argc,char **argv){
  ZtoEEanalyzer();
}


float DeltaR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta2 - eta1;
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  float DELTAR = sqrt(pow(dphi,2)+pow(deta,2))*1.0;
  return DELTAR;
}
float DeltaPhi(float phi1, float phi2)
{
  float dphi = phi2 - phi1;
  if (fabs(dphi) > 3.14) dphi = 6.28 - fabs(dphi);
  return dphi;
}
