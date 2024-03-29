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
/*double c_may10[10] ={1.0117,0.994,1.006,0.9976,0.9967,0.9963,0.9963,1.0028,0.9950,1.0126};//elec correction factors//
double c_promtv4[10] ={1.0474,1.0018,1.0100,1.0007,1.0032,1.0004,0.9990,1.0109,1.0036,1.0597}; 
double c_aug5[10] ={0.9937,0.9951,1.0142,1.0047,1.0014,1.0008,1.0023,1.0182,0.9880,0.9884};
double c_promtv6[10] ={0.9713,0.9786,1.0181,1.0064,1.0041,1.0025,1.0046,1.0167,0.9752,0.9659};
*/
/*
double par0[10]={394.199,353.677,414.954,430.286,397.128,424.611,500,500,500,500,};
double par1[10]={-3.9303,-31.3253,-7.96372,2.11221,-6.87922,1.66457,20,20,20,20,};
double par2[10]={35.5977,37.9435,35.1436,34.1395,35.6599,34.4009,30,30,30,30,};
double par3[10]={38.024,20.9857,29.5735,37.962,34.7899,38.268,50,50,50,50,};
*/

double par0[10]={409.489,350.422,419.937,428.576,399.475,435.435,500,500,500,500,};
double par1[10]={5.73172,-32.7205,-4.86986,2.66733,-6.01504,7.09596,20,20,20,20,};
double par2[10]={34.7649,38.1048,34.8069,34.1828,35.5456,33.7855,30,30,30,30,};
double par3[10]={43.9355,20.0796,32.3097,38.8749,35.2414,41.9533,50,50,50,50,};



  TChain myTree("demo/Tree");
















/*
myTree.Add("/tmp/bbilin/2011A_aug5_11_03_2012.root");
myTree.Add("/tmp/bbilin/2011A_may10_11_03_2012.root");
myTree.Add("/tmp/bbilin/2011A_promtv4_11_03_2012.root");
myTree.Add("/tmp/bbilin/2011A_promtv6_11_03_2012.root");
*/
//TFile *theFile = new TFile("06_03_2011A_4.root","RECREATE");

myTree.Add("/tmp/bbilin/data_2011B_Nov19_15_03_2012_44X.root");
TFile *theFile = new TFile("16_03_2011B_44X.root","RECREATE");

//myTree.Add("/tmp/bbilin/mc_mdgrph_04_03_2012.root");
//TFile *theFile = new TFile("06_03_mc_4.root","RECREATE");

//myTree.Add("/tmp/bbilin/WW_bg_12_03_2012.root");
//TFile *theFile = new TFile("12_03_bg_WW.root","RECREATE");

//myTree.Add("/tmp/bbilin/ZZ_bg_12_03_2012.root");
//TFile *theFile = new TFile("12_03_bg_ZZ.root","RECREATE");

//myTree.Add("/tmp/bbilin/WZ_bg_12_03_2012.root");
//TFile *theFile = new TFile("12_03_bg_WZ.root","RECREATE");

//myTree.Add("/tmp/bbilin/Ztautau_bg_12_03_2012.root");
//TFile *theFile = new TFile("12_03_bg_Ztautau.root","RECREATE");

//myTree.Add("/tmp/bbilin/Wjets_bg_12_03_2012.root");
//TFile *theFile = new TFile("12_03_bg_Wjets.root","RECREATE");

//myTree.Add("/tmp/bbilin/tt_bg_12_03_2012.root");
//TFile *theFile = new TFile("12_03_bg_tt.root","RECREATE");

TH1::AddDirectory(true);

  int event,run,lumi,bxnumber,realdata;
//  int hlt_trigger_fired;
  int Elecindex;
  int ElecIsEE[50],ElecIsEB[50];
  float ElecPt[50], ElecEta[50], ElecPhi[50],ElecPx[50], ElecPy[50], ElecPz[50], ElecM[50], ElecE[50], Elecdr03TkSumPt[50], Elecdr03EcalRecHitSumEt[50], Elecdr03HcalTowerSumEt[50], ElecscSigmaIEtaIEta[50], ElecdeltaPhiSuperClusterTrackAtVtx[50], ElecdeltaEtaSuperClusterTrackAtVtx[50], ElechadronicOverEm[50],ElecgsfTrack_numberOfLostHits[50],ElecGsfTrk_d0[50],ElecDcotTheta[50];


  //particle information
  int par_index, mom[50], daug[50];
  float ParticlePt[50], ParticleEta[50], ParticlePhi[50], ParticlePx[50], ParticlePy[50], ParticlePz[50], ParticleE[50], ParticleM[50];
  int ParticleId[50], ParticleStatus[50], ParticleMother[50][10], ParticleDaughter[50][10];
  //  int id_elec[20];
  int ElecCharge[50];
  //reco jets 
  int CaloNjets;

  int techTrigger[50],nGoodVertices;

 
  float CalojetEta[50], CalojetPhi[50],CalojetPt[50],CalojetPx[50],CalojetPy[50],CalojetPz[50],CalojetE[50],CaloHadEHF[100],CaloEmEHF[100];
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
//  myTree.SetBranchAddress("hlt_trigger_fired",&hlt_trigger_fired);
  myTree.SetBranchAddress("Elecindex",&Elecindex);
  myTree.SetBranchAddress("ElecPt",ElecPt);
  myTree.SetBranchAddress("ElecGsfTrk_d0",ElecGsfTrk_d0);
  myTree.SetBranchAddress("ElecDcotTheta",ElecDcotTheta);
  myTree.SetBranchAddress("ElecEta",ElecEta);
  myTree.SetBranchAddress("ElecPhi",ElecPhi);
  myTree.SetBranchAddress("ElecPx",ElecPx);
  myTree.SetBranchAddress("ElecPy",ElecPy);
  myTree.SetBranchAddress("ElecPz",ElecPz);
  myTree.SetBranchAddress("ElecE",ElecE);
  myTree.SetBranchAddress("ElecM",ElecM);
  myTree.SetBranchAddress("ElecCharge",ElecCharge);
  myTree.SetBranchAddress("Elecdr03TkSumPt", Elecdr03TkSumPt);
  myTree.SetBranchAddress("Elecdr03EcalRecHitSumEt", Elecdr03EcalRecHitSumEt);
  myTree.SetBranchAddress("Elecdr03HcalTowerSumEt", Elecdr03HcalTowerSumEt);
  myTree.SetBranchAddress("ElecscSigmaIEtaIEta", ElecscSigmaIEtaIEta);
  myTree.SetBranchAddress("ElecdeltaPhiSuperClusterTrackAtVtx", ElecdeltaPhiSuperClusterTrackAtVtx);
  myTree.SetBranchAddress("ElecdeltaEtaSuperClusterTrackAtVtx", ElecdeltaEtaSuperClusterTrackAtVtx);
  myTree.SetBranchAddress("ElechadronicOverEm", ElechadronicOverEm);
  myTree.SetBranchAddress("ElecgsfTrack_numberOfLostHits", ElecgsfTrack_numberOfLostHits);
  myTree.SetBranchAddress("ElecIsEB", ElecIsEB);
  myTree.SetBranchAddress("ElecIsEE", ElecIsEE);

  myTree.SetBranchAddress("techTrigger",techTrigger);

  
  myTree.SetBranchAddress("par_index", &par_index);
  myTree.SetBranchAddress("ParticlePt", ParticlePt);
  myTree.SetBranchAddress("ParticleEta", ParticleEta);
  myTree.SetBranchAddress("ParticlePhi", ParticlePhi);
  myTree.SetBranchAddress("ParticlePx", ParticlePx);
  myTree.SetBranchAddress("ParticlePy", ParticlePy);
  myTree.SetBranchAddress("ParticlePz", ParticlePz);
  myTree.SetBranchAddress("ParticleE", ParticleE);
  myTree.SetBranchAddress("ParticleM", ParticleM);  
  myTree.SetBranchAddress("ParticleId", ParticleId);
  myTree.SetBranchAddress("mom", mom);
  myTree.SetBranchAddress("daug", daug);
  myTree.SetBranchAddress("ParticleStatus", ParticleStatus);
  myTree.SetBranchAddress("ParticleMother", ParticleMother);
  myTree.SetBranchAddress("ParticleDaughter", ParticleDaughter);

myTree.SetBranchAddress("nGoodVertices", &nGoodVertices);
  
  
  //recojets

 



  myTree.SetBranchAddress("CaloNjets",&CaloNjets);
  myTree.SetBranchAddress("CalojetEta",CalojetEta);
  myTree.SetBranchAddress("CalojetPhi",CalojetPhi);
  myTree.SetBranchAddress("CalojetPt",CalojetPt);
  myTree.SetBranchAddress("CalojetPx",CalojetPx);
  myTree.SetBranchAddress("CalojetPy",CalojetPy);
  myTree.SetBranchAddress("CalojetPz",CalojetPz);
  myTree.SetBranchAddress("CalojetE",CalojetE);
  myTree.SetBranchAddress("CaloHadEHF",CaloHadEHF);
  myTree.SetBranchAddress("CaloEmEHF",CaloEmEHF);



  myTree.SetBranchAddress("nVertices",&nVertices);

//TFile *theFile = new TFile("06_03_2011A_4.root","RECREATE");

//TFile *theFile = new TFile("06_03_mc_4.root","RECREATE");

//TFile *theFile = new TFile("12_03_bg_WW.root","RECREATE");

//TFile *theFile = new TFile("12_03_bg_ZZ.root","RECREATE");

//TFile *theFile = new TFile("12_03_bg_WZ.root","RECREATE");

//TFile *theFile = new TFile("12_03_bg_Ztautau.root","RECREATE");

//TFile *theFile = new TFile("12_03_bg_Wjets.root","RECREATE");

//TFile *theFile = new TFile("12_03_bg_tt.root","RECREATE");

  theFile->cd();


TH1F *h_deltaphi_z_jet=  new TH1F ("h_deltaphi_z_jet","h_deltaphi_z_jet",35,0,3.5);
TH1F *h_elec_pt=  new TH1F ("h_elec_pt","h_elec_pt",70,0,200);
TH1F *h_elec_eta=  new TH1F ("h_elec_eta","h_elec_eta",16,-8,8);
TH1F *h_leading_jet_pt=  new TH1F ("h_leading_jet_pt","h_leading_jet_pt",80,0,200);
TH1F *h_leading_jet_eta=  new TH1F ("h_leading_jet_eta","h_leading_jet_eta",40,-6,6);
TH1F *h_mz=  new TH1F ("h_mz","h_mz",50,70,110);
TH1F *h_dielecphi=  new TH1F ( "h_dielecphi", "h_dielecphi",20,-3,3);
TH1F *h_dielec_PT=  new TH1F ( "h_dielec_PT", "h_dielec_PT",100,0,200);
TH1F *h_dielec_rapidity=  new TH1F ( "h_dielec_rapidity", "h_dielec_rapidity",50,-5,5);

TH1F *h_deltar_elec_Calojet=  new TH1F ("h_deltar_elec_Calojet","h_deltar_elec_Calojet",100,0,2);
TH1F *h_deltar_parton_jet=  new TH1F ("h_deltar_parton_jet","h_deltar_parton_jet",100,0,2);
TH1F *h_parton_pt=  new TH1F ("h_parton_pt","h_parton_pt",100,0,200);
TH1F *h_parton_eta=  new TH1F ("h_parton_eta","h_parton_eta",50,-3.,3.);
TH1F *h_parton_phi=  new TH1F ("h_parton_phi","h_parton_phi",50,-3.,3.);

TH1F *h_quark_pt=  new TH1F ("h_quark_pt","h_quark_pt",100,0,200);
TH1F *h_quark_eta=  new TH1F ("h_quark_eta","h_quark_eta",50,-3.,3.);
TH1F *h_quark_phi=  new TH1F ("h_quark_phi","h_quark_phi",50,-3.,3.);

TH1F *h_gluon_pt=  new TH1F ("h_gluon_pt","h_gluon_pt",100,0,200);
TH1F *h_gluon_eta=  new TH1F ("h_gluon_eta","h_gluon_eta",50,-3.,3.);
TH1F *h_gluon_phi=  new TH1F ("h_gluon_phi","h_gluon_phi",50,-3.,3.);

TH1F * h_jet_parton=  new TH1F ("h_jet_parton","P_{T}^{Jet}/P_{T}^{Parton}",50,-0.5,3.5);
TH1F * h_jet_gluon=  new TH1F ("h_jet_gluon","P_{T}^{Jet}/P_{T}^{Gluon}",50,-0.5,3.5);
TH1F * h_jet_quark=  new TH1F ("h_jet_quark","P_{T}^{Jet}/P_{T}^{Quark}",50,-0.5,3.5);
TProfile* h_ratio_jetpt_partonpt = new TProfile("h_ratio_jetpt_partonpt","",10,0,200);
TProfile* h_ratio_gjetpt_partonpt = new TProfile("h_ratio_gjetpt_partonpt","",10,0,200);
TProfile* h_ratio_qjetpt_partonpt = new TProfile("h_ratio_qjetpt_partonpt","",10,0,200);
TProfile* h_ratio_corrjetpt_partonpt = new TProfile("h_ratio_corrjetpt_partonpt","",10,0,200);
TProfile* h_ratio_corrgjetpt_partonpt = new TProfile("h_ratio_corrgjetpt_partonpt","",10,0,200);
TProfile* h_ratio_corrqjetpt_partonpt = new TProfile("h_ratio_corrqjetpt_partonpt","",10,0,200);




TH1F * h_number_of_st3_prt =  new TH1F ("h_number_of_st3_prt","h_number_of_st3_prt",30,-0.5,29.5);

TH1F * h_nlo_jet_pt =  new TH1F ("h_nlo_jet_pt","nlo_jet_pt",120,0,300);

TH1F * h_mz_parton=  new TH1F ("h_mz_parton","mz_parton",120,0,300);

TH1F * n_particle_after_sellecting_goodjet = new TH1F ("n_particle_after_sellecting_goodjet","n_particle_after_sellecting_goodjet",20,-0.5,19.5);

 int NbinsX = 10;
  int NbinsY = 10;
  float etamin = 0.;
  float etamax = 5.191;
  
TH2F* h_ptz_etajet = new TH2F("h_ptz_etajet","",10,0,200,10,0,5.191);
  h_ptz_etajet->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptz_etajet->GetYaxis()->SetTitle("#eta(j)");
TH2F* h_ptjet_etajet = new TH2F("h_ptjet_etajet","",10,0,200,10,0,5.191);
  h_ptjet_etajet->GetXaxis()->SetTitle("p_{T}(j) (GeV)");
  h_ptjet_etajet->GetYaxis()->SetTitle("#eta(j)");
TH2F* h_ptz_etajet_g = new TH2F("h_ptz_etajet_g","",10,0,200,10,0,5.191);
  h_ptz_etajet_g->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptz_etajet_g->GetYaxis()->SetTitle("#eta(j)");
TH2F* h_ptz_etajet_q = new TH2F("h_ptz_etajet_q","",10,0,200,10,0,5.191);
  h_ptz_etajet_q->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptz_etajet_q->GetYaxis()->SetTitle("#eta(j)");
 TProfile2D* h_ptZ_etajet__average_PTZ = new TProfile2D("h_ptZ_etajet__average_PTZ","<P_{T}(Z)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_PTZ->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_PTZ->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_PTZ_g = new TProfile2D("h_ptZ_etajet__average_PTZ_g","<P_{T}(Z)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_PTZ_g->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_PTZ_g->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_PTZ_q = new TProfile2D("h_ptZ_etajet__average_PTZ_q","<P_{T}(Z)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_PTZ_q->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_PTZ_q->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_ETraw = new TProfile2D("h_ptZ_etajet__average_ETraw","<E_{T}(j)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_ETraw->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_ETraw->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_ETraw_g = new TProfile2D("h_ptZ_etajet__average_ETraw_g","<E_{T}(j)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_ETraw_g->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_ETraw_g->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_ETraw_q = new TProfile2D("h_ptZ_etajet__average_ETraw_q","<E_{T}(j)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_ETraw_q->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_ETraw_q->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_ratio = new TProfile2D("h_ptZ_etajet__average_ratio","<E_{T}(j)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_ratio->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_ratio->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptjet_etajet__average_ratio = new TProfile2D("h_ptjet_etajet__average_ratio","<E_{T}(j)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptjet_etajet__average_ratio->GetXaxis()->SetTitle("p_{T}(jet) (GeV)");
  h_ptjet_etajet__average_ratio->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_ratio_g = new TProfile2D("h_ptZ_etajet__average_ratio_g","<E_{T}(j)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_ratio_g->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_ratio_g->GetYaxis()->SetTitle("#eta(j)");

  TProfile2D* h_ptZ_etajet__average_ratio_q = new TProfile2D("h_ptZ_etajet__average_ratio_q","<E_{T}(j)>",NbinsX,0,200,NbinsY,etamin,etamax);
  h_ptZ_etajet__average_ratio_q->GetXaxis()->SetTitle("p_{T}(Z) (GeV)");
  h_ptZ_etajet__average_ratio_q->GetYaxis()->SetTitle("#eta(j)");


TH1F *h_ratio_m[100];

TH1F *h_ratio_p[100];

TH1F *h_ratio_all_m[100];

TH1F *h_ratio_all_p[100];

TH1F *h_invratio_m[100];

TH1F *h_invratio_p[100];

TH1F *h_invratio_all_m[100];

TH1F *h_invratio_all_p[100];

TH1F *h_ratio_parton_m[100];

TH1F *h_ratio_parton_p[100];

TH1F *h_ratio_parton_all_m[100];

TH1F *h_ratio_parton_all_p[100];

TH1F *h_invratio_parton_m[100];

TH1F *h_invratio_parton_p[100];

TH1F *h_invratio_parton_all_m[100];

TH1F *h_invratio_parton_all_p[100];

char name1[100];
char name2[100];
char name3[100];
char name4[100];
char name5[100];
char name6[100];
char name7[100];
char name8[100];
char name9[100];
char name10[100];
char name11[100];
char name12[100];
char name13[100];
char name14[100];
char name15[100];
char name16[100];

for(int i = 0; i<3;i++){

sprintf(name1,"h_ratio_m%i",i);
sprintf(name2,"h_ratio_p%i",i);
sprintf(name5,"h_invratio_m%i",i);
sprintf(name6,"h_invratio_p%i",i);
sprintf(name9,"h_ratio_m%i_parton",i);
sprintf(name10,"h_ratio_p%i_parton",i);
sprintf(name11,"h_invratio_m%i_parton",i);
sprintf(name12,"h_invratio_p%i_parton",i);
h_ratio_m[i] = new TH1F (name1,"",60,-0.5,5.5);
h_ratio_p[i] = new TH1F (name2,"",60,-0.5,5.5);
h_invratio_m[i] = new TH1F (name5,"",60,-0.5,5.5);
h_invratio_p[i] = new TH1F (name6,"",60,-0.5,5.5);
h_ratio_parton_m[i] = new TH1F (name9,"",60,-0.5,5.5);
h_ratio_parton_p[i] = new TH1F (name10,"",60,-0.5,5.5);
h_invratio_parton_m[i] = new TH1F (name11,"",60,-0.5,5.5);
h_invratio_parton_p[i] = new TH1F (name12,"",60,-0.5,5.5);
}

for(int i = 0; i<13;i++){

sprintf(name3,"h_ratio_all_m%i",i);
sprintf(name4,"h_ratio_all_p%i",i);
sprintf(name7,"h_invratio_all_m%i",i);
sprintf(name8,"h_invratio_all_p%i",i);
sprintf(name13,"h_ratio_m%i_parton_all",i);
sprintf(name14,"h_ratio_p%i_parton_all",i);
sprintf(name15,"h_invratio_m%i_parton_all",i);
sprintf(name16,"h_invratio_p%i_parton_all",i);
h_ratio_all_m[i] = new TH1F (name3,"",60,-0.5,5.5);
h_ratio_all_p[i] = new TH1F (name4,"",60,-0.5,5.5);
h_invratio_all_m[i] = new TH1F (name7,"",60,-0.5,5.5);
h_invratio_all_p[i] = new TH1F (name8,"",60,-0.5,5.5);
h_ratio_parton_all_m[i] = new TH1F (name13,"",60,-0.5,5.5);
h_ratio_parton_all_p[i] = new TH1F (name14,"",60,-0.5,5.5);
h_invratio_parton_all_m[i] = new TH1F (name15,"",60,-0.5,5.5);
h_invratio_parton_all_p[i] = new TH1F (name16,"",60,-0.5,5.5);
}

  double xbinplus[11] = {2.853,2.946,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538};

  TH1* CorrVsEtaPlus  = new TProfile("CorrVsEtaPlus","pT(Z)/pT(jet) for HF+",10,xbinplus,0,5);
  CorrVsEtaPlus->Sumw2();
  TH1* CorrVsEtaMinus  = new TProfile("CorrVsEtaMinus","pT(Z)/pT(jet) for HF-",10,xbinplus,0,5);
  CorrVsEtaMinus->Sumw2();





 Int_t nevent = myTree.GetEntries();
 //Int_t  nevent = 4000000;
  for (Int_t iev=0;iev<nevent;iev++) {




 	  if (iev%100000 == 0) cout<<iev<<"/"<<nevent<<endl;
    	myTree.GetEntry(iev); 

if(realdata && tr_1!=1 && tr_2!=1 && tr_3!=1) continue;




int eCharge[50]={};
int eIsEE[50]={};
int eIsEB[50]= {};
double ePt[50] = {};
double eEta[50] = {};
double ePhi[50]= {};
double ePx[50]={};
double ePy[50]={};
double ePz[50] = {};
double eM[50]= {}; 
double eE[50]={}; 
double edr03TkSumPt[50]={}; 
double edr03EcalRecHitSumEt[50]={}; 
double edr03HcalTowerSumEt[50]= {}; 
double escSigmaIEtaIEta[50]={}; 
double edeltaPhiSuperClusterTrackAtVtx[50]={};
double edeltaEtaSuperClusterTrackAtVtx[50]= {}; 
double ehadronicOverEm[50]={} ;
double egsfTrack_numberOfLostHits[50]={} ;
double eGsfTrk_d0[50]={} ;
double eDcotTheta[50]={};
int e_cut_fired[50]={};



//cout<<"new"<<endl;
for(int ii = 0 ; ii < Elecindex ; ii++)
{

eCharge[ii]=ElecCharge[ii];
eIsEE[ii]=ElecIsEE[ii];
eIsEB[ii]= ElecIsEB[ii];
//cout<<ePt[ii]<<endl;
eEta[ii] = ElecEta[ii];
ePhi[ii]= ElecPhi[ii];
//if(!realdata){

ePt[ii] = ElecPt[ii];
ePx[ii]= ElecPx[ii];
ePy[ii]=  ElecPy[ii];
ePz[ii] =ElecPz[ii];

//}

/*
if(realdata){
for(int i = 0 ; i < 9 ; i++){
if((i-5)*0.5<eEta[ii] && eEta[ii]<(i-4)*0.5){
//cout<<eEta[ii]<<"     "<<i<<"    "<< c_promtv4[i]<<endl;


if(160430<run && run<163870){
ePt[ii] = c_may10[i]*ElecPt[ii];
ePx[ii]= c_may10[i]*ElecPx[ii];
ePy[ii]= c_may10[i]* ElecPy[ii];
ePz[ii] = c_may10[i]* ElecPz[ii];
//cout<<"1"<<endl;
}

if(165087<run && run<167914){
ePt[ii] = c_promtv4[i]*ElecPt[ii];
ePx[ii]= c_promtv4[i]*ElecPx[ii];
ePy[ii]= c_promtv4[i]* ElecPy[ii];
ePz[ii] = c_promtv4[i]* ElecPz[ii];
//cout<<"2"<<endl;
}

if(170825<run && run<172620){
ePt[ii] = c_aug5[i]*ElecPt[ii];
ePx[ii]= c_aug5[i]*ElecPx[ii];
ePy[ii]= c_aug5[i]* ElecPy[ii];
ePz[ii] = c_aug5[i]* ElecPz[ii];
//cout<<"3"<<endl;
}

if(172619<run && run<173693){
ePt[ii] = c_promtv6[i]*ElecPt[ii];
ePx[ii]= c_promtv6[i]*ElecPx[ii];
ePy[ii]= c_promtv6[i]* ElecPy[ii];
ePz[ii] = c_promtv6[i]* ElecPz[ii];
//cout<<"4"<<endl;
}
}}}
*/

eM[ii]= ElecM[ii]; 
eE[ii]= ElecE[ii]; 
edr03TkSumPt[ii]= Elecdr03TkSumPt[ii]; 
edr03EcalRecHitSumEt[ii]= Elecdr03EcalRecHitSumEt[ii]; 
edr03HcalTowerSumEt[ii]= Elecdr03HcalTowerSumEt[ii]; 
escSigmaIEtaIEta[ii]= ElecscSigmaIEtaIEta[ii]; 
edeltaPhiSuperClusterTrackAtVtx[ii]=ElecdeltaPhiSuperClusterTrackAtVtx[ii];
edeltaEtaSuperClusterTrackAtVtx[ii]= ElecdeltaEtaSuperClusterTrackAtVtx[ii]; 
ehadronicOverEm[ii]= ElechadronicOverEm[ii];
egsfTrack_numberOfLostHits[ii]= ElecgsfTrack_numberOfLostHits[ii];
eGsfTrk_d0[ii]= ElecGsfTrk_d0[ii];
eDcotTheta[ii]=ElecDcotTheta[ii];





if(
egsfTrack_numberOfLostHits[ii]<1 && ( fabs(eGsfTrk_d0[ii])>0.02 || fabs(eDcotTheta[ii])>0.02 ) &&  ePt[ii] > 20.  && fabs(eEta[ii]) <2.4 &&

((eIsEB[ii] ==1 && eIsEB[ii]==1 && edr03TkSumPt[ii]/ePt[ii]<0.09 && edr03EcalRecHitSumEt[ii]/ePt[ii]<0.07 && edr03HcalTowerSumEt[ii]/ePt[ii]<0.10 && escSigmaIEtaIEta[ii]<0.01 && fabs(edeltaPhiSuperClusterTrackAtVtx[ii])<0.06 && fabs(edeltaEtaSuperClusterTrackAtVtx[ii])<0.004 && ehadronicOverEm[ii]<0.04 
) ||
(eIsEE[ii] ==1 && edr03TkSumPt[ii]/ePt[ii]<0.04 && edr03EcalRecHitSumEt[ii]/ePt[ii]<0.05 && edr03HcalTowerSumEt[ii]/ePt[ii]<0.025 && escSigmaIEtaIEta[ii]<0.03 && fabs(edeltaPhiSuperClusterTrackAtVtx[ii])<0.03 && fabs(edeltaEtaSuperClusterTrackAtVtx[ii])<0.007 && ehadronicOverEm[ii]<0.025 ))

){e_cut_fired[ii]=1;}

else {e_cut_fired[ii]=0;}

//cout<<e_cut_fired[ii]<<endl;
//if(e_cut_fired[ii]==1) cout<<"  no of lost hits  "<<egsfTrack_numberOfLostHits[ii]<<"  d0   "<<fabs(eGsfTrk_d0[ii])<<"   dcottheta   "<<fabs(eDcotTheta[ii])<<"   ept   "<<ePt[ii]<<"   eeta   "<<fabs(eEta[ii])<<endl;
}
//cout<<"aaa"<<endl;


double jetEta[50]={};
double jetPhi[50]={};
double jetPt[50]={};
double jetPx[50]={};
double jetPy[50]={};
double jetPz[50]={};
double jetE[50]={};
double EmEHF[50]={};
double HadEHF[50]={};
int sorted_jet_index[50]={};

double jPt[50]={};

for(int ii=0;ii<CaloNjets;ii++){
jPt[ii]=0;
jetPt[ii]=0;
 jetEta[ii]=0;
 jetPhi[ii]=0;
 jetPt[ii]=0;
 jetPx[ii]=0;
 jetPy[ii]=0;
 jetPz[ii]=0;
 jetE[ii]=0;
 EmEHF[ii]=0;
 HadEHF[ii]=0;
}

for(int ii=0;ii<CaloNjets;ii++){
jPt[ii]=fabs(CalojetPt[ii]);
}

double corr_fact=-999.;
int sort_jet = 0;
int ind=0;
TMath::Sort(CaloNjets,jPt,sorted_jet_index);

for(int ii = 0 ; ii < CaloNjets ; ii++)
{

 sort_jet = sorted_jet_index[ind];

jetEta[ind]=CalojetEta[sort_jet]  ;
jetPhi[ind]=CalojetPhi[sort_jet] ;
jetPt[ind]= jPt[sort_jet];
if(jetPt[ind]<-10)cout<<"  "<<jetPt[0]<<"  "<<jetPt[1]<<"  "<<jetPt[ind]<<endl;
jetPx[ind]= CalojetPx[sort_jet];
jetPy[ind]= CalojetPy[sort_jet];
jetPz[ind]= CalojetPz[sort_jet];
jetE[ind]= CalojetE[sort_jet];
HadEHF[ind]=CaloHadEHF[sort_jet];
EmEHF[ind]=CaloEmEHF[sort_jet];
ind++;
}



    float particle_pt[50] = {};
    float particle_eta[50] = {};
    float particle_phi[50] = {};
    float particle_px[50] = {};
    float particle_py[50] = {};
    float particle_pz[50] = {};
    float particle_m[50] = {};
    float particle_e[50] = {};
    int particle_id[50] = {};


//elecs


    int nparticle_gen =0;



    for (int ppo = 0;ppo<50;++ppo) particle_id[ppo] = 9999; 
    for (int j = 0;j<par_index;++j){
      if (ParticleStatus[j] != 3) continue;
      particle_pt[nparticle_gen] = ParticlePt[j];
      particle_eta[nparticle_gen] = ParticleEta[j];
      particle_phi[nparticle_gen] = ParticlePhi[j];
      particle_id[nparticle_gen] = ParticleId[j];
      particle_px[nparticle_gen] = ParticlePx[j];
      particle_py[nparticle_gen] = ParticlePy[j];
      particle_pz[nparticle_gen] = ParticlePz[j];
      particle_m[nparticle_gen] = ParticleM[j];
      particle_e[nparticle_gen] = ParticleE[j];

//	cout<<nparticle_gen<<"   "<<particle_id[nparticle_gen]<<endl;

nparticle_gen++ ;//9 or 10 only
}
if((nparticle_gen==10 )){

h_parton_pt->Fill(particle_pt[7]);
h_parton_eta->Fill(particle_eta[7]);
h_parton_phi->Fill(particle_phi[7]);

if(fabs(particle_id[7])<6){
h_quark_pt->Fill(particle_pt[7]);
h_quark_eta->Fill(particle_eta[7]);
h_quark_phi->Fill(particle_phi[7]);
}
if(fabs(particle_id[7])==21){
h_gluon_pt->Fill(particle_pt[7]);
h_gluon_eta->Fill(particle_eta[7]);
h_gluon_phi->Fill(particle_phi[7]);
}
}

//cout<<"neww"<<endl;
int pr_elecindex_1=-999;
int pr_elecindex_2=-999;
for(int ii=0;ii<nparticle_gen;ii++){
if(particle_id[ii]==11)pr_elecindex_1=ii;
if(particle_id[ii]==-11)pr_elecindex_2=ii;
//cout<<nparticle_gen<<"   "<<particle_id[ii]<<endl;
}

double p1dp2 = (pow(particle_px[pr_elecindex_1]+particle_px[pr_elecindex_2],2) + pow(particle_py[pr_elecindex_1]+particle_py[pr_elecindex_2],2) + pow(particle_pz[pr_elecindex_1]+particle_pz[pr_elecindex_2],2));
double mZ_particle = fabs(sqrt(pow(particle_e[pr_elecindex_1]+particle_e[pr_elecindex_2],2)-p1dp2 ));
h_mz_parton->Fill(mZ_particle);

h_number_of_st3_prt->Fill(nparticle_gen);
//cout<<nparton_gen;

//float ptcut = 20.;
    float p1dotp2 = -99.;
    float dielecpt = -999.;
    float dielecrapidity = -999.;
    float dielecphi = -999.;
    float MZelec  = -999. ;
	int select=0;
	int index1 =0;
	int index2 =0;

    for (int j = 0; j < Elecindex; ++j){  
      for (int jk = j; jk < Elecindex  ; ++jk){

	if (j == jk) continue;
if ((fabs(eEta[j]) >1.4442 && fabs(eEta[j])<1.560)|| (fabs(eEta[jk]) >1.4442 && fabs(eEta[jk])<1.560)) continue;

 	if (
 (eCharge[j]*eCharge[jk] == -1) && e_cut_fired[j]==1 && e_cut_fired[jk]==1){

if(select==1)continue;

index1=j;
index2=jk;
	h_elec_pt->Fill(ePt[j]);
	h_elec_pt->Fill(ePt[jk]);
	h_elec_eta->Fill(eEta[j]);
	h_elec_eta->Fill(eEta[jk]);

	       p1dotp2 = ePx[j]*ePx[jk]+ePy[j]*ePy[jk]+ePz[j]*ePz[jk];
               MZelec = eM[j]+eM[jk]+2*(eE[j]*eE[jk]-p1dotp2);
		
               if (MZelec < 0.) MZelec = -1.*sqrt(fabs(MZelec));
               if (MZelec > 0.) MZelec = sqrt(MZelec);

		h_mz->Fill(MZelec);
 //              h_elec_leading_PT->Fill(TMath::Max(ePt[j],ePt[jk]));

                  dielecpt = sqrt(pow((ePx[j]+ePx[jk]),2) +pow((ePy[j]+ePy[jk]),2));
//cout<<dielecpt<<endl;
                  dielecrapidity = 0.5*log((eE[j]+eE[jk]+ePz[j]+ePz[jk])/(eE[j]+eE[jk]-ePz[j]-ePz[jk]));
                  dielecphi = atan2(ePy[j]+ePy[jk],ePx[j]+ePx[jk]);
//		 h_elec_leading_eta->Fill(eEta[j]);

		select=1;
	    }
      }
   }

if (select==0) continue;


//      h_Mass_Z->Fill(MZelec);
//      h_Mass_Z_range->Fill(MZelec);
      h_dielecphi->Fill(dielecphi);
      h_dielec_PT->Fill(dielecpt);
      h_dielec_rapidity->Fill(dielecrapidity);


//      if (MZelec < 60 || MZelec > 120) continue;

  if (MZelec < 70 || MZelec > 110/* || dielecpt<20.*/ ) continue;
  
      
            int numberofCalojets = 0;
      int jetindex1 =-999; 
//      int in=0;
      int nlo_index=0;
      int markgoodev=0;
      int check2=0;
      double deltaphi = 999.0;
      double ratio0 = -999.0;
      double ratio1 = -999.0;
      double ratio2 = -999.0;
      double ratio3 = -999.0;

      int markCalojet;
      for (int Calo=0;Calo < CaloNjets ;Calo++){
           markCalojet = 0;

	  float deltar = fabs(DeltaR(jetEta[Calo],eEta[index1],jetPhi[Calo],ePhi[index1]));
	  float deltar2 = fabs(DeltaR(jetEta[Calo],eEta[index2],jetPhi[Calo],ePhi[index2]));


	h_deltar_elec_Calojet -> Fill(deltar);
	h_deltar_elec_Calojet -> Fill(deltar2); //test the code for deltar
	  if (deltar<0.5) markCalojet = 1; 
	  if (deltar2<0.5) markCalojet = 1; 
	if (markCalojet != 1){
	if (jetPt[Calo] > 30.){
	  numberofCalojets += 1;}		

	if(nlo_index==1 && check2==0){
		check2=1;			
		if(jetPt[Calo]/dielecpt <0.2 && dielecpt>5.)markgoodev=1;
//		if(jetPt[Calo] <15.)markgoodev=1;
		h_nlo_jet_pt->Fill(jetPt[Calo]);
				}
	if(nlo_index==0){
jetindex1=Calo;
}
	nlo_index=1;
}}


double cjetE=0.;
double cjetPt=0.;
double cjetPx=0.;
double cjetPy=0.;
double cjetPz=0.;

for(int i=0; i<10;i++){//eta index
if(fabs(jetEta[jetindex1])>0.519*i && fabs(jetEta[jetindex1])<0.519*(i+1) ){
 corr_fact= par0[i]/(par1[i] + par2[i]* log(par3[i]*jetPt[jetindex1]));
//cout<<corr_fact<<endl;
cjetE=corr_fact*jetE[jetindex1];
cjetPx=corr_fact*jetPx[jetindex1];
cjetPy=corr_fact*jetPy[jetindex1];
cjetPz=corr_fact*jetPz[jetindex1];
cjetPt=corr_fact*jetPt[jetindex1];
}
}


//cout<< markgoodev<<endl;

//h_numberofCalojets->Fill(numberofCalojets);

if (markgoodev!=1) continue;
n_particle_after_sellecting_goodjet->Fill(nparticle_gen);


//cout<< jetindex1<<endl;

double eta_index[]= {2.853,2.946,3.139,3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191};
double eta_index_m[]= {-2.853,-2.946,-3.139,-3.314,-3.489,-3.664,-3.839,-4.013,-4.191,-4.363,-4.538,-4.716,-4.889,-5.191}; 
//double eta_index2[]= {3.314,3.489,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889,5.191}; 
//double eta_index_2m[]= {-3.314,-3.489,-3.664,-3.839,-4.013,-4.191,-4.363,-4.538,-4.716,-4.889,-5.191} ;
double eta_index_cb[]= {2.853,3.139,3.489,5.191} ;
double eta_index_cbm[]= {-2.853,-3.139,-3.489,-5.191} ;
//	double y[]={1,1.157,1.063,1.040,1.081,1.046,1.050,1.023,1.059,0.993,1.025,1,1};
//	double ym[]={1,1.157,1.063,1.040,1.081,1.046,1.050,1.023,1.059,0.993,1.025,1,1};




deltaphi = fabs(DeltaPhi(dielecphi,jetPhi[jetindex1]));
 
h_deltaphi_z_jet->Fill(deltaphi);

if (deltaphi<2.8) continue;

	h_leading_jet_pt->Fill(jetPt[jetindex1]);
	h_leading_jet_eta->Fill(jetEta[jetindex1]);

if(!(realdata) && nparticle_gen!=10 ) continue;

    ratio0 = dielecpt/jetPt[jetindex1];
    ratio1 = dielecpt/particle_pt[7];
    ratio2 =jetPt[jetindex1]/particle_pt[7];
    ratio3 =cjetPt/particle_pt[7];
//cout<<cjetPt[jetindex1]<<endl; 


double deltar_parton_jet=DeltaR(jetEta[jetindex1],particle_eta[7],jetPhi[jetindex1],particle_phi[7]);
h_deltar_parton_jet->Fill(deltar_parton_jet);
    h_ptz_etajet->Fill(dielecpt,fabs(jetEta[jetindex1]));
    h_ptjet_etajet->Fill(jetPt[jetindex1],fabs(jetEta[jetindex1]));
    h_ptZ_etajet__average_PTZ->Fill(dielecpt,fabs(jetEta[jetindex1]),dielecpt);
    h_ptZ_etajet__average_ETraw->Fill(dielecpt,fabs(jetEta[jetindex1]),jetPt[jetindex1]);
    h_ptZ_etajet__average_ratio->Fill(dielecpt,fabs(jetEta[jetindex1]),ratio0);
    h_ptjet_etajet__average_ratio->Fill(jetPt[jetindex1],fabs(jetEta[jetindex1]),ratio0);
    h_jet_parton->Fill(ratio2);
    h_ratio_jetpt_partonpt->Fill(particle_pt[7],ratio2);
    h_ratio_corrjetpt_partonpt->Fill(particle_pt[7],ratio3);


h_parton_pt->Fill(particle_pt[7]);
 if(fabs(particle_id[7])<6){//quark jets
h_quark_pt->Fill(particle_pt[7]);
h_jet_quark->Fill(ratio2);
h_ptz_etajet_q->Fill(dielecpt,fabs(jetEta[jetindex1]));
 h_ptZ_etajet__average_PTZ_q->Fill(dielecpt,fabs(jetEta[jetindex1]),dielecpt);
    h_ptZ_etajet__average_ETraw_q->Fill(dielecpt,fabs(jetEta[jetindex1]),jetPt[jetindex1]);
    h_ptZ_etajet__average_ratio_q->Fill(dielecpt,fabs(jetEta[jetindex1]),ratio0);
        h_ratio_qjetpt_partonpt->Fill(particle_pt[7],ratio2);
	h_ratio_corrqjetpt_partonpt->Fill(particle_pt[7],ratio3);	
}

 if(fabs(particle_id[7])==21){//gluon jets
h_gluon_pt->Fill(particle_pt[7]);
h_jet_gluon->Fill(ratio2);
h_ptz_etajet_g->Fill(dielecpt,fabs(jetEta[jetindex1]));
 h_ptZ_etajet__average_PTZ_g->Fill(dielecpt,fabs(jetEta[jetindex1]),dielecpt);
    h_ptZ_etajet__average_ETraw_g->Fill(dielecpt,fabs(jetEta[jetindex1]),jetPt[jetindex1]);
    h_ptZ_etajet__average_ratio_g->Fill(dielecpt,fabs(jetEta[jetindex1]),ratio0);
    h_ratio_gjetpt_partonpt->Fill(particle_pt[7],ratio2);
    h_ratio_corrgjetpt_partonpt->Fill(particle_pt[7],ratio3);
}


/* 
for(int i = 0; i<13;i++){
  if(jetEta[jetindex1]>eta_index[i] && jetEta[jetindex1]<eta_index[i+1]){
corr =  ((HadEHF[jetindex1]*0.5+EmEHF[jetindex1]) * y[i] +0.5*HadEHF[jetindex1])/jetE[jetindex1] ;
	jetPt[jetindex1] = corr *jetPt[jetindex1];
}
if(jetEta[jetindex1]<eta_index_m[i] && jetEta[jetindex1]>eta_index_m[i+1])){
corr =  ((HadEHF[jetindex1]*0.5+EmEHF[jetindex1]) * ym[i] +0.5*HadEHF[jetindex1])/jetE[jetindex1] ;
	jetPt[jetindex1] = corr *jetPt[jetindex1];
}}
for(int i = 0; i<10;i++){
if(jetEta[jetindex1]>eta_index2[i] && jetEta[jetindex1]<eta_index2[i+1]){
corr =  ((HadEHF[jetindex1]*0.5+EmEHF[jetindex1]) *(ratio0 - y[i+3]) +0.5*ratio0*HadEHF[jetindex1])/(0.5*HadEHF[jetindex1]) ;
}
if(jetEta[jetindex1]<eta_index2m[i] && jetEta[jetindex1]>eta_index2m[i+1]){
corr =  ((HadEHF[jetindex1]*0.5+EmEHF[jetindex1]) *(ratio0 - ym[i+3]) +0.5*ratio0*HadEHF[jetindex1])/(0.5*HadEHF[jetindex1]) ;
}}
*/



//cout<<"NEW"<<endl;




for(int i = 0; i<3;i++){

  if(jetEta[jetindex1]>eta_index_cb[i] && jetEta[jetindex1]<eta_index_cb[i+1]){
//cout<<"   "<<i<<"   "<<"jet eta    "<<jetEta[jetindex1]<<endl;
        h_ratio_p[i]->Fill(ratio0);
 	h_invratio_p[i]->Fill(1/ratio0);
	h_ratio_parton_p[i]->Fill(ratio1);
 	h_invratio_parton_p[i]->Fill(1/ratio1);

}
  if(jetEta[jetindex1]<eta_index_cbm[i] && jetEta[jetindex1]>eta_index_cbm[i+1]){
//cout<<"   "<<i<<"   "<<"jet eta    "<<jetEta[jetindex1]<<endl;
        h_ratio_m[i]->Fill(ratio0);
	h_invratio_m[i]->Fill(1/ratio0);
	h_ratio_parton_m[i]->Fill(ratio1);
 	h_invratio_parton_m[i]->Fill(1/ratio1);

}} 

for(int i = 0; i<13;i++){
  if(jetEta[jetindex1]>eta_index[i] && jetEta[jetindex1]<eta_index[i+1]){
        h_ratio_all_p[i]->Fill(ratio0);
	h_invratio_all_p[i]->Fill(1/ratio0);
	h_ratio_parton_all_p[i]->Fill(ratio1);
 	h_invratio_parton_all_p[i]->Fill(1/ratio1);
}

if(jetEta[jetindex1]<eta_index_m[i] && jetEta[jetindex1]>eta_index_m[i+1]){
 h_ratio_all_m[i]->Fill(ratio0);
 h_invratio_all_m[i]->Fill(1/ratio0);
	h_ratio_parton_all_m[i]->Fill(ratio1);
 	h_invratio_parton_all_m[i]->Fill(1/ratio1);
}
} 

//	  if(jetEta[jetindex1]>-5.191 && jetEta[jetindex1]<-2.853) CorrVsEtaMinus->Fill(fabs(jetEta[jetindex1]),ratio0);
//          if(jetEta[jetindex1]>2.853  && jetEta[jetindex1]<5.191) CorrVsEtaPlus->Fill(jetEta[jetindex1],ratio0);	
//     	  if(fabs(jetEta[jetindex1])>3 && fabs(jetEta[jetindex1])<5) h_z_pt_Et_jet_HF_all_1jet->Fill(ratio0);
//	  if(jetEta[jetindex1]>3 && jetEta[jetindex1]<5) h_z_pt_Et_jet_HF_plus_1jet->Fill(ratio0);
//	  if(jetEta[jetindex1]>-5 && jetEta[jetindex1]<-3) h_z_pt_Et_jet_HF_minus_1jet->Fill(ratio0);



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
