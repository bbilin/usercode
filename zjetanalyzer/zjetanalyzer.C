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

myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/data_2012A.root");
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/data_2012B.root");
TFile *theFile = new TFile("rootfiles/16_08_data.root","RECREATE");
 dataname = "DATA";


//MC
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_tarball.root");
 ismc =1;
TFile *theFile = new TFile("rootfiles/16_08_mc.root","RECREATE");
 dataname = "MC";
*/
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_10_50.root");
 ismc =1;
TFile *theFile = new TFile("rootfiles/16_08_mc_low.root","RECREATE");
 dataname = "MC_LOW";
*/

//TTbar
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_tt.root");
TFile *theFile = new TFile("rootfiles/16_08_bg_tt.root","RECREATE");
 dataname = "TTBAR";
*/

//WJets
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_WJets.root");
TFile *theFile = new TFile("rootfiles/16_08_bg_wjets.root","RECREATE");
 dataname = "WJets";
*/

//ZTauTau
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_tarball.root");
 istau =1;
TFile *theFile = new TFile("rootfiles/16_08_bg_ztautau.root","RECREATE");
 dataname = "ZTauTau";
*/

//WW
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_ww.root");
TFile *theFile = new TFile("rootfiles/16_08_ww.root","RECREATE");
 dataname = "WW";
*/

//WZ
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_wz.root");
TFile *theFile = new TFile("rootfiles/16_08_wz.root","RECREATE");
 dataname = "WZ";
*/

//ZZ
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/bg_zz.root");
TFile *theFile = new TFile("rootfiles/16_08_zz.root","RECREATE");
 dataname = "ZZ";
*/


TH1::AddDirectory(true);

  int event,run,lumi,bxnumber,realdata;
//  int hlt_trigger_fired;
  int Elecindex, HFElecindex;
  float ElecPt[50], ElecEta[50], ElecPhi[50],ElecPx[50], ElecPy[50], ElecPz[50], ElecM[50], ElecE[50],Elecdr03TkSumPt[50], Elecdr03EcalRecHitSumEt[50], Elecdr03HcalTowerSumEt[50], ElecscSigmaIEtaIEta[50], ElecdeltaPhiSuperClusterTrackAtVtx[50], ElecdeltaEtaSuperClusterTrackAtVtx[50], ElechadronicOverEm[50],ElecgsfTrack_numberOfLostHits[50],ElecGsfTrk_d0[50],ElecDcotTheta[50],ElecMVATrigId[50],ElecMVANonTrigId[50],hfelec_Pt[50],hfelec_Phi[50],hfelec_Eta[50], hfelec_L9[50], hfelec_S9[50], hfelec_L25[50], hfelec_L1[50];

double MyWeight;

  //particle information
  int par_index, mom[50], daug[50];
 
  int ParticleId[50], ParticleStatus[50], ParticleMother[50][10], ParticleDaughter[50][10];
  //  int id_elec[20];
  int ElecCharge[50];
  //reco jets 
  int pfNjets;

  int techTrigger[50],nGoodVertices;

 
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
  myTree.SetBranchAddress("Elecdr03TkSumPt", Elecdr03TkSumPt);
  myTree.SetBranchAddress("Elecdr03EcalRecHitSumEt", Elecdr03EcalRecHitSumEt);
  myTree.SetBranchAddress("Elecdr03HcalTowerSumEt", Elecdr03HcalTowerSumEt);
  myTree.SetBranchAddress("ElecscSigmaIEtaIEta", ElecscSigmaIEtaIEta);
  myTree.SetBranchAddress("ElecdeltaPhiSuperClusterTrackAtVtx", ElecdeltaPhiSuperClusterTrackAtVtx);
  myTree.SetBranchAddress("ElecdeltaEtaSuperClusterTrackAtVtx", ElecdeltaEtaSuperClusterTrackAtVtx);
  myTree.SetBranchAddress("ElechadronicOverEm", ElechadronicOverEm);
  myTree.SetBranchAddress("ElecgsfTrack_numberOfLostHits", ElecgsfTrack_numberOfLostHits);


myTree.SetBranchAddress("HFElecindex",&HFElecindex );
myTree.SetBranchAddress("hfelec_Pt",hfelec_Pt );

myTree.SetBranchAddress("hfelec_Eta",hfelec_Eta );
myTree.SetBranchAddress("hfelec_Phi",hfelec_Phi );
myTree.SetBranchAddress("hfelec_L9", hfelec_L9);
myTree.SetBranchAddress("hfelec_S9",hfelec_S9 );
myTree.SetBranchAddress("hfelec_L25", hfelec_L25);
myTree.SetBranchAddress("hfelec_L1",hfelec_L1 );


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

TH1F* h_nelec_fire_mva = new TH1F ("h_nelec_fire_mva","h_nelec_fire_mva",50,-0.5,49.5);
TH1F* h_no_good_vtx = new TH1F ("h_no_good_vtx","h_no_good_vtx",60,-0.5,59.5);
TH1F* h_no_good_vtx_noWeight = new TH1F ("h_no_good_vtx_noWeight","h_no_good_vtx_noWeight",60,-0.5,59.5);

TH1F *h_elec_mva=  new TH1F ("h_elec_mva","h_elec_mva",90,-1.5,1.5);
TH1F *h_elec_mva_pt25=  new TH1F ("h_elec_mva_pt25","h_elec_mva_pt25",90,-1.5,1.5);
TH1F *h_elec_mva_l_cut=  new TH1F ("h_elec_mva_l_cut","h_elec_mva_l_cut",90,-1.5,1.5);
TH1F *h_elec_mva_nl_cut=  new TH1F ("h_elec_mva_nl_cut","h_elec_mva_nl_cut",90,-1.5,1.5);
TH1F *h_elec_mva_l_mva=  new TH1F ("h_elec_mva_l","h_elec_mva_l_mva",90,-1.5,1.5);
TH1F *h_elec_mva_nl_mva=  new TH1F ("h_elec_mva_nl_mva","h_elec_mva_nl_mva",90,-1.5,1.5);

TH1F *h_hfeL1=  new TH1F ("h_hfeL1","h_hfeL1",120,0,1000);
TH1F *h_hfeL9=  new TH1F ("h_hfeL9","h_hfeL9",120,0,1000);
TH1F *h_hfeS9=  new TH1F ("h_hfeS9","h_hfeS9",120,0,1000);
TH1F *h_hfeL25=  new TH1F ("h_hfeL25","h_hfeL25",90,0,1000);
TH1F *h_hfecut1=  new TH1F ("h_hfecut1","h_hfecut1",90,0,10);
TH1F *h_hfecut2=  new TH1F ("h_hfecut2","h_hfecut2",90,0,10);

TH1F *h_hfePt=  new TH1F ("h_hfe_Pt","h_hfe_Pt",70,0,200);
TH1F *h_hfeEta=  new TH1F ("h_hfe_Eta","h_hfe_Eta",80,-6,6);

TH1F*h_no_good_vtx_zevents[12]; TH1F*h_no_good_vtx_noWeight_zevents[12]; TH1F *h_deltaphi_z_jet[12];TH1F *h_elec_pt[12];TH1F *h_elec_eta[12];TH1F *h_mz[12]; TH1F *h_mz_same_sign[12];TH1F *h_dielecphi[12];TH1F *h_dielec_PT[12];TH1F *h_dielec_rapidity[12];TH1F *h_mz_1j_all[12];TH1F *h_dielecphi_1j_all[12];TH1F *h_dielec_PT_1j_all[12];TH1F *h_dielec_rapidity_1j_all[12];TH1F *h_mz_1j_c[12];TH1F *h_dielecphi_1j_c[12];TH1F *h_dielec_PT_1j_c[12];TH1F *h_dielec_rapidity_1j_c[12];TH1F *h_mz_1j_hf[12];TH1F *h_dielecphi_1j_hf[12];TH1F *h_dielec_PT_1j_hf[12];TH1F *h_dielec_rapidity_1j_hf[12];TH1F *h_deltar_elec_PFjet[12];TH1F *h_deltar_elec_PFjet_hf[12];TH1F *h_numberofPFjets[12];TH1F *h_numberofPFjets_hbhe[12];TH1F *h_numberofPFjets_hf[12];TH1F *h_leading_jet_pt[12];TH1F *h_leading_jet_eta[12];TH1F *h_jet_pt_all[12];TH1F *h_jet_eta_all[12];TH1F *h_jet_pt_hbhe[12];TH1F *h_jet_eta_hbhe[12];TH1F *h_jet_pt_hf[12];TH1F *h_jet_eta_hf[12];

//char name[100][100];
//char name2[100][100][100];
TString name2[100];

//sprintf(name2[0],"h_no_good_vtx_zevents_%i",i);
name2[0] = "h_no_good_vtx_zevents_";
name2[1] ="h_no_good_vtx_noWeight_zevents_";
name2[2] ="h_deltaphi_z_jet_";
name2[3] ="h_elec_pt_";
name2[4] ="h_elec_eta_";
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
name2[30] ="h_jet_eta_all";
name2[31] ="h_jet_pt_hbhe_";
name2[32] ="h_jet_eta_hbhe";
name2[33] ="h_jet_pt_hf_";
name2[34] ="h_jet_eta_hf";
for(int i=0;i<12;i++){
for(int k=0;k<35;k++){
name2[k]+=(i);
}
h_no_good_vtx_zevents[i]=new TH1F (name2[0],name2[0],60,-0.5,59.5);
h_no_good_vtx_noWeight_zevents[i]=new TH1F (name2[1],name2[1],60,-0.5,59.5);
h_deltaphi_z_jet[i]=  new TH1F (name2[2],name2[2],35,0,3.5);
h_elec_pt[i]=  new TH1F (name2[3],name2[3],70,0,200);
h_elec_eta[i]=  new TH1F (name2[4],name2[4],80,-5,5);
h_mz[i]=  new TH1F (name2[5],name2[5],80,0,200);
h_mz_same_sign[i]=  new TH1F (name2[6],name2[6],80,0,200);
h_dielecphi[i]=  new TH1F ( name2[7], name2[7],20,-3,3);
h_dielec_PT[i]=  new TH1F ( name2[8], name2[8],100,0,200);
h_dielec_rapidity[i]=  new TH1F ( name2[9], name2[9],50,-5,5);
h_mz_1j_all[i]=  new TH1F ( name2[10], name2[10],80,0,200);
h_dielecphi_1j_all[i]=  new TH1F ( name2[11], name2[11],20,-3,3);
h_dielec_PT_1j_all[i]=  new TH1F ( name2[12], name2[12],100,0,200);
h_dielec_rapidity_1j_all[i]=  new TH1F ( name2[13], name2[13],50,-5,5);
h_mz_1j_c[i]=  new TH1F ( name2[14], name2[14],80,0,200);
h_dielecphi_1j_c[i]=  new TH1F ( name2[15], name2[15],20,-3,3);
h_dielec_PT_1j_c[i]=  new TH1F ( name2[16], name2[16],100,0,200);
h_dielec_rapidity_1j_c[i]=  new TH1F ( name2[17], name2[17],50,-5,5);
h_mz_1j_hf[i]=  new TH1F ( name2[18], name2[18],80,0,200);
h_dielecphi_1j_hf[i]=  new TH1F ( name2[19], name2[19],20,-3,3);
h_dielec_PT_1j_hf[i]=  new TH1F ( name2[20], name2[20],100,0,200);
h_dielec_rapidity_1j_hf[i]=  new TH1F ( name2[21], name2[21],50,-5,5);
h_deltar_elec_PFjet[i]=  new TH1F (name2[22], name2[22],100,0,2);
h_deltar_elec_PFjet_hf[i]=  new TH1F (name2[23], name2[23],100,0,2);
h_numberofPFjets[i] = new TH1F (name2[24], name2[24],14,-0.5,13.5);
h_numberofPFjets_hbhe[i] = new TH1F (name2[25], name2[25],14,-0.5,13.5);
h_numberofPFjets_hf[i] = new TH1F (name2[26], name2[26],14,-0.5,13.5);
h_leading_jet_pt[i]=  new TH1F (name2[27], name2[27],80,0,200);
h_leading_jet_eta[i]=  new TH1F (name2[28], name2[28],40,-6,6);
h_jet_pt_all[i]=  new TH1F (name2[29], name2[29],80,0,200);
h_jet_eta_all[i]=  new TH1F (name2[30], name2[30],40,-6,6);
h_jet_pt_hbhe[i]=  new TH1F (name2[31], name2[31],80,0,200);
h_jet_eta_hbhe[i]=  new TH1F (name2[32], name2[32],40,-6,6);
h_jet_pt_hf[i]=  new TH1F (name2[33], name2[33],80,0,200);
h_jet_eta_hf[i]=  new TH1F (name2[34], name2[34],40,-6,6);
}

 TH1F *h_short=  new TH1F ( "h_short", "h_short",100,0,200);
 TH1F *h_long=  new TH1F ( "h_long", "h_long",100,0,200);







 Int_t nevent = myTree.GetEntries();
//Int_t nevent =1000000;  

for (Int_t iev=0;iev<nevent;iev++) {

 	  if (iev%100000 == 0) cout<<dataname<<" ===>>> "<<iev<<"/"<<nevent<<endl;

    	myTree.GetEntry(iev); 

if(realdata && tr_3!=1) continue;

if (realdata) MyWeight=1; 

h_no_good_vtx-> Fill(nGoodVertices,MyWeight);
h_no_good_vtx_noWeight-> Fill(nGoodVertices);


int elec_ind_mva[50],elec_ind_cut[50],elec_ind_hf[50];

for(int i=0; i<50; i++){elec_ind_mva[i]=0;elec_ind_cut[i]=0; elec_ind_hf[i]=0;}

int eCharge[50]={},nelec_fire_mva=0,nelec_hf_fire_cut=0, e_cut_fired[50]={},nelec_fire_cut=0;
double ePt[50] = {}, eEta[50] = {}, ePhi[50]= {}, ePx[50]={}, ePy[50]={}, ePz[50] = {}, eM[50]= {}, eE[50]={}, eMVATrigId[50]={},edr03TkSumPt[50]={}, edr03EcalRecHitSumEt[50]={}, edr03HcalTowerSumEt[50]= {}, escSigmaIEtaIEta[50]={}, edeltaPhiSuperClusterTrackAtVtx[50]={}, edeltaEtaSuperClusterTrackAtVtx[50]= {}, ehadronicOverEm[50]={} , egsfTrack_numberOfLostHits[50]={} , eGsfTrk_d0[50]={} , eDcotTheta[50]={}; /*,eMVANonTrigId[50]={}*/; 
double jetEta[50]={}, jetPhi[50]={}, jetPt[50]={},jPt[50]={}, HadEHF[200]={}, EmEHF[200]={};
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
double cut1= 1/eE[ii] -(1/(sqrt(pow(ePt[ii],2)+pow(ePz[ii],2) ) ) );
if(
egsfTrack_numberOfLostHits[ii]<1 &&  /*fabs(eGsfTrk_d0[ii])>0.02 && fabs(eDcotTheta[ii])>0.02 && */  ePt[ii] > 20.  && fabs(eEta[ii]) <2.4 && (fabs(eEta[ii]) <1.4442 || fabs(eEta[ii])>1.566) &&(fabs(cut1)<0.05)&&
(
( fabs(eEta[ii]) <1.4442 &&escSigmaIEtaIEta[ii]<0.01 && fabs(edeltaPhiSuperClusterTrackAtVtx[ii])<0.06 && fabs(edeltaEtaSuperClusterTrackAtVtx[ii])<0.004 && ehadronicOverEm[ii]<0.12 
) ||
(fabs(eEta[ii])>1.566 &&escSigmaIEtaIEta[ii]<0.03 && fabs(edeltaPhiSuperClusterTrackAtVtx[ii])<0.03 && fabs(edeltaEtaSuperClusterTrackAtVtx[ii])<0.007 && ehadronicOverEm[ii]<0.10 
)
)
){e_cut_fired[ii]=1;
elec_ind_cut[nelec_fire_cut]=ii;
nelec_fire_cut++;
}
//eMVANonTrigId[ii]=ElecMVANonTrigId[ii];
h_elec_mva->Fill(eMVATrigId[ii],MyWeight);
if(ePt[ii] >25.)h_elec_mva_pt25->Fill(eMVATrigId[ii],MyWeight);
if(eMVATrigId[ii]>0 &&ePt[ii] >25.&&(fabs(eEta[ii]) <1.4442 || fabs(eEta[ii])>1.560)){
elec_ind_mva[nelec_fire_mva]=ii;
nelec_fire_mva++;
}


}

if(nelec_fire_cut==2 &&nelec_fire_mva==2 &&(elec_ind_cut[0]!=elec_ind_mva[0]) )cout<<nelec_fire_cut<<"    "<<nelec_fire_mva<<"   "<<elec_ind_cut[0]<<"  "<<elec_ind_cut[1]<<"   "<<elec_ind_mva[0]<<"  "<<elec_ind_mva[1] <<endl;

h_nelec_fire_mva->Fill(nelec_fire_mva);

double hfePt[50]={}, hfeEta[50]={},hfePhi[50]={},  hfeL9[50]={},  hfeS9[50]={},  hfeL25[50]={},  hfeL1[50]={};
//
for(int ii = 0 ; ii < HFElecindex ; ii++)
{
hfePt[ii]=hfelec_Pt[ii];

hfeEta[ii]=hfelec_Eta[ii];
hfePhi[ii]=hfelec_Phi[ii];
hfeL9[ii]= hfelec_L9[ii];
hfeS9[ii]= hfelec_S9[ii];
hfeL25[ii]=  hfelec_L25[ii];
hfeL1[ii]=  hfelec_L1[ii];
h_hfeL1->Fill(hfeL1[ii]);
h_hfeL9->Fill(hfeL9[ii]);
h_hfeS9->Fill(hfeS9[ii]);
h_hfeL25->Fill(hfeL25[ii]);
h_hfecut1->Fill(hfeL9[ii]/hfeL25[ii]);
double cut2 = (hfeL1[ii]/hfeL9[ii]) -(1.125 * (hfeS9[ii]/hfeL9[ii]));
h_hfecut2->Fill(cut2);
if (
(hfeL9[ii]/hfeL25[ii]>0.96)&& (cut2>0.4)&&hfePt[ii]>20.) {
h_hfePt->Fill(hfePt[ii],MyWeight);
h_hfeEta->Fill(hfeEta[ii],MyWeight);

elec_ind_hf[nelec_hf_fire_cut]=ii;
nelec_hf_fire_cut++;
}}


    float p1dotp2[2] ={};
    float dielecpt[2] ={};
    float dielecrapidity[2] ={};
    float dielecphi[2] ={};
    float MZelec[2] ={};

	int select_mva=0;
	int select_cut=0;
	int index1[2] ={};
	int index2[2] ={};
	int iii=-99;
	int mva_cut=0;

////////////////////////
//MVA BASED ELECTRONS //
////////////////////////

if(nelec_fire_mva==2){
	int iiii=-99;
			if(eMVATrigId[index1[mva_cut]]<0.5 && eMVATrigId[index2[mva_cut]]<0.5) iiii =0;
			if(eMVATrigId[index1[mva_cut]]>0.5 && eMVATrigId[index2[mva_cut]]>0.5) iiii =1;
			if((eMVATrigId[index1[mva_cut]]<0.5 && eMVATrigId[index2[mva_cut]]>0.5)||(eMVATrigId[index1[mva_cut]]>0.5 && eMVATrigId[index2[mva_cut]]<0.5))iii=2;
			if(eMVATrigId[index1[mva_cut]]>0.5 && eMVATrigId[index2[mva_cut]]<0.5)iiii=4;
			if(eMVATrigId[index1[mva_cut]]<0.5 && eMVATrigId[index2[mva_cut]]>0.5)iiii=5;

int j=elec_ind_mva[0], jk = elec_ind_mva[1];

	       p1dotp2[0] = ePx[j]*ePx[jk]+ePy[j]*ePy[jk]+ePz[j]*ePz[jk];
               MZelec[0] = fabs(sqrt(pow(eM[j],2)+pow(eM[jk],2)+2*(eE[j]*eE[jk]-p1dotp2[0])));
//		float Mz = fabs(sqrt(pow(eE[j]+eE[jk],2)- (pow(ePx[j]+ePx[jk],2)+pow(ePy[j]+ePy[jk],2)+pow(ePz[j]+ePz[jk],2) )));


if (eCharge[j]*eCharge[jk] != -1){
h_mz_same_sign[iiii]->Fill(MZelec[0],MyWeight); h_mz_same_sign[3]->Fill(MZelec[0],MyWeight);
}

if (eCharge[j]*eCharge[jk] == -1) {

h_no_good_vtx_zevents[iiii]-> Fill(nGoodVertices,MyWeight); h_no_good_vtx_zevents[3]-> Fill(nGoodVertices,MyWeight);
h_no_good_vtx_noWeight_zevents[iiii]-> Fill(nGoodVertices); h_no_good_vtx_noWeight_zevents[3]-> Fill(nGoodVertices);

h_elec_mva_l_mva->Fill(eMVATrigId[j],MyWeight);
h_elec_mva_nl_mva->Fill(eMVATrigId[jk],MyWeight);


	
                  dielecpt[0] = sqrt(pow((ePx[j]+ePx[jk]),2) +pow((ePy[j]+ePy[jk]),2));

                  dielecrapidity[0] = 0.5*log((eE[j]+eE[jk]+ePz[j]+ePz[jk])/(eE[j]+eE[jk]-ePz[j]-ePz[jk]));
                  dielecphi[0] = atan2(ePy[j]+ePy[jk],ePx[j]+ePx[jk]);

	h_elec_pt[iiii]->Fill(ePt[j],MyWeight); h_elec_pt[3]->Fill(ePt[j],MyWeight);
	h_elec_eta[iiii]->Fill(eEta[j],MyWeight); h_elec_eta[3]->Fill(eEta[j],MyWeight);



      h_mz[iiii]->Fill(MZelec[0],MyWeight); h_mz[3]->Fill(MZelec[0],MyWeight);
      h_dielecphi[iiii]->Fill(dielecphi[0],MyWeight); h_dielecphi[3]->Fill(dielecphi[0],MyWeight);
      h_dielec_PT[iiii]->Fill(dielecpt[0],MyWeight); h_dielec_PT[3]->Fill(dielecpt[0],MyWeight);
      h_dielec_rapidity[iiii]->Fill(dielecrapidity[0],MyWeight); h_dielec_rapidity[3]->Fill(dielecrapidity[0],MyWeight);

		index1[0]=j;
		index2[0]=jk;
		select_mva=1;	    
}
}
////////////////////////
//CUT BASED ELECTRONS //
////////////////////////
if(nelec_fire_cut==2){

//cout<<"test"<<endl;
int j=elec_ind_cut[0], jk = elec_ind_cut[1];
	       p1dotp2[1] = ePx[j]*ePx[jk]+ePy[j]*ePy[jk]+ePz[j]*ePz[jk];
               MZelec[1] = fabs(sqrt(pow(eM[j],2)+pow(eM[jk],2)+2*(eE[j]*eE[jk]-p1dotp2[1])));
//		float Mz = fabs(sqrt(pow(eE[j]+eE[jk],2)- (pow(ePx[j]+ePx[jk],2)+pow(ePy[j]+ePy[jk],2)+pow(ePz[j]+ePz[jk],2) )));

if (eCharge[j]*eCharge[jk] != -1){
h_mz_same_sign[9]->Fill(MZelec[1],MyWeight); 
}

if (eCharge[j]*eCharge[jk] == -1) {

h_no_good_vtx_zevents[6]-> Fill(nGoodVertices,MyWeight); 
h_no_good_vtx_noWeight_zevents[6]-> Fill(nGoodVertices);

h_elec_mva_l_cut->Fill(eMVATrigId[j],MyWeight);
h_elec_mva_nl_cut->Fill(eMVATrigId[jk],MyWeight);


	
                  dielecpt[1] = sqrt(pow((ePx[j]+ePx[jk]),2) +pow((ePy[j]+ePy[jk]),2));

                  dielecrapidity[1] = 0.5*log((eE[j]+eE[jk]+ePz[j]+ePz[jk])/(eE[j]+eE[jk]-ePz[j]-ePz[jk]));
                  dielecphi[1] = atan2(ePy[j]+ePy[jk],ePx[j]+ePx[jk]);

	h_elec_pt[6]->Fill(ePt[j],MyWeight); 
	h_elec_eta[6]->Fill(eEta[j],MyWeight);



      h_mz[6]->Fill(MZelec[1],MyWeight);
      h_dielecphi[6]->Fill(dielecphi[1],MyWeight); 
      h_dielec_PT[6]->Fill(dielecpt[1],MyWeight); 
      h_dielec_rapidity[6]->Fill(dielecrapidity[1],MyWeight); 

		index1[1]=j;
		index2[1]=jk;
		select_cut=1;	    
}
}
if (select_mva==0 && select_cut==0) continue;
int good_mva_evt=0;
int good_cut_evt=0;
 // if (MZelec[0] < 60 || MZelec[0] > 120 ) continue;
  if (select_mva==1 &&MZelec[0] > 60 && MZelec[0] < 120 ) good_mva_evt=1;
  if (select_cut==1 &&MZelec[1] > 60 && MZelec[1] < 120 ) good_cut_evt=1;
if(good_mva_evt==0 && good_cut_evt==0) continue;
if(good_cut_evt==1 && good_mva_evt==0){iii=6; mva_cut=1;}

if(good_mva_evt==1 && good_cut_evt==0){
			mva_cut=0;
			if(eMVATrigId[index1[mva_cut]]<0.5 && eMVATrigId[index2[mva_cut]]<0.5) iii =0;
			if(eMVATrigId[index1[mva_cut]]>0.5 && eMVATrigId[index2[mva_cut]]>0.5) iii =1;
			if((eMVATrigId[index1[mva_cut]]<0.5 && eMVATrigId[index2[mva_cut]]>0.5)||(eMVATrigId[index1[mva_cut]]>0.5 && eMVATrigId[index2[mva_cut]]<0.5))iii=2;
			if(eMVATrigId[index1[mva_cut]]>0.5 && eMVATrigId[index2[mva_cut]]<0.5)iii=4;
			if(eMVATrigId[index1[mva_cut]]<0.5 && eMVATrigId[index2[mva_cut]]>0.5)iii=5;

}

if(good_cut_evt==1 && good_mva_evt==1){
if(elec_ind_cut[0]!=elec_ind_mva[0] || elec_ind_cut[1]!=elec_ind_mva[1] )continue;
iii=7; mva_cut=0;
}

      double deltaphi = 999.0;
      int j_in=0;
      int markPFjet;
      int markPFjet_hf;
           int numberofPFjets = 0;
           int jetindex[30] ={};
	   double Long[50]={};
	   double Short[50]={};
           int numberofPFjets_hbhe = 0;
           int numberofPFjets_hf = 0;

      for (int PF=0;PF < pfNjets ;PF++){
	if (jetPt[PF] < 30.) continue;
           markPFjet = 0;
	   markPFjet_hf = 0;
if(fabs(jetPt[PF])>2.853){
		for(int ii = 0 ; ii < nelec_hf_fire_cut ; ii++){
        float deltar_hf = fabs(DeltaR(jetEta[PF],hfeEta[elec_ind_hf[ii]],jetPhi[PF],hfePhi[elec_ind_hf[ii]]));
	h_deltar_elec_PFjet_hf[iii] -> Fill(deltar_hf,MyWeight); 
if(good_mva_evt==1 && good_cut_evt==0){
h_deltar_elec_PFjet_hf[3] -> Fill(deltar_hf,MyWeight);
}
	        if (deltar_hf<0.5) markPFjet_hf = 1; 
		}
}

	  float deltar = fabs(DeltaR(jetEta[PF],eEta[index1[mva_cut]],jetPhi[PF],ePhi[index1[mva_cut]]));
	  float deltar2 = fabs(DeltaR(jetEta[PF],eEta[index2[mva_cut]],jetPhi[PF],ePhi[index2[mva_cut]]));

	h_deltar_elec_PFjet[iii] -> Fill(deltar,MyWeight); 
	h_deltar_elec_PFjet[iii] -> Fill(deltar2,MyWeight); 

if(good_mva_evt==1 && good_cut_evt==0){
h_deltar_elec_PFjet[3] -> Fill(deltar,MyWeight);
h_deltar_elec_PFjet[3] -> Fill(deltar2,MyWeight);
}

	  if (deltar<0.5) markPFjet = 1; 
	  if (deltar2<0.5) markPFjet = 1; 
	  if ((markPFjet != 1&& fabs(jetEta[PF])<2.853)|| (markPFjet_hf != 1 && fabs(jetEta[PF])>2.853)){
	  numberofPFjets += 1;
	  jetindex[j_in]=PF;
                  if (fabs(jetEta[PF])<2.853){
	  			numberofPFjets_hbhe += 1;
         				     }

	                if(fabs(jetEta[PF])>2.853){
			        numberofPFjets_hf += 1;
				Long[j_in]= 0.5 * HadEHF[PF] + EmEHF[PF];
				Short[j_in]= 0.5 * HadEHF[PF];
				h_short->Fill(Short[j_in]);
				h_long->Fill(Long[j_in]);
						  }
j_in++;
}			
}

int jetindex1= jetindex[0];
deltaphi = fabs(DeltaPhi(dielecphi[mva_cut],jetPhi[jetindex1]));
h_numberofPFjets[iii]->Fill(numberofPFjets,MyWeight); 
h_numberofPFjets_hf[iii]->Fill(numberofPFjets_hf,MyWeight); 
h_numberofPFjets_hbhe[iii]->Fill(numberofPFjets_hbhe,MyWeight);
if(numberofPFjets==1)h_deltaphi_z_jet[iii]->Fill(deltaphi,MyWeight); 

if(good_mva_evt==1 && good_cut_evt==0){
h_numberofPFjets[3]->Fill(numberofPFjets,MyWeight);
h_numberofPFjets_hf[3]->Fill(numberofPFjets_hf,MyWeight);
h_numberofPFjets_hbhe[3]->Fill(numberofPFjets_hbhe,MyWeight);
if(numberofPFjets==1)h_deltaphi_z_jet[3]->Fill(deltaphi,MyWeight);
}


	h_leading_jet_pt[iii]->Fill(jetPt[jetindex1],MyWeight); 
	h_leading_jet_eta[iii]->Fill(jetEta[jetindex1],MyWeight);
if(good_mva_evt==1 && good_cut_evt==0){
h_leading_jet_pt[3]->Fill(jetPt[jetindex1],MyWeight);
 h_leading_jet_eta[3]->Fill(jetEta[jetindex1],MyWeight);
}

if(numberofPFjets==1){

   h_jet_pt_all[iii]->Fill(jetPt[jetindex1]);   
   h_jet_eta_all[iii]->Fill(jetEta[jetindex1]);   
   h_mz_1j_all[iii]->Fill(MZelec[mva_cut],MyWeight); 
   h_dielecphi_1j_all[iii]->Fill(dielecphi[mva_cut],MyWeight); 
   h_dielec_PT_1j_all[iii]->Fill(dielecpt[mva_cut],MyWeight); 
   h_dielec_rapidity_1j_all[iii]->Fill(dielecrapidity[mva_cut],MyWeight); 
if(good_mva_evt==1 && good_cut_evt==0){
h_jet_pt_all[3]->Fill(jetPt[jetindex1]);
h_jet_eta_all[3]->Fill(jetEta[jetindex1]);
h_mz_1j_all[3]->Fill(MZelec[mva_cut],MyWeight);
h_dielecphi_1j_all[3]->Fill(dielecphi[mva_cut],MyWeight);
h_dielec_PT_1j_all[3]->Fill(dielecpt[mva_cut],MyWeight);
h_dielec_rapidity_1j_all[3]->Fill(dielecrapidity[mva_cut],MyWeight);
}
if (fabs(jetEta[jetindex1])<2.853){
if(numberofPFjets_hbhe!=1) cout<<"DANGER THERE IS AN ERROR 1"<<endl;
  h_jet_pt_hbhe[iii]->Fill(jetPt[jetindex1]);   
   h_jet_eta_hbhe[iii]->Fill(jetEta[jetindex1]);   
  h_mz_1j_c[iii]->Fill(MZelec[mva_cut],MyWeight); 
  h_dielecphi_1j_c[iii]->Fill(dielecphi[mva_cut],MyWeight); 
  h_dielec_PT_1j_c[iii]->Fill(dielecpt[mva_cut],MyWeight); 
  h_dielec_rapidity_1j_c[iii]->Fill(dielecrapidity[mva_cut],MyWeight); 
if(good_mva_evt==1 && good_cut_evt==0){
h_jet_pt_hbhe[3]->Fill(jetPt[jetindex1]);
h_jet_eta_hbhe[3]->Fill(jetEta[jetindex1]);
h_mz_1j_c[3]->Fill(MZelec[mva_cut],MyWeight);
h_dielecphi_1j_c[3]->Fill(dielecphi[mva_cut],MyWeight);
h_dielec_PT_1j_c[3]->Fill(dielecpt[mva_cut],MyWeight);
h_dielec_rapidity_1j_c[3]->Fill(dielecrapidity[mva_cut],MyWeight);
}
}

if (fabs(jetEta[jetindex1])>2.853){
if(numberofPFjets_hf!=1) cout<<"DANGER THERE IS AN ERROR 2"<<endl;
  h_jet_pt_hf[iii]->Fill(jetPt[jetindex1]);   
  h_jet_eta_hf[iii]->Fill(jetEta[jetindex1]);   
  h_mz_1j_hf[iii]->Fill(MZelec[mva_cut],MyWeight); 
  h_dielecphi_1j_hf[iii]->Fill(dielecphi[mva_cut],MyWeight); 
  h_dielec_PT_1j_hf[iii]->Fill(dielecpt[mva_cut],MyWeight); 
  h_dielec_rapidity_1j_hf[iii]->Fill(dielecrapidity[mva_cut],MyWeight); 
if(good_mva_evt==1 && good_cut_evt==0){
h_jet_eta_hf[3]->Fill(jetEta[jetindex1]);
h_jet_pt_hf[3]->Fill(jetPt[jetindex1]);
h_dielecphi_1j_hf[3]->Fill(dielecphi[mva_cut],MyWeight);
h_dielec_PT_1j_hf[3]->Fill(dielecpt[mva_cut],MyWeight);
h_dielec_rapidity_1j_hf[3]->Fill(dielecrapidity[mva_cut],MyWeight);
}
}

}





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
