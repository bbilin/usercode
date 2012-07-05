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


myTree.Add("/tmp/bbilin/mc_09_04_2012.root");
TFile *theFile = new TFile("17_04_mc.root","RECREATE");


TH1::AddDirectory(true);

  int event,run,lumi,bxnumber,realdata;
//  int hlt_trigger_fired;
  int Elecindex, HFElecindex;
  int ElecIsEE[50],ElecIsEB[50];
  float ElecPt[50], ElecEta[50], ElecPhi[50],ElecPx[50], ElecPy[50], ElecPz[50], ElecM[50], ElecE[50], Elecdr03TkSumPt[50], Elecdr03EcalRecHitSumEt[50], Elecdr03HcalTowerSumEt[50], ElecscSigmaIEtaIEta[50], ElecdeltaPhiSuperClusterTrackAtVtx[50], ElecdeltaEtaSuperClusterTrackAtVtx[50], ElechadronicOverEm[50],ElecgsfTrack_numberOfLostHits[50],ElecGsfTrk_d0[50],ElecDcotTheta[50],hfelec_Pt[50],hfelec_Px[50],hfelec_Py[50],hfelec_Pz[50],hfelec_Phi[50],hfelec_E[50],hfelec_M[50],hfelec_Eta[50], hfelec_L9[50], hfelec_S9[50], hfelec_L25[50], hfelec_L1[50], hhfclustereta[50], hhfclusterphi[50],hfelec_charge[50],ElecMVATrigId[50],ElecMVANonTrigId[50];

double MyWeight;

  //particle information
  int par_index, mom[50], daug[50];
  float ParticlePt[50], ParticleEta[50], ParticlePhi[50], ParticlePx[50], ParticlePy[50], ParticlePz[50], ParticleE[50], ParticleM[50];
  int ParticleId[50], ParticleStatus[50], ParticleMother[50][10], ParticleDaughter[50][10];
  //  int id_elec[20];
  int ElecCharge[50];
  //reco jets 
  int pfNjets;

  int techTrigger[50],nGoodVertices;

 
  float PFjetEta[50], PFjetPhi[50],PFCorrjetPt[50],PFjetPx[50],PFjetPy[50],PFjetPz[50],PFjetE[50],PFHadEHF[100],PFEmEHF[100];
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
  myTree.SetBranchAddress("ElecMVATrigId",ElecMVATrigId);
  myTree.SetBranchAddress("ElecMVANonTrigId",ElecMVANonTrigId);


myTree.SetBranchAddress("HFElecindex",&HFElecindex );
myTree.SetBranchAddress("hfelec_Pt",hfelec_Pt );
myTree.SetBranchAddress("hfelec_Px",hfelec_Px );
myTree.SetBranchAddress("hfelec_Py",hfelec_Py );
myTree.SetBranchAddress("hfelec_Pz",hfelec_Pz );
myTree.SetBranchAddress("hfelec_E",hfelec_E );
myTree.SetBranchAddress("hfelec_M",hfelec_M );
myTree.SetBranchAddress("hfelec_Eta",hfelec_Eta );
myTree.SetBranchAddress("hfelec_Phi",hfelec_Phi );
myTree.SetBranchAddress("hfelec_L9", hfelec_L9);
myTree.SetBranchAddress("hfelec_S9",hfelec_S9 );
myTree.SetBranchAddress("hfelec_L25", hfelec_L25);
myTree.SetBranchAddress("hfelec_L1",hfelec_L1 );
myTree.SetBranchAddress("hhfclustereta",hhfclustereta );
myTree.SetBranchAddress("hhfclusterphi",hhfclusterphi );
myTree.SetBranchAddress("hfelec_charge",hfelec_charge );

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
 
myTree.SetBranchAddress("MyWeight", &MyWeight );
  
  //recojets

 



  myTree.SetBranchAddress("pfNjets",&pfNjets);
  myTree.SetBranchAddress("PFjetEta",PFjetEta);
  myTree.SetBranchAddress("PFjetPhi",PFjetPhi);
  myTree.SetBranchAddress("PFCorrjetPt",PFCorrjetPt);
  myTree.SetBranchAddress("PFjetPx",PFjetPx);
  myTree.SetBranchAddress("PFjetPy",PFjetPy);
  myTree.SetBranchAddress("PFjetPz",PFjetPz);
  myTree.SetBranchAddress("PFjetE",PFjetE);
  myTree.SetBranchAddress("PFHadEHF",PFHadEHF);
  myTree.SetBranchAddress("PFEmEHF",PFEmEHF);



  myTree.SetBranchAddress("nVertices",&nVertices);

  theFile->cd();
TH1F* h_no_good_vtx = new TH1F ("h_no_good_vtx","h_no_good_vtx",20,-0.5,19.5);

TH1F* h_no_good_vtx_noWeight = new TH1F ("h_no_good_vtx_noWeight","h_no_good_vtx_noWeight",20,-0.5,19.5);

TH1F *h_deltaphi_z_jet_c=  new TH1F ("h_deltaphi_z_jet_c","h_deltaphi_z_jet_c",35,0,3.5);

//TH1F *h_deltaphi_z_jet_1hf=  new TH1F ("h_deltaphi_z_jet_1hf","h_deltaphi_z_jet_1hf",35,0,3.5);

//TH1F *h_deltaphi_z_jet_2hf=  new TH1F ("h_deltaphi_z_jet_2hf","h_deltaphi_z_jet_2hf",35,0,3.5);

TH1F *h_elec_pt=  new TH1F ("h_elec_pt","h_elec_pt",70,0,200);
TH1F *h_elec_eta=  new TH1F ("h_elec_eta","h_elec_eta",80,-5,5);

TH1F *h_ept_hf=  new TH1F ("h_ept_hf","h_ept_hf",70,0,200);
TH1F *h_eeta_hf=  new TH1F ("h_eeta_hf","h_eeta_hf",80,-6,6);

TH1F *h_hfePt=  new TH1F ("h_hfe_Pt","h_hfe_Pt",70,0,200);
TH1F *h_hfeEta=  new TH1F ("h_hfe_Eta","h_hfe_Eta",80,-6,6);



TH1F *h_leading_jet_pt_c=  new TH1F ("h_leading_jet_pt_c","h_leading_jet_pt_c",80,0,200);
TH1F *h_leading_jet_eta_c=  new TH1F ("h_leading_jet_eta_c","h_leading_jet_eta_c",40,-6,6);

//TH1F *h_leading_jet_pt_1hf=  new TH1F ("h_leading_jet_pt_1hf","h_leading_jet_pt_1hf",80,0,200);
//TH1F *h_leading_jet_eta_1hf=  new TH1F ("h_leading_jet_eta_1hf","h_leading_jet_eta_1hf",40,-6,6);

//TH1F *h_leading_jet_pt_2hf=  new TH1F ("h_leading_jet_pt_2hf","h_leading_jet_pt_2hf",80,0,200);
//TH1F *h_leading_jet_eta_2hf=  new TH1F ("h_leading_jet_eta_2hf","h_leading_jet_eta_2hf",40,-6,6);



TH1F *h_mz_c=  new TH1F ("h_mz_c","h_mz_c",80,0,200);
TH1F *h_mz_1hf=  new TH1F ("h_mz_1hf","h_mz_1hf",50,0,200);

TH1F *h_mz_hf_e1=  new TH1F ("h_mz_hf_e1","h_mz_hf_e1",50,0,200 );
TH1F *h_mz_hf_e2=  new TH1F ("h_mz_hf_e2","h_mz_hf_e2",50,0,200 );
TH1F *h_mz_hf_e3=  new TH1F ("h_mz_hf_e3","h_mz_hf_e3",50,0,200 );
//TH1F *h_mz_2hf=  new TH1F ("h_mz_2hf","h_mz_2hf",50,0,200);


TH1F *h_dielecphi_c=  new TH1F ( "h_dielecphi_c", "h_dielecphi_c",20,-3,3);
TH1F *h_dielec_PT_c=  new TH1F ( "h_dielec_PT_c", "h_dielec_PT_c",100,0,200);
TH1F *h_dielec_rapidity_c=  new TH1F ( "h_dielec_rapidity_c", "h_dielec_rapidity_c",50,-5,5);

TH1F *h_dielecphi_1hf=  new TH1F ( "h_dielecphi_1hf", "h_dielecphi_1hf",20,-3,3);
TH1F *h_dielec_PT_1hf=  new TH1F ( "h_dielec_PT_1hf", "h_dielec_PT_1hf",100,0,200);
TH1F *h_dielec_rapidity_1hf=  new TH1F ( "h_dielec_rapidity_1hf", "h_dielec_rapidity_1hf",50,-5,5);

//TH1F *h_dielecphi_2hf=  new TH1F ( "h_dielecphi_2hf", "h_dielecphi_2hf",20,-3,3);
//TH1F *h_dielec_PT_2hf=  new TH1F ( "h_dielec_PT_2hf", "h_dielec_PT_2hf",100,0,200);
//TH1F *h_dielec_rapidity_2hf=  new TH1F ( "h_dielec_rapidity_2hf", "h_dielec_rapidity_2hf",50,-5,5);


TH1F *h_deltar_elec_PFjet_c=  new TH1F ("h_deltar_elec_PFjet_c","h_deltar_elec_PFjet_c",100,0,2);
//TH1F *h_deltar_elec_PFjet_1hf=  new TH1F ("h_deltar_elec_PFjet_1hf","h_deltar_elec_PFjet_1hf",100,0,2);
//TH1F *h_deltar_elec_PFjet_2hf=  new TH1F ("h_deltar_elec_PFjet_2hf","h_deltar_elec_PFjet_2hf",100,0,2);

TH1F *h_numberofPFjets = new TH1F ("h_numberofPFjets","h_numberofPFjets",14,-0.5,13.5);

TH2F *h_elec_events = new TH2F ("h_elec_events","h_elec_events",4,0,4,4,0,4);

 Int_t nevent = myTree.GetEntries();
 //Int_t  nevent = 4000000;
  for (Int_t iev=0;iev<nevent;iev++) {

//if(iev<8)continue;


 	  if (iev%100000 == 0) cout<<iev<<"/"<<nevent<<endl;

    	myTree.GetEntry(iev); 

//if(realdata && tr_1!=1 && tr_2!=1 && tr_3!=1) continue;
if (realdata) MyWeight=1;
//cout<<"ahmetmehmet"<<MyWeight<<endl;
h_no_good_vtx-> Fill(nGoodVertices,MyWeight);

h_no_good_vtx_noWeight-> Fill(nGoodVertices);

//cout<<"testtttttt"<<endl;
int elec_ind[50], elec_ind_l[50],elec_ind_hf[50];
	for(int i=0; i<50; i++)
			{elec_ind[i]=0;
			elec_ind_hf[i]=0;
			elec_ind_l[i]=0;}

int eCharge[50]={}, eIsEE[50]={}, eIsEB[50]= {},e_cut_fired[50]={},nelec_fire_cut=0, e_cut_fired_l[50]={},nelec_fire_cut_l=0,nelec_hf_fire_cut=0;

double ePt[50] = {}, eEta[50] = {}, ePhi[50]= {}, ePx[50]={}, ePy[50]={}, ePz[50] = {}, eM[50]= {}, eE[50]={}, edr03TkSumPt[50]={}, edr03EcalRecHitSumEt[50]={}, edr03HcalTowerSumEt[50]= {}, escSigmaIEtaIEta[50]={}, edeltaPhiSuperClusterTrackAtVtx[50]={}, edeltaEtaSuperClusterTrackAtVtx[50]= {}, ehadronicOverEm[50]={} , egsfTrack_numberOfLostHits[50]={} , eGsfTrk_d0[50]={} , eDcotTheta[50]={}; 

double jetEta[50]={}, jetPhi[50]={}, jetPt[50]={}, jetPx[50]={}, jetPy[50]={}, jetPz[50]={}, jetE[50]={}, EmEHF[50]={}, HadEHF[50]={},jPt[50]={};

int sorted_jet_index[50]={};



for(int ii=0;ii<pfNjets;ii++){
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

for(int ii=0;ii<pfNjets;ii++){
jPt[ii]=fabs(PFCorrjetPt[ii]);
}

//double corr_fact=-999.;
int sort_jet = 0;
int ind=0;
TMath::Sort(pfNjets,jPt,sorted_jet_index);

for(int ii = 0 ; ii < pfNjets ; ii++)
{

 sort_jet = sorted_jet_index[ind];

jetEta[ind]=PFjetEta[sort_jet]  ;
jetPhi[ind]=PFjetPhi[sort_jet] ;
jetPt[ind]= jPt[sort_jet];
if(jetPt[ind]<-10)cout<<"  "<<jetPt[0]<<"  "<<jetPt[1]<<"  "<<jetPt[ind]<<endl;
jetPx[ind]= PFjetPx[sort_jet];
jetPy[ind]= PFjetPy[sort_jet];
jetPz[ind]= PFjetPz[sort_jet];
jetE[ind]= PFjetE[sort_jet];
HadEHF[ind]=PFHadEHF[sort_jet];
EmEHF[ind]=PFEmEHF[sort_jet];
ind++;
}


//cout<<"new"<<endl;
for(int ii = 0 ; ii < Elecindex ; ii++){

eCharge[ii]=ElecCharge[ii];
eIsEE[ii]=ElecIsEE[ii];
eIsEB[ii]= ElecIsEB[ii];

eEta[ii] = ElecEta[ii];
ePhi[ii]= ElecPhi[ii];


ePt[ii] = ElecPt[ii];
ePx[ii]= ElecPx[ii];
ePy[ii]=  ElecPy[ii];
ePz[ii] =ElecPz[ii];


if ((fabs(eEta[ii]))>2.853)  {
//cout<<"existence   "<<eCharge[ii]<<endl;
h_ept_hf->Fill(ePt[ii]);
h_eeta_hf->Fill(eEta[ii]);
} 

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

}

    float p1dotp2 = -99.;
    float dielecpt = -999.;
    float dielecrapidity = -999.;
    float dielecphi = -999.;
    float MZelec  = -999. ;

	int select=0;
	int index1 =0;
	int index2 =0;

if(nelec_fire_cut==2){
//cout<<"test"<<endl;
int j=elec_ind[0], jk = elec_ind[1];
if (eCharge[j]*eCharge[jk] != -1) continue;
	       p1dotp2 = ePx[j]*ePx[jk]+ePy[j]*ePy[jk]+ePz[j]*ePz[jk];
               MZelec = fabs(sqrt(eM[j]+eM[jk]+2*(eE[j]*eE[jk]-p1dotp2)));
	

                  dielecpt = sqrt(pow((ePx[j]+ePx[jk]),2) +pow((ePy[j]+ePy[jk]),2));

                  dielecrapidity = 0.5*log((eE[j]+eE[jk]+ePz[j]+ePz[jk])/(eE[j]+eE[jk]-ePz[j]-ePz[jk]));
                  dielecphi = atan2(ePy[j]+ePy[jk],ePx[j]+ePx[jk]);
index1=j;
index2=jk;
	h_elec_pt->Fill(ePt[j],MyWeight);
//	h_elec_pt->Fill(ePt[jk]);
	h_elec_eta->Fill(eEta[j],MyWeight);
//	h_elec_eta->Fill(eEta[jk]);

      h_mz_c->Fill(MZelec,MyWeight);
      h_dielecphi_c->Fill(dielecphi,MyWeight);
      h_dielec_PT_c->Fill(dielecpt,MyWeight);
      h_dielec_rapidity_c->Fill(dielecrapidity,MyWeight);

		select=1;	    
}


if(nelec_hf_fire_cut==1 &&nelec_fire_cut==1){

int j=elec_ind[0], jk = elec_ind_hf[0];

	       p1dotp2 = ePx[j]*hfePx[jk]+ePy[j]*hfePy[jk]+ePz[j]*hfePz[jk];
               MZelec = fabs(sqrt(eM[j]+hfeM[jk]+2*(eE[j]*hfeE[jk]-p1dotp2)));

                  dielecpt = sqrt(pow((ePx[j]+hfePx[jk]),2) +pow((ePy[j]+hfePy[jk]),2));

                  dielecrapidity = 0.5*log((eE[j]+hfeE[jk]+ePz[j]+hfePz[jk])/(eE[j]+hfeE[jk]-ePz[j]-hfePz[jk]));
                  dielecphi = atan2(ePy[j]+hfePy[jk],ePx[j]+hfePx[jk]);

      h_mz_1hf->Fill(MZelec,MyWeight);
      h_dielecphi_1hf->Fill(dielecphi,MyWeight);
      h_dielec_PT_1hf->Fill(dielecpt,MyWeight);
      h_dielec_rapidity_1hf->Fill(dielecrapidity,MyWeight);

index1=j;
index2=jk;

//		select=1;	    
}

float mztest1h=0,mztest2h=0,mztest3h=0, mztest12=0, mztest13=0, mztest23=0, p1dotp2_1h=0,p1dotp2_2h=0,p1dotp2_3h=0, p1dotp2_12=0, p1dotp2_13=0, p1dotp2_23=0;

if(nelec_fire_cut>0){ 
int j=elec_ind[0], jk = elec_ind[1] ,jkl=elec_ind[2] , hf= elec_ind_hf[0] ;

		p1dotp2_12 = ePx[j]*ePx[jk]+ePy[j]*ePy[jk]+ePz[j]*ePz[jk];
		mztest12 = fabs(sqrt(eM[j]+eM[jk]+2*(eE[j]*eE[jk]-p1dotp2_12)));

	       p1dotp2_13 = ePx[j]*ePx[jkl]+ePy[j]*ePy[jkl]+ePz[j]*ePz[jkl];
               mztest13 = fabs(sqrt(eM[j]+eM[jkl]+2*(eE[j]*eE[jkl]-p1dotp2_13)));

	       p1dotp2_23 = ePx[jk]*ePx[jkl]+ePy[jk]*ePy[jkl]+ePz[jk]*ePz[jkl];
               mztest23 = fabs(sqrt(eM[jk]+eM[jkl]+2*(eE[jk]*eE[jkl]-p1dotp2_23)));

	       p1dotp2_1h = ePx[j]*hfePx[hf]+ePy[j]*hfePy[hf]+ePz[j]*hfePz[hf];
               mztest1h = fabs(sqrt(eM[j]+hfeM[hf]+2*(eE[j]*hfeE[hf]-p1dotp2_1h)));

	       p1dotp2_2h = ePx[jk]*hfePx[hf]+ePy[jk]*hfePy[hf]+ePz[jk]*hfePz[hf];
               mztest2h = fabs(sqrt(eM[jk]+hfeM[hf]+2*(eE[jk]*hfeE[hf]-p1dotp2_2h)));

	       p1dotp2_3h = ePx[jkl]*hfePx[hf]+ePy[jkl]*hfePy[hf]+ePz[jkl]*hfePz[hf];
               mztest3h = fabs(sqrt(eM[jkl]+hfeM[hf]+2*(eE[jkl]*hfeE[hf]-p1dotp2_3h)));

if(nelec_hf_fire_cut==1){
if(nelec_fire_cut>0)h_mz_hf_e1->Fill(mztest1h);
if(nelec_fire_cut>1)h_mz_hf_e2->Fill(mztest2h);
if(nelec_fire_cut>2)h_mz_hf_e3->Fill(mztest3h);
}

if(nelec_hf_fire_cut==1 && mztest1h>80. && mztest1h<100.)h_elec_events->Fill(0.,3.);

if(nelec_fire_cut>1){
if(mztest12>80. && mztest12<100.)h_elec_events->Fill(0.,1.);
if(nelec_hf_fire_cut==1 && mztest2h>80. && mztest2h<100.)h_elec_events->Fill(1.,3.);
}

if(nelec_fire_cut>2){
if(mztest13>80. && mztest13<100.)h_elec_events->Fill(0.,2.);
if(mztest23>80. && mztest23<100.)h_elec_events->Fill(1.,2.);
if(nelec_hf_fire_cut==1 && mztest3h>80. && mztest3h<100.)h_elec_events->Fill(2.,3.);
}
}

if (select==0) continue;

//      if (MZelec < 60 || MZelec > 120) continue;

  if (MZelec < 80 || MZelec > 100/* || dielecpt<20.*/ ) continue;
  
      
            int numberofPFjets = 0;
      int jetindex[30] ={}; 

      double deltaphi = 999.0;

      int j_in=0;
      int markPFjet;
      for (int PF=0;PF < pfNjets ;PF++){
           markPFjet = 0;

	  float deltar = fabs(DeltaR(jetEta[PF],eEta[index1],jetPhi[PF],ePhi[index1]));
	  float deltar2 = fabs(DeltaR(jetEta[PF],eEta[index2],jetPhi[PF],ePhi[index2]));


h_deltar_elec_PFjet_c -> Fill(deltar,MyWeight);
	h_deltar_elec_PFjet_c -> Fill(deltar2,MyWeight);



	  if (deltar<0.5) markPFjet = 1; 
	  if (deltar2<0.5) markPFjet = 1; 
	if (markPFjet != 1){
	if (jetPt[PF] > 30.){
	  numberofPFjets += 1;
	  jetindex[j_in]=PF;

j_in++;
}		
}}
h_numberofPFjets->Fill(numberofPFjets,MyWeight);
if (numberofPFjets!=1)continue;
int jetindex1=jetindex[0];


//cout<< markgoodev<<endl;

deltaphi = fabs(DeltaPhi(dielecphi,jetPhi[jetindex1]));

h_deltaphi_z_jet_c->Fill(deltaphi,MyWeight);
if (deltaphi<2.8) continue;

	h_leading_jet_pt_c->Fill(jetPt[jetindex1],MyWeight);
	h_leading_jet_eta_c->Fill(jetEta[jetindex1],MyWeight);





//cout<<cjetPt[jetindex1]<<endl; 


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
