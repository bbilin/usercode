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
//2012 DATA 

myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/data_2012A.root");
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/data_2012B.root");
TFile *theFile = new TFile("05_07_data.root","RECREATE");


//MC
/*
myTree.Add("root://eoscms.cern.ch//eos/cms/store/user/bbilin/ntuples/mc_Madgraph_tarball.root");
TFile *theFile = new TFile("05_07_mc.root","RECREATE");
*/

TH1::AddDirectory(true);

  int event,run,lumi,bxnumber,realdata;
//  int hlt_trigger_fired;
  int Elecindex;
  float ElecPt[50], ElecEta[50], ElecPhi[50],ElecPx[50], ElecPy[50], ElecPz[50], ElecM[50], ElecE[50],ElecMVATrigId[50],ElecMVANonTrigId[50];

double MyWeight;

  //particle information
  int par_index, mom[50], daug[50];
 
  int ParticleId[50], ParticleStatus[50], ParticleMother[50][10], ParticleDaughter[50][10];
  //  int id_elec[20];
  int ElecCharge[50];
  //reco jets 
  int pfNjets;

  int techTrigger[50],nGoodVertices;

 
  float PFjetEta[50], PFjetPhi[50],PFCorrjetPt[50];
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
 

  myTree.SetBranchAddress("nVertices",&nVertices);

  theFile->cd();

TH1F* h_no_good_vtx = new TH1F ("h_no_good_vtx","h_no_good_vtx",50,-0.5,49.5);
TH1F* h_no_good_vtx_noWeight = new TH1F ("h_no_good_vtx_noWeight","h_no_good_vtx_noWeight",50,-0.5,49.5);

TH1F *h_deltaphi_z_jet=  new TH1F ("h_deltaphi_z_jet","h_deltaphi_z_jet",35,0,3.5);
TH1F *h_elec_mva_l=  new TH1F ("h_elec_mva_l","h_elec_mva_l",30,-1.5,1.5);
TH1F *h_elec_mva_nl=  new TH1F ("h_elec_mva_nl","h_elec_mva_nl",30,-1.5,1.5);

TH1F *h_elec_pt=  new TH1F ("h_elec_pt","h_elec_pt",70,0,200);
TH1F *h_elec_eta=  new TH1F ("h_elec_eta","h_elec_eta",80,-5,5);

TH1F *h_leading_jet_pt=  new TH1F ("h_leading_jet_pt","h_leading_jet_pt",80,0,200);
TH1F *h_leading_jet_eta=  new TH1F ("h_leading_jet_eta","h_leading_jet_eta",40,-6,6);

TH1F *h_mz=  new TH1F ("h_mz","h_mz",80,0,200);

TH1F *h_dielecphi=  new TH1F ( "h_dielecphi", "h_dielecphi",20,-3,3);
TH1F *h_dielec_PT=  new TH1F ( "h_dielec_PT", "h_dielec_PT",100,0,200);
TH1F *h_dielec_rapidity=  new TH1F ( "h_dielec_rapidity", "h_dielec_rapidity",50,-5,5);

      TH1F *h_mz_1j_all=  new TH1F ( "h_mz_1j_all", "h_mz_1j_all",80,0,200);
      TH1F *h_dielecphi_1j_all=  new TH1F ( "h_dielecphi_1j_all", "h_dielecphi_1j_all",20,-3,3);
      TH1F *h_dielec_PT_1j_all=  new TH1F ( "h_dielec_PT_1j_all", "h_dielec_PT_1j_all",100,0,200);
      TH1F *h_dielec_rapidity_1j_all=  new TH1F ( "h_dielec_rapidity_1j_all", "h_dielec_rapidity_1j_all",50,-5,5);

      TH1F *h_mz_1j_c=  new TH1F ( "h_mz_1j_c", "h_mz_1j_c",80,0,200);
      TH1F *h_dielecphi_1j_c=  new TH1F ( "h_dielecphi_1j_c", "h_dielecphi_1j_c",20,-3,3);
      TH1F *h_dielec_PT_1j_c=  new TH1F ( "h_dielec_PT_1j_c", "h_dielec_PT_1j_c",100,0,200);
      TH1F *h_dielec_rapidity_1j_c=  new TH1F ( "h_dielec_rapidity_1j_c", "h_dielec_rapidity_1j_c",50,-5,5);

      TH1F *h_mz_1j_hf=  new TH1F ( "h_mz_1j_hf", "h_mz_1j_hf",80,0,200);
      TH1F *h_dielecphi_1j_hf=  new TH1F ( "h_dielecphi_1j_hf", "h_dielecphi_1j_hf",20,-3,3);
      TH1F *h_dielec_PT_1j_hf=  new TH1F ( "h_dielec_PT_1j_hf", "h_dielec_PT_1j_hf",100,0,200);
      TH1F *h_dielec_rapidity_1j_hf=  new TH1F ( "h_dielec_rapidity_1j_hf", "h_dielec_rapidity_1j_hf",50,-5,5);



TH1F *h_deltar_elec_PFjet=  new TH1F ("h_deltar_elec_PFjet","h_deltar_elec_PFjet",100,0,2);

TH1F *h_numberofPFjets = new TH1F ("h_numberofPFjets","h_numberofPFjets",14,-0.5,13.5);

 Int_t nevent = myTree.GetEntries();
  for (Int_t iev=0;iev<nevent;iev++) {



 	  if (iev%100000 == 0) cout<<iev<<"/"<<nevent<<endl;

    	myTree.GetEntry(iev); 

if(realdata && tr_3!=1) continue;

/*if (realdata)*/ MyWeight=1; //no PU Reweighting NOW!!

h_no_good_vtx-> Fill(nGoodVertices,MyWeight);

h_no_good_vtx_noWeight-> Fill(nGoodVertices);

int elec_ind[50];

for(int i=0; i<50; i++){elec_ind[i]=0;}

int eCharge[50]={},nelec_fire_cut=0;
double ePt[50] = {}, eEta[50] = {}, ePhi[50]= {}, ePx[50]={}, ePy[50]={}, ePz[50] = {}, eM[50]= {}, eE[50]={}, eMVATrigId[50]={}/*,eMVANonTrigId[50]={}*/; 
double jetEta[50]={}, jetPhi[50]={}, jetPt[50]={},jPt[50]={};
int sorted_jet_index[50]={};

for(int ii=0;ii<pfNjets;ii++){
jPt[ii]=0; jetPt[ii]=0; jetEta[ii]=0; jetPhi[ii]=0; jetPt[ii]=0; 
}



    int particle_id[50] = {};

    int nparticle_gen =0;


int mcIsElec=0;
//int mcIsTau=0;

    for (int ppo = 0;ppo<50;++ppo) particle_id[ppo] = 9999; 
    for (int j = 0;j<par_index;++j){
      if (ParticleStatus[j] != 3) continue;
      particle_id[nparticle_gen] = ParticleId[j];
	if (particle_id[nparticle_gen]==11) mcIsElec=1;
//	if (particle_id[nparticle_gen]==15) mcIsTau=1;

nparticle_gen++ ;//9 or 10 only
}

if(!realdata && mcIsElec==0) continue;

for(int ii=0;ii<pfNjets;ii++){
jPt[ii]=fabs(PFCorrjetPt[ii]);
}

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
//eMVANonTrigId[ii]=ElecMVANonTrigId[ii];



if(eMVATrigId[ii]>0 &&ePt[ii] >25.){
elec_ind[nelec_fire_cut]=ii;
nelec_fire_cut++;
}
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

h_elec_mva_l->Fill(eMVATrigId[j]);

h_elec_mva_nl->Fill(eMVATrigId[jk]);

	       p1dotp2 = ePx[j]*ePx[jk]+ePy[j]*ePy[jk]+ePz[j]*ePz[jk];
               MZelec = fabs(sqrt(eM[j]+eM[jk]+2*(eE[j]*eE[jk]-p1dotp2)));
	
                  dielecpt = sqrt(pow((ePx[j]+ePx[jk]),2) +pow((ePy[j]+ePy[jk]),2));

                  dielecrapidity = 0.5*log((eE[j]+eE[jk]+ePz[j]+ePz[jk])/(eE[j]+eE[jk]-ePz[j]-ePz[jk]));
                  dielecphi = atan2(ePy[j]+ePy[jk],ePx[j]+ePx[jk]);
index1=j;
index2=jk;

	h_elec_pt->Fill(ePt[j],MyWeight);
	h_elec_eta->Fill(eEta[j],MyWeight);

      h_mz->Fill(MZelec,MyWeight);
      h_dielecphi->Fill(dielecphi,MyWeight);
      h_dielec_PT->Fill(dielecpt,MyWeight);
      h_dielec_rapidity->Fill(dielecrapidity,MyWeight);

		select=1;	    
}

if (select==0) continue;

  if (MZelec < 80 || MZelec > 100 ) continue;
  
           int numberofPFjets = 0;
           int jetindex[30] ={}; 

      double deltaphi = 999.0;

      int j_in=0;
      int markPFjet;
      for (int PF=0;PF < pfNjets ;PF++){
           markPFjet = 0;

	  float deltar = fabs(DeltaR(jetEta[PF],eEta[index1],jetPhi[PF],ePhi[index1]));
	  float deltar2 = fabs(DeltaR(jetEta[PF],eEta[index2],jetPhi[PF],ePhi[index2]));

	h_deltar_elec_PFjet -> Fill(deltar,MyWeight);
	h_deltar_elec_PFjet -> Fill(deltar2,MyWeight);

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


deltaphi = fabs(DeltaPhi(dielecphi,jetPhi[jetindex1]));

h_deltaphi_z_jet->Fill(deltaphi,MyWeight);

if (deltaphi<2.8) continue;

	h_leading_jet_pt->Fill(jetPt[jetindex1],MyWeight);
	h_leading_jet_eta->Fill(jetEta[jetindex1],MyWeight);

   
if(numberofPFjets==1){

      h_mz_1j_all->Fill(MZelec,MyWeight);
      h_dielecphi_1j_all->Fill(dielecphi,MyWeight);
      h_dielec_PT_1j_all->Fill(dielecpt,MyWeight);
      h_dielec_rapidity_1j_all->Fill(dielecrapidity,MyWeight);

if (fabs(jetEta[jetindex1])<2.8){
    h_mz_1j_c->Fill(MZelec,MyWeight);
      h_dielecphi_1j_c->Fill(dielecphi,MyWeight);
      h_dielec_PT_1j_c->Fill(dielecpt,MyWeight);
      h_dielec_rapidity_1j_c->Fill(dielecrapidity,MyWeight);
}

if (fabs(jetEta[jetindex1])>2.8){
    h_mz_1j_hf->Fill(MZelec,MyWeight);
      h_dielecphi_1j_hf->Fill(dielecphi,MyWeight);
      h_dielec_PT_1j_hf->Fill(dielecpt,MyWeight);
      h_dielec_rapidity_1j_hf->Fill(dielecrapidity,MyWeight);
}

}


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
