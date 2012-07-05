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
#include "TLatex.h"
using namespace std;
void tree1r()
{

  gStyle->SetPalette(1);
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  float DeltaPhi(float phi1, float phi2);

  TChain myTree("analyzeBasicPat/MuonTree");

  
  myTree.Add("/tmp/bbilin/WWtoAnything_V3.root");

  

  TH1::AddDirectory(true);
/*	
	//quark correction factors driven from DATA & MC/////
	TFile *theFile = new TFile ("WW_corr_mc.root","RECREATE");
double par0[10]={320.265,310.384,302.332,294.038,300.43,319.707,299.17,593.3,298.747,300};
double par1[10]={48.6647,32.1916,21.2444,15.866,18.3176,39.7669,20.5297,455.153,-70.0774,20};
double par2[10]={27.2194,28.5505,29.6814,30.8269,29.972,27.2802,30.0906,-19.185,30.0779,30};
double par3[10]={16.7014,14.199,9.92609,8.32889,18.9738,24.479,15.5766,39.5681,50.7993,50};

//double par0[10]={480.572,470.528,470.773,485.117,503.587,563.889,642.667,5499.84,1.84432e+07,402.617,};
//double par1[10]={174.514,152.078,137.116,157.416,215.085,286.139,394.117,9689.66,8.38169e+07,18.8427,};
//double par2[10]={22.0171,23.0618,23.1514,21.7883,19.513,13.6485,5.72768,-606.01,-3.81839e+06,29.7577,};
//double par3[10]={60.1158,55.5369,53.7107,58.638,69.0377,87.6053,114.281,2461.65,2.22927e+07,10.389,};

//double par0[10]={471.741,463.677,475.437,489.5,481.099,547.817,618.151,1576.88,9.97124e+06,500};
//double par1[10]={76.2915,59.4876,56.8093,73.6539,97.5527,169.189,254.539,1699.21,5.2654e+07,30};
//double par2[10]={30.7679,31.5864,31.0351,29.8814,29.7059,24.3482,18.3419,-79.1475,-2.39582e+06,30};
//double par3[10]={79.7822,70.5934,70.0326,80.545,90.9784,132.09,179.83,1029.97,3.07823e+07,50};
*/
double MAX_CALOjetRawPt[10]={};
//DATA///////////////////////////////////////////////////////////////

 TFile *theFile = new TFile ("WW_corr_data.root","RECREATE");
MAX_CALOjetRawPt[0]=  146.096;
MAX_CALOjetRawPt[1]=  148.437;
MAX_CALOjetRawPt[2]=  138.441;
MAX_CALOjetRawPt[3]=  144.366;
MAX_CALOjetRawPt[4]=  151.072;
MAX_CALOjetRawPt[5]=  113.335;
MAX_CALOjetRawPt[6]=  64.931;
MAX_CALOjetRawPt[7]=  39.8758;
MAX_CALOjetRawPt[8]=  0;
MAX_CALOjetRawPt[9]=  0;
/*
double par0[10]={362.232,338.668,330.291,279.676,370.942,324.357,272.108,492.757,500,500,};
double par1[10]={-21.0363,-43.8899,-46.4916,-78.6747,-16.2441,-48.3755,-71.8572,33.4904,30,30,};
double par2[10]={37.9513,39.5706,40.6289,43.2247,37.5116,40.8895,43.2215,30.4557,30,30,};
double par3[10]={24.8149,16.7491,9.5523,7.19913,26.9194,9.1779,6.45521,54.3038,50,50,};

double par0[10]={357.186,405.356,401.468,369.805,361.445,392.61,293.991,384.242,343.989,500,};
double par1[10]={-25.3164,-0.36385,-11.0435,-26.703,-23.4656,-11.4211,-69.0918,-20.4608,-42.7395,30,};
double par2[10]={37.2946,35.4146,36.577,38.5169,37.4113,36.7487,43.1843,38.2994,45.2739,30,};
double par3[10]={26.5531,33.9905,22.5771,16.8856,26.5541,25.3713,10.1463,14.372,4.74744,50,};

double par0[10]={366.827,337.897,327.453,267.843,385.753,320.293,351.524,493.956,500,500,};
double par1[10]={-18.3195,-46.5987,-49.0848,-85.2418,-7.89028,-52.9679,-33.7025,32.9547,30,30,};
double par2[10]={37.735,39.664,40.8675,44.0132,36.7172,41.2938,39.1151,30.3838,30,30,};
double par3[10]={25.7211,16.6791,8.81087,5.96181,31.0787,7.85709,14.3269,53.7737,50,50,};


double par0[10]={362.232,338.668,330.291,279.676,370.942,324.357,272.108,492.757,500,500,};
double par1[10]={-21.0363,-43.8899,-46.4916,-78.6747,-16.2441,-48.3755,-71.8572,33.4904,30,30,};
double par2[10]={37.9513,39.5706,40.6289,43.2247,37.5116,40.8895,43.2215,30.4557,30,30,};
double par3[10]={24.8149,16.7491,9.5523,7.19913,26.9194,9.1779,6.45521,54.3038,50,50,};

double par0[10]={363.027,403.546,406.083,369.461,361.008,404.168,298.632,500,500,500,};
double par1[10]={-21.9681,-1.89841,-7.91237,-24.6082,-23.4216,-5.09109,-73.3075,30,30,30,};
double par2[10]={36.9547,35.482,36.2604,38.583,37.3889,36.0237,42.8583,30,30,30,};
double par3[10]={28.581,33.7846,23.5066,15.7097,26.5694,29.0248,12.329,50,50,50,};
*/
double par0[10]={363.316,403.081,406.134,368.809,364.087,408.667,300.238,500,500,500,};
double par1[10]={-21.7272,-2.18084,-7.89329,-24.9499,-21.2873,-2.24861,-72.9198,30,30,30,};
double par2[10]={36.9416,35.5082,36.2573,38.62,37.2425,35.7492,42.758,30,30,30,};
double par3[10]={28.5793,33.6476,23.5205,15.6399,26.9212,30.0747,12.5792,50,50,50,};

/////////////////////////////////////////////////////////////////////

//MC/////////////////////////////////////////////////////////////////

/* TFile *theFile = new TFile ("WW_corr_mc.root","RECREATE");
MAX_CALOjetRawPt[0]=  154.767;
MAX_CALOjetRawPt[1]=  153.49;
MAX_CALOjetRawPt[2]=  144.698;
MAX_CALOjetRawPt[3]=  147.328;
MAX_CALOjetRawPt[4]=  145.954;
MAX_CALOjetRawPt[5]=  99.5864;
MAX_CALOjetRawPt[6]=  62.1746;
MAX_CALOjetRawPt[7]=  42.64;
MAX_CALOjetRawPt[8]=  20.4477;
MAX_CALOjetRawPt[9]=  10.1101;
double par0[10]={488.114,487.208,500.337,505.511,501.288,567.476,639.514,1465.69,6.05824e+06,500,};
double par1[10]={93.4051,82.4448,79.9232,89.7417,118.804,190.393,280.699,1513.15,2.95991e+07,30,};
double par2[10]={29.5189,29.8185,29.1204,28.6193,28.1688,22.8348,16.569,-67.1972,-1.37328e+06,30,};
double par3[10]={89.6642,84.4986,85.4173,90.8975,102.823,144.086,194.731,940.192,1.761e+07,50,};

*/

///////////////////////////////////////////////////////////////////////
  int event,realdata;
  int nVertices,nGoodVertices;
  float vtxZ[50],vtxZerr[50];
  float vtxY[50],vtxYerr[50];
  float vtxX[50],vtxXerr[50];
  int vtxisValid[50],vtxisFake[50];
  ///pfjets
 

  //calojets
  int caloNjets;
  float CALOjetMom[100];
  float CALOjetEta[100], CALOjetPhi[100],CALOjetPt[100],CALOjetPx[100],CALOjetPy[100],CALOjetPz[100],CALOjetE[100],CALOjetRawPt[100],CALOjetRawPx[100],CALOjetRawPy[100],CALOjetRawPz[100],CALOjetRawE[100],CALOHadEHF[100],CALOEmEHF[100];

  //calo7jets
 

  ///Gen particle
  int GenInx;
  int GenStatus[50],GenID[50];

 
  myTree.SetBranchAddress("event",  &event);
  myTree.SetBranchAddress("realdata",&realdata);
  myTree.SetBranchAddress("nVertices",&nVertices);
  myTree.SetBranchAddress("nGoodVertices",&nGoodVertices);
  myTree.SetBranchAddress("vtxX",vtxX);
  myTree.SetBranchAddress("vtxY",vtxY);
  myTree.SetBranchAddress("vtxZ",vtxZ);
  myTree.SetBranchAddress("vtxXerr",vtxXerr);
  myTree.SetBranchAddress("vtxYerr",vtxYerr);
  myTree.SetBranchAddress("vtxZerr",vtxZerr);
  myTree.SetBranchAddress("vtxisValid",vtxisValid);
  myTree.SetBranchAddress("vtxisFake",vtxisFake);


  myTree.SetBranchAddress("caloNjets",&caloNjets);
  myTree.SetBranchAddress("CALOjetEta",CALOjetEta);
  myTree.SetBranchAddress("CALOjetPhi",CALOjetPhi);
  myTree.SetBranchAddress("CALOjetPt",CALOjetPt);
  myTree.SetBranchAddress("CALOjetPx",CALOjetPx);
  myTree.SetBranchAddress("CALOjetPy",CALOjetPy);
  myTree.SetBranchAddress("CALOjetPz",CALOjetPz);
  myTree.SetBranchAddress("CALOjetE",CALOjetE);
  myTree.SetBranchAddress("CALOjetRawPt",CALOjetRawPt);
  myTree.SetBranchAddress("CALOjetRawPx",CALOjetRawPx);
  myTree.SetBranchAddress("CALOjetRawPy",CALOjetRawPy);
  myTree.SetBranchAddress("CALOjetRawPz",CALOjetRawPz);
  myTree.SetBranchAddress("CALOjetRawE",CALOjetRawE);
  myTree.SetBranchAddress("CALOjetMom",CALOjetMom);
  myTree.SetBranchAddress("CALOHadEHF",CALOHadEHF);
  myTree.SetBranchAddress("CALOEmEHF",CALOEmEHF);



  myTree.SetBranchAddress("GenInx", &GenInx);
  myTree.SetBranchAddress("GenStatus", GenStatus);
  myTree.SetBranchAddress("GenID", GenID);
  
  ////Histogram booking///////////






  theFile->cd();

TH1F* Wmass_UnCorr = new TH1F("Wmass_UnCorr","Wmass with Uncorrected jets",80,0,300);
Wmass_UnCorr->Sumw2(); 

 TH1F* Wmass_corr = new TH1F("Wmass_corr","Wmass corrected by quark driven constants",80,0,300);
  Wmass_corr->Sumw2();

TH1F* Wmass_UnCorr_hbhe = new TH1F("Wmass_UnCorr_hbhe","Wmass with Uncorrected jets in HBHE",80,0,300);
Wmass_UnCorr_hbhe->Sumw2(); 

 TH1F* Wmass_corr_hbhe = new TH1F("Wmass_corr_hbhe","Wmass corrected by quark driven constants in HBHE",80,0,300);
  Wmass_corr_hbhe->Sumw2();

TH1F* Wmass_UnCorr_hf = new TH1F("Wmass_UnCorr_hf","Wmass with Uncorrected jets in HF",80,0,300);
Wmass_UnCorr_hf->Sumw2(); 

 TH1F* Wmass_corr_hf = new TH1F("Wmass_corr_hf","Wmass corrected by quark driven constants in HF",80,0,300);
  Wmass_corr_hf->Sumw2();


  Int_t nevent = myTree.GetEntries();
  for (Int_t iev=0;iev<nevent;iev++) {
    if (iev%100000 == 0) cout<<iev<<"/"<<nevent<<endl;
    myTree.GetEntry(iev);      






  int common24 = 0;
   int index1[10],index2[10];
   for (int j = 0; j < caloNjets; ++j){  
     for (int jk = j; jk < caloNjets  ; ++jk){
	if (j == jk) continue;
	if(CALOjetMom[j] == 24 && CALOjetMom[jk] == 24 )
	  {
	    index1[common24] = j;
	    index2[common24] = jk;
	    common24++;
	}
     }
   }


   int ind1 = -99;
   int ind2 = -99;
    for (int gg = 0;gg<common24;++gg){
      ind1 = index1[gg];
      ind2 = index2[gg];
    }
  int common = 0;
      int index3[10],index4[10];
   for (int j = 0; j < caloNjets; ++j){  
     for (int jk = j; jk < caloNjets  ; ++jk){
	if (j == jk) continue;
	if(CALOjetMom[j] == -24 && CALOjetMom[jk] == -24 
	)
	  {
	    index3[common] = j;
	    index4[common] = jk;
	    common++;
	}
     }
   }
double corr_fact=-999.;
   int ind3 = -99;
   int ind4 = -99;
   int corrjetE1= -99;
   int corrjetPx1= -99; 
   int corrjetPy1= -99; 
   int corrjetPz1= -99; 

   int corrjetE2= -99;
   int corrjetPx2= -99; 
   int corrjetPy2= -99; 
   int corrjetPz2= -99; 

   int corrjetE3= -99;
   int corrjetPx3= -99; 
   int corrjetPy3= -99; 
   int corrjetPz3= -99; 

   int corrjetE4= -99;
   int corrjetPx4= -99; 
   int corrjetPy4= -99; 
   int corrjetPz4= -99; 


    for (int gg = 0;gg<common;++gg){
      ind3 = index3[gg];
      ind4 = index4[gg];
    }

if (CALOjetRawPt[ind1]<0 ||CALOjetRawPt[ind2]<0)continue;
//cout<<CALOjetRawPt[ind1]<<endl;


//double eta_ind[11]={0.,0.5,1.1,1.6,2.1,2.6,3.1,3.6,4.2,4.7,5.2};

   if ((ind1 != -99 && ind2 != -99)){
//cout<<"new"<<endl;
for(int i=0; i<10;i++){//eta index
//if(fabs(CALOjetEta[ind1])>2.8||fabs(CALOjetEta[ind2])>3.2/*||fabs(CALOjetRawPt[ind1])>60.||fabs(CALOjetRawPt[ind2])>60.*/) continue;
//if(c[i][j]==0)continue;
if(fabs(CALOjetEta[ind1])>0.519*i && fabs(CALOjetEta[ind1])<0.519*(i+1) ){
/*if (CALOjetRawPt[ind1]<MAX_CALOjetRawPt[i])*/ corr_fact= par0[i]/(par1[i] + par2[i]* log(par3[i]*CALOjetRawPt[ind1]));
//if (CALOjetRawPt[ind1]>MAX_CALOjetRawPt[i]) corr_fact= par0[i]/(par1[i] + par2[i]* log(par3[i]*MAX_CALOjetRawPt[i]));
corrjetE1=corr_fact*CALOjetRawE[ind1];
corrjetPx1=corr_fact*CALOjetRawPx[ind1];
corrjetPy1=corr_fact*CALOjetRawPy[ind1];
corrjetPz1=corr_fact*CALOjetRawPz[ind1];
}

if( fabs(CALOjetEta[ind1])>0.519*i && fabs(CALOjetEta[ind1])<0.519*(i+1)){
/*if (CALOjetRawPt[ind2]<MAX_CALOjetRawPt[i])*/ corr_fact= par0[i]/(par1[i] + par2[i]* log(par3[i]*CALOjetRawPt[ind2]));
//if (CALOjetRawPt[ind2]>MAX_CALOjetRawPt[i]) corr_fact= par0[i]/(par1[i] + par2[i]* log(par3[i]*MAX_CALOjetRawPt[i]));
corrjetE2=corr_fact*CALOjetRawE[ind2];
corrjetPx2=corr_fact*CALOjetRawPx[ind2];
corrjetPy2=corr_fact*CALOjetRawPy[ind2];
corrjetPz2=corr_fact*CALOjetRawPz[ind2];
}
}
}

double WrawMass1=-99;
double WrawMass2=-99;
double WcorrMass1=-99;
double WcorrMass2=-99;


if((ind1!=-99 && ind2!=-99&&corrjetE1!=-99&&corrjetE2!=-99&&corrjetPx1!=-99&&corrjetPx2!=-99)){

WrawMass1 = sqrt(fabs(pow(CALOjetRawE[ind1]+CALOjetRawE[ind2],2)-pow(CALOjetRawPx[ind1]+CALOjetRawPx[ind2],2)-pow(CALOjetRawPy[ind1]+CALOjetRawPy[ind2],2)-pow(CALOjetRawPz[ind1]+CALOjetRawPz[ind2],2)));

WcorrMass1 = sqrt(fabs(pow(corrjetE1+corrjetE2,2)-pow(corrjetPx1+corrjetPx2,2)-pow(corrjetPy1+corrjetPy2,2)-pow(corrjetPz1+corrjetPz2,2)));

Wmass_UnCorr->Fill(WrawMass1);

Wmass_corr->Fill(WcorrMass1);


if(CALOjetEta[ind1]<2.8 &&CALOjetEta[ind2]<2.8){
Wmass_UnCorr_hbhe->Fill(WrawMass1);
Wmass_corr_hbhe->Fill(WcorrMass1);
}

if(CALOjetEta[ind1]>2.8 &&CALOjetEta[ind2]>2.8){
Wmass_UnCorr_hf->Fill(WrawMass1);
Wmass_corr_hf->Fill(WcorrMass1);
}
}


if((ind3 != -99 && ind4 != -99)){

WrawMass2 = sqrt(fabs(pow(CALOjetRawE[ind3]+CALOjetRawE[ind4],2)-pow(CALOjetRawPx[ind3]+CALOjetRawPx[ind4],2)-pow(CALOjetRawPy[ind3]+CALOjetRawPy[ind4],2)-pow(CALOjetRawPz[ind3]+CALOjetRawPz[ind4],2)));

WcorrMass2 = sqrt(fabs(pow(corrjetE3+corrjetE4,2)-pow(corrjetPx3+corrjetPx4,2)-pow(corrjetPy3+corrjetPy4,2)-pow(corrjetPz3+corrjetPz4,2)));



}



}
  theFile->Write();
  theFile->Close();

}//end void  

int main(int argc,char **argv){
  tree1r();
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
