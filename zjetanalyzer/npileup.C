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
void npilup()
{


  gStyle->SetPalette(1);
 



  TChain myTree("demo/Tree");






myTree.Add("/tmp/bbilin/mc_30_03_2012.root");
TFile *theFile = new TFile("npileup.root","RECREATE");



TH1::AddDirectory(true);

  int event,run,lumi,bxnumber,realdata;
//  int hlt_trigger_fired;


 
  float TrueNumInteractions;
  //met 


  int nVertices, tr_1, tr_2, tr_3;
  
  myTree.SetBranchAddress("event",  &event);
  myTree.SetBranchAddress("run", &run);
  myTree.SetBranchAddress("lumi", &lumi);
  myTree.SetBranchAddress("TrueNumInteractions", &TrueNumInteractions);
  theFile->cd();



TH1D *npil = new TH1D ("npil","npil",50,0,50);

 Int_t nevent = myTree.GetEntries();

  for (Int_t iev=0;iev<nevent;iev++) {




 	  if (iev%100000 == 0) cout<<iev<<"/"<<nevent<<endl;
    	myTree.GetEntry(iev); 

npil->Fill(TrueNumInteractions);







}//end of event loop


  theFile->Write();
  theFile->Close();


}//end void  

int main(int argc,char **argv){
  npilup();
}

