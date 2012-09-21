{

 //  gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

TFile * file[16];

file[0] = TFile::Open("rootfiles_new/21_09_data.root");

file[1] = TFile::Open("rootfiles_new/21_09_mc.root");

file[2] = TFile::Open("rootfiles_new/21_09_bg_ztautau.root");

file[3] = TFile::Open("rootfiles_new/21_09_bg_tt.root");

file[4] = TFile::Open("rootfiles_new/21_09_bg_wjets.root");

file[5] = TFile::Open("rootfiles_new/21_09_ww.root");
file[6] = TFile::Open("rootfiles_new/21_09_wz.root");
file[7] = TFile::Open("rootfiles_new/21_09_zz.root");

file[8] = TFile::Open("rootfiles_new/21_09_mc_low.root");

TH1F *h_[200][200];
TH1F *h1_[200][200];
TH1F *h2_[200][200];
TH1F *h12_[200][200];
int numberoffiles = 9;
int number_of_histos = 68;

for(int ii=0;ii<numberoffiles;ii++){
h2_[ii][0] =(TH1F*) file[ii] ->Get("h_numberofPFjets_1_0");
h_[ii][0]= (TH1F*)h2_[ii][0]->Clone();
h_[ii][0]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity of all electrons with mva>0 with cuts");
h_[ii][0]->GetYaxis()->SetTitle("Entries");

h2_[ii][1] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_1_0");
h_[ii][1]= (TH1F*)h2_[ii][1]->Clone();
h_[ii][1]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE of all electrons with mva>0 with cuts");
h_[ii][1]->GetYaxis()->SetTitle("Entries");

h2_[ii][2] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_1_0");
h_[ii][2]= (TH1F*)h2_[ii][2]->Clone();
h_[ii][2]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF of all electrons with mva>0 with cuts");
h_[ii][2]->GetYaxis()->SetTitle("Entries");

h2_[ii][3] =(TH1F*) file[ii] ->Get("h_mz_1_0");
h_[ii][3]= (TH1F*)h2_[ii][3]->Clone();
h_[ii][3]->GetXaxis()->SetTitle("Z Mass of all electrons with mva>0 with cuts");
h_[ii][3]->GetYaxis()->SetTitle("Entries");
h_[ii][3]->GetXaxis()->SetRangeUser(0.,200.);

h2_[ii][4] =(TH1F*) file[ii] ->Get("h_mz_1j_all_1_0");
h_[ii][4]= (TH1F*)h2_[ii][4]->Clone();
h_[ii][4]->GetXaxis()->SetTitle("Z mass for excl. Z+1 jet events of all electrons with mva>0 with cuts");
h_[ii][4]->GetYaxis()->SetTitle("Entries");

h2_[ii][5] =(TH1F*) file[ii] ->Get("h_mz_1j_c_1_0");
h_[ii][5]= (TH1F*)h2_[ii][5]->Clone();
h_[ii][5]->GetXaxis()->SetTitle("Z mass for excl. Z+1 central jet events of all electrons with mva>0 with cuts");
h_[ii][5]->GetYaxis()->SetTitle("Entries");

h2_[ii][6] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_1_0");
h_[ii][6]= (TH1F*)h2_[ii][6]->Clone();
h_[ii][6]->GetXaxis()->SetTitle("Z mass for excl. Z+1 HF jet events of all electrons with mva>0 with cuts");
h_[ii][6]->GetYaxis()->SetTitle("Entries");

h2_[ii][7] =(TH1F*) file[ii] ->Get("h_dielec_PT_1_0");
h_[ii][7]= (TH1F*)h2_[ii][7]->Clone();
h_[ii][7]->GetXaxis()->SetTitle("Z P_{T}  of all electrons with mva>0 with cuts");
h_[ii][7]->GetYaxis()->SetTitle("Entries");

h2_[ii][8] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_1_0");
h_[ii][8]= (TH1F*)h2_[ii][8]->Clone();
h_[ii][8]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events of all electrons with mva>0 with cuts");
h_[ii][8]->GetYaxis()->SetTitle("Entries");

h2_[ii][9] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_1_0");
h_[ii][9]= (TH1F*)h2_[ii][9]->Clone();
h_[ii][9]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events of all electrons with mva>0 with cuts");
h_[ii][9]->GetYaxis()->SetTitle("Entries");

h2_[ii][10] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_1_0");
h_[ii][10]= (TH1F*)h2_[ii][10]->Clone();
h_[ii][10]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events of all electrons with mva>0 with cuts");
h_[ii][10]->GetYaxis()->SetTitle("Entries");

h2_[ii][11] =(TH1F*) file[ii] ->Get("h_dielecphi_1_0");
h_[ii][11]= (TH1F*)h2_[ii][11]->Clone();
h_[ii][11]->GetXaxis()->SetTitle("Z #phi of all electrons with mva>0 with cuts");
h_[ii][11]->GetYaxis()->SetTitle("Entries");

h2_[ii][12] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_1_0");
h_[ii][12]= (TH1F*)h2_[ii][12]->Clone();
h_[ii][12]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events of all electrons with mva>0 with cuts");
h_[ii][12]->GetYaxis()->SetTitle("Entries");

h2_[ii][13] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_1_0");
h_[ii][13]= (TH1F*)h2_[ii][13]->Clone();
h_[ii][13]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events of all electrons with mva>0 with cuts");
h_[ii][13]->GetYaxis()->SetTitle("Entries");

h2_[ii][14] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_1_0");
h_[ii][14]= (TH1F*)h2_[ii][14]->Clone();
h_[ii][14]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events of all electrons with mva>0 with cuts");
h_[ii][14]->GetYaxis()->SetTitle("Entries");

h2_[ii][15] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1_0");
h_[ii][15]= (TH1F*)h2_[ii][15]->Clone();
h_[ii][15]->GetXaxis()->SetTitle("Z Y of all electrons with mva>0 with cuts");
h_[ii][15]->GetYaxis()->SetTitle("Entries");

h2_[ii][16] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_1_0");
h_[ii][16]= (TH1F*)h2_[ii][16]->Clone();
h_[ii][16]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events of all electrons with mva>0 with cuts");
h_[ii][16]->GetYaxis()->SetTitle("Entries");

h2_[ii][17] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_1_0");
h_[ii][17]= (TH1F*)h2_[ii][17]->Clone();
h_[ii][17]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events of all electrons with mva>0 with cuts");
h_[ii][17]->GetYaxis()->SetTitle("Entries");

h2_[ii][18] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_1_0");
h_[ii][18]= (TH1F*)h2_[ii][18]->Clone();
h_[ii][18]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events of all electrons with mva>0 with cuts");
h_[ii][18]->GetYaxis()->SetTitle("Entries");

h2_[ii][19] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_0");
h_[ii][19]= (TH1F*)h2_[ii][19]->Clone();
h_[ii][19]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets of all electrons with mva>0 with cuts");
h_[ii][19]->GetYaxis()->SetTitle("Entries");

h2_[ii][20] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_0");
h_[ii][20]= (TH1F*)h2_[ii][20]->Clone();
h_[ii][20]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets of all electrons with mva>0 with cuts");
h_[ii][20]->GetYaxis()->SetTitle("Entries");

h2_[ii][21] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_1_0");
h_[ii][21]= (TH1F*)h2_[ii][21]->Clone();
h_[ii][21]->GetXaxis()->SetTitle("# of good vertices of Z events of all electrons with mva>0 with cuts");
h_[ii][21]->GetYaxis()->SetTitle("Entries");

h2_[ii][22] =(TH1F*) file[ii] ->Get("h_mz_same_sign_1_0");
h_[ii][22]= (TH1F*)h2_[ii][22]->Clone();
h_[ii][22]->GetXaxis()->SetTitle("same sign electron invariant mass of all electrons with mva>0 with cuts");
h_[ii][22]->GetYaxis()->SetTitle("Entries");

h2_[ii][23] =(TH1F*) file[ii] ->Get("h_jet_pt_all_1_0");
h_[ii][23]= (TH1F*)h2_[ii][23]->Clone();
h_[ii][23]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector of all electrons with mva>0 with cuts");
h_[ii][23]->GetYaxis()->SetTitle("Entries");

h2_[ii][24] =(TH1F*) file[ii] ->Get("h_jet_eta_all_1_0");
h_[ii][24]= (TH1F*)h2_[ii][24]->Clone();
h_[ii][24]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector of all electrons with mva>0 with cuts");
h_[ii][24]->GetYaxis()->SetTitle("Entries");

h2_[ii][25] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_1_0");
h_[ii][25]= (TH1F*)h2_[ii][25]->Clone();
h_[ii][25]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  of all electrons with mva>0 with cuts");
h_[ii][25]->GetYaxis()->SetTitle("Entries");

h2_[ii][26] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_1_0");
h_[ii][26]= (TH1F*)h2_[ii][26]->Clone();
h_[ii][26]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  of all electrons with mva>0 with cuts");
h_[ii][26]->GetYaxis()->SetTitle("Entries");

h2_[ii][27] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_1_0");
h_[ii][27]= (TH1F*)h2_[ii][27]->Clone();
h_[ii][27]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  of all electrons with mva>0 with cuts");
h_[ii][27]->GetYaxis()->SetTitle("Entries");

h2_[ii][28] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_1_0");
h_[ii][28]= (TH1F*)h2_[ii][28]->Clone();
h_[ii][28]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  of all electrons with mva>0 with cuts");
h_[ii][28]->GetYaxis()->SetTitle("Entries");

h2_[ii][29] =(TH1F*) file[ii] ->Get("h_eEta_mva0_do02_pfiso01_mhits0");
h_[ii][29]= (TH1F*)h2_[ii][29]->Clone();
h_[ii][29]->GetXaxis()->SetTitle("Electron Eta for electrons with mva>0 with cuts");
h_[ii][29]->GetYaxis()->SetTitle("Entries");

h2_[ii][30] =(TH1F*) file[ii] ->Get("h_ePhi_mva0_do02_pfiso01_mhits0");
h_[ii][30]= (TH1F*)h2_[ii][30]->Clone();
h_[ii][30]->GetXaxis()->SetTitle("Electron Phi for electrons with mva>0 with cuts");
h_[ii][30]->GetYaxis()->SetTitle("Entries");

h2_[ii][31] =(TH1F*) file[ii] ->Get("h_ePt_mva0_do02_pfiso01_mhits0");
h_[ii][31]= (TH1F*)h2_[ii][31]->Clone();
h_[ii][31]->GetXaxis()->SetTitle("Electron Pt for electrons with mva>0 with cuts");
h_[ii][31]->GetYaxis()->SetTitle("Entries");

h2_[ii][32] =(TH1F*) file[ii] ->Get("h_ePx_mva0_do02_pfiso01_mhits0");
h_[ii][32]= (TH1F*)h2_[ii][32]->Clone();
h_[ii][32]->GetXaxis()->SetTitle("Electron Px for electrons with mva>0 with cuts");
h_[ii][32]->GetYaxis()->SetTitle("Entries");

h2_[ii][33] =(TH1F*) file[ii] ->Get("h_ePy_mva0_do02_pfiso01_mhits0");
h_[ii][33]= (TH1F*)h2_[ii][33]->Clone();
h_[ii][33]->GetXaxis()->SetTitle("Electron Py for electrons with mva>0 with cuts");
h_[ii][33]->GetYaxis()->SetTitle("Entries");

h2_[ii][34] =(TH1F*) file[ii] ->Get("h_ePz_mva0_do02_pfiso01_mhits0");
h_[ii][34]= (TH1F*)h2_[ii][34]->Clone();
h_[ii][34]->GetXaxis()->SetTitle("Electron Pz for electrons with mva>0 with cuts");
h_[ii][34]->GetYaxis()->SetTitle("Entries");

h2_[ii][35] =(TH1F*) file[ii] ->Get("h_eM_mva0_do02_pfiso01_mhits0");
h_[ii][35]= (TH1F*)h2_[ii][35]->Clone();
h_[ii][35]->GetXaxis()->SetTitle("Electron M for electrons with mva>0 with cuts");
h_[ii][35]->GetYaxis()->SetTitle("Entries");

h2_[ii][36] =(TH1F*) file[ii] ->Get("h_eE_mva0_do02_pfiso01_mhits0");
h_[ii][36]= (TH1F*)h2_[ii][36]->Clone();
h_[ii][36]->GetXaxis()->SetTitle("Electron Energy for electrons with mva>0 with cuts");
h_[ii][36]->GetYaxis()->SetTitle("Entries");

h2_[ii][37] =(TH1F*) file[ii] ->Get("h_eMVATrigId_mva0_do02_pfiso01_mhits0");
h_[ii][37]= (TH1F*)h2_[ii][37]->Clone();
h_[ii][37]->GetXaxis()->SetTitle("Electron MVA Trig Id for electrons with mva>0 with cuts");
h_[ii][37]->GetYaxis()->SetTitle("Entries");

h2_[ii][38] =(TH1F*) file[ii] ->Get("h_escSigmaIEtaIEta_mva0_do02_pfiso01_mhits0");
h_[ii][38]= (TH1F*)h2_[ii][38]->Clone();
h_[ii][38]->GetXaxis()->SetTitle("Electron SuperCluster #sigma i #eta i #etafor electrons with mva>0 with cuts");
h_[ii][38]->GetYaxis()->SetTitle("Entries");

h2_[ii][39] =(TH1F*) file[ii] ->Get("h_edeltaPhiSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0");
h_[ii][39]= (TH1F*)h2_[ii][39]->Clone();
h_[ii][39]->GetXaxis()->SetTitle("Electron SuperCluster #Delta #phi for electrons with mva>0 with cuts");
h_[ii][39]->GetYaxis()->SetTitle("Entries");

h2_[ii][40] =(TH1F*) file[ii] ->Get("h_edeltaEtaSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0");
h_[ii][40]= (TH1F*)h2_[ii][40]->Clone();
h_[ii][40]->GetXaxis()->SetTitle("Electron SuperCluster #Delta #eta for electrons with mva>0 with cuts");
h_[ii][40]->GetYaxis()->SetTitle("Entries");

h2_[ii][41] =(TH1F*) file[ii] ->Get("h_ehadronicOverEm_mva0_do02_pfiso01_mhits0");
h_[ii][41]= (TH1F*)h2_[ii][41]->Clone();
h_[ii][41]->GetXaxis()->SetTitle("Electron E/M for electrons with mva>0 with cuts");
h_[ii][41]->GetYaxis()->SetTitle("Entries");

h2_[ii][42] =(TH1F*) file[ii] ->Get("h_egsfTrack_numberOfLostHits_mva0_do02_pfiso01_mhits0");
h_[ii][42]= (TH1F*)h2_[ii][42]->Clone();
h_[ii][42]->GetXaxis()->SetTitle("Electron number of lost hits for electrons with mva>0 with cuts");
h_[ii][42]->GetYaxis()->SetTitle("Entries");

h2_[ii][43] =(TH1F*) file[ii] ->Get("h_ed0vtx_mva0_do02_pfiso01_mhits0");
h_[ii][43]= (TH1F*)h2_[ii][43]->Clone();
h_[ii][43]->GetXaxis()->SetTitle("Electron d0 for electrons with mva>0 with cuts");
h_[ii][43]->GetYaxis()->SetTitle("Entries");

h2_[ii][44] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0");
h_[ii][44]= (TH1F*)h2_[ii][44]->Clone();
h_[ii][44]->GetXaxis()->SetTitle("Electron dz for electrons with mva>0 with cuts");
h_[ii][44]->GetYaxis()->SetTitle("Entries");

h2_[ii][45] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0");
h_[ii][45]= (TH1F*)h2_[ii][45]->Clone();
h_[ii][45]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0");
h_[ii][45]->GetYaxis()->SetTitle("Entries");

h2_[ii][46] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0");
h_[ii][46]= (TH1F*)h2_[ii][46]->Clone();
h_[ii][46]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0");
h_[ii][46]->GetYaxis()->SetTitle("Entries");

h2_[ii][47] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0");
h_[ii][47]= (TH1F*)h2_[ii][47]->Clone();
h_[ii][47]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0");
h_[ii][47]->GetYaxis()->SetTitle("Entries");

h2_[ii][48] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0");
h_[ii][48]= (TH1F*)h2_[ii][48]->Clone();
h_[ii][48]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0");
h_[ii][48]->GetYaxis()->SetTitle("Entries");

h2_[ii][49] =(TH1F*) file[ii] ->Get("h_erelIso_mva0_do02_pfiso01_mhits0");
h_[ii][49]= (TH1F*)h2_[ii][49]->Clone();
h_[ii][49]->GetXaxis()->SetTitle("Electron PF iso for electrons with mva>0 with cuts");
h_[ii][49]->GetYaxis()->SetTitle("Entries");

h2_[ii][50] =(TH1F*) file[ii] ->Get("h_erelIsodb_mva0_do02_pfiso01_mhits0");
h_[ii][50]= (TH1F*)h2_[ii][50]->Clone();
h_[ii][50]->GetXaxis()->SetTitle("db corrected Electron PF iso for electrons with mva>0 with cuts");
h_[ii][50]->GetYaxis()->SetTitle("Entries");

h2_[ii][51] =(TH1F*) file[ii] ->Get("h_erelIsorho_mva0_do02_pfiso01_mhits0");
h_[ii][51]= (TH1F*)h2_[ii][51]->Clone();
h_[ii][51]->GetXaxis()->SetTitle("EA corrected Electron PF iso for electrons with mva>0 with cuts");
h_[ii][51]->GetYaxis()->SetTitle("Entries");

h2_[ii][52] =(TH1F*) file[ii] ->Get("h_efMVAVar_fbrem_mva0_do02_pfiso01_mhits0");
h_[ii][52]= (TH1F*)h2_[ii][52]->Clone();
h_[ii][52]->GetXaxis()->SetTitle("Electron fbrem(mva variable) for electrons with mva>0 with cuts");
h_[ii][52]->GetYaxis()->SetTitle("Entries");

h2_[ii][53] =(TH1F*) file[ii] ->Get("h_efMVAVar_kfchi2_mva0_do02_pfiso01_mhits0");
h_[ii][53]= (TH1F*)h2_[ii][53]->Clone();
h_[ii][53]->GetXaxis()->SetTitle("Electron kfchi^2(mva variable) for electrons with mva>0 with cuts");
h_[ii][53]->GetYaxis()->SetTitle("Entries");

h2_[ii][54] =(TH1F*) file[ii] ->Get("h_efMVAVar_kfhits_mva0_do02_pfiso01_mhits0");
h_[ii][54]= (TH1F*)h2_[ii][54]->Clone();
h_[ii][54]->GetXaxis()->SetTitle("Electron kfhits(mva variable) for electrons with mva>0 with cuts");
h_[ii][54]->GetYaxis()->SetTitle("Entries");

h2_[ii][55] =(TH1F*) file[ii] ->Get("h_efMVAVar_gsfchi2_mva0_do02_pfiso01_mhits0");
h_[ii][55]= (TH1F*)h2_[ii][55]->Clone();
h_[ii][55]->GetXaxis()->SetTitle("Electron gsfchi^2(mva variable) for electrons with mva>0 with cuts");
h_[ii][55]->GetYaxis()->SetTitle("Entries");

h2_[ii][56] =(TH1F*) file[ii] ->Get("h_efMVAVar_detacalo_mva0_do02_pfiso01_mhits0");
h_[ii][56]= (TH1F*)h2_[ii][56]->Clone();
h_[ii][56]->GetXaxis()->SetTitle("Electron d #eta Calo(mva variable) for electrons with mva>0 with cuts");
h_[ii][56]->GetYaxis()->SetTitle("Entries");

h2_[ii][57] =(TH1F*) file[ii] ->Get("h_efMVAVar_see_mva0_do02_pfiso01_mhits0");
h_[ii][57]= (TH1F*)h2_[ii][57]->Clone();
h_[ii][57]->GetXaxis()->SetTitle("Electron  #sigma i #eta i #eta(mva variable) for electrons with mva>0 with cuts");
h_[ii][57]->GetYaxis()->SetTitle("Entries");

h2_[ii][58] =(TH1F*) file[ii] ->Get("h_efMVAVar_spp_mva0_do02_pfiso01_mhits0");
h_[ii][58]= (TH1F*)h2_[ii][58]->Clone();
h_[ii][58]->GetXaxis()->SetTitle("Electron  #sigma i #phi i #phi(mva variable) for electrons with mva>0 with cuts");
h_[ii][58]->GetYaxis()->SetTitle("Entries");

h2_[ii][59] =(TH1F*) file[ii] ->Get("h_efMVAVar_etawidth_mva0_do02_pfiso01_mhits0");
h_[ii][59]= (TH1F*)h2_[ii][59]->Clone();
h_[ii][59]->GetXaxis()->SetTitle("Electron eta width(mva variable) for electrons with mva>0 with cuts");
h_[ii][59]->GetYaxis()->SetTitle("Entries");

h2_[ii][60] =(TH1F*) file[ii] ->Get("h_efMVAVar_phiwidth_mva0_do02_pfiso01_mhits0");
h_[ii][60]= (TH1F*)h2_[ii][60]->Clone();
h_[ii][60]->GetXaxis()->SetTitle("Electron phi width(mva variable) for electrons with mva>0 with cuts");
h_[ii][60]->GetYaxis()->SetTitle("Entries");

h2_[ii][61] =(TH1F*) file[ii] ->Get("h_efMVAVar_e1x5e5x5_mva0_do02_pfiso01_mhits0");
h_[ii][61]= (TH1F*)h2_[ii][61]->Clone();
h_[ii][61]->GetXaxis()->SetTitle("Electron e1x5/e5x5(mva variable) for electrons with mva>0 with cuts");
h_[ii][61]->GetYaxis()->SetTitle("Entries");

h2_[ii][62] =(TH1F*) file[ii] ->Get("h_efMVAVar_R9_mva0_do02_pfiso01_mhits0");
h_[ii][62]= (TH1F*)h2_[ii][62]->Clone();
h_[ii][62]->GetXaxis()->SetTitle("Second Leading Electron MVA_mva0_do02_pfiso01_mhits0");
h_[ii][62]->GetYaxis()->SetTitle("Entries");

h2_[ii][63] =(TH1F*) file[ii] ->Get("h_efMVAVar_EoP_mva0_do02_pfiso01_mhits0");
h_[ii][63]= (TH1F*)h2_[ii][63]->Clone();
h_[ii][63]->GetXaxis()->SetTitle("Electron scE/p(mva variable) for electrons with mva>0 with cuts");
h_[ii][63]->GetYaxis()->SetTitle("Entries");

h2_[ii][64] =(TH1F*) file[ii] ->Get("h_efMVAVar_IoEmIoP_mva0_do02_pfiso01_mhits0");
h_[ii][64]= (TH1F*)h2_[ii][64]->Clone();
h_[ii][64]->GetXaxis()->SetTitle("Electron 1/E - 1/p(mva variable) for electrons with mva>0 with cuts");
h_[ii][64]->GetYaxis()->SetTitle("Entries");

h2_[ii][65] =(TH1F*) file[ii] ->Get("h_efMVAVar_eleEoPout_mva0_do02_pfiso01_mhits0");
h_[ii][65]= (TH1F*)h2_[ii][65]->Clone();
h_[ii][65]->GetXaxis()->SetTitle("Electron sc E/P out(mva variable) for electrons with mva>0 with cuts");
h_[ii][65]->GetYaxis()->SetTitle("Entries");

h2_[ii][66] =(TH1F*) file[ii] ->Get("h_efMVAVar_PreShowerOverRaw_mva0_do02_pfiso01_mhits0");
h_[ii][66]= (TH1F*)h2_[ii][66]->Clone();
h_[ii][66]->GetXaxis()->SetTitle("Electron PreShowerOverRaw(mva variable) for electrons with mva>0 with cuts");
h_[ii][66]->GetYaxis()->SetTitle("Entries");

h2_[ii][67] =(TH1F*) file[ii] ->Get("h_MZe_l_mva0_do02_pfiso01_mhits0");
h_[ii][67]= (TH1F*)h2_[ii][67]->Clone();
h_[ii][67]->GetXaxis()->SetTitle("Di-Electron invariant mass for electrons with mva>0 with cuts");
h_[ii][67]->GetYaxis()->SetTitle("Entries");
h_[ii][67]->GetXaxis()->SetRangeUser(50.,130.);
}

char name1[200];
char name2[200];
// TH1F *hsum[100];
TCanvas *c[200];
TH1F*qcd_bckg_0;
TH1F*qcd_bckg_1;
TH1F*qcd_bckg_4;
TH1F*qcd_bckg_5;
TH1F *Mont[200];
for(int ii = 0; ii<number_of_histos;ii++){

if(ii==20) continue;
//if(ii!=3 &&ii!=67)continue;
if(ii==58||ii==57||(ii<49 && ii>44)||ii==40||(ii<29&&ii>18 &&ii!=22)||ii==11||ii==7)continue;

bool logx=false;
bool logy=false;
bool integral_norm = false;
bool same_sign_0 = false;

if((ii>58&&ii<68&& ii!=62)||ii==56||ii==55||ii==54||ii==53||ii==51||ii==50||ii==49 ||(ii>30&&ii<45&& ii!=35)||ii==22||ii==29||(ii<19&&ii>=0&&(ii!=7||ii!=11||ii!=15))  )logy=true;

if(ii==66||ii==65||ii==63||ii==60||ii==59||ii==55||ii==51||ii==50||ii==49||ii==41)logx=true;

if(ii==21) same_sign_0=true;



//cout<<ii<<same_sign_0<<endl;
c[ii] = new TCanvas();
c[ii]->Divide(1,2,0.01,0);


h_[0][ii]->Sumw2();
if(!integral_norm){
h_[1][ii]->Sumw2();
h_[1][ii]->Scale((3503.75*5097.)/30461028.);
h_[2][ii]->Sumw2();
h_[2][ii]->Scale((3503.75*5097.)/30461028.);
h_[3][ii]->Sumw2();
h_[3][ii]->Scale((225.197*5097.)/6675889.);
h_[4][ii]->Sumw2();
h_[4][ii]->Scale((36257.2*5097.)/13562800.);
h_[5][ii]->Sumw2();
h_[5][ii]->Scale((54.838*5097.)/9931470.);
h_[6][ii]->Sumw2();
h_[6][ii]->Scale((32.31*5097.)/9988254.);
h_[7][ii]->Sumw2();
h_[7][ii]->Scale((17.654*5097.)/9652754.);

h_[8][ii]->Scale((1000.25*5097.)/6858991.);
}


if(integral_norm){
h_[1][ii]->Scale(h_[0][ii]->Integral()/h_[1][ii]->Integral());
}
if(!integral_norm){
TH1F *diboson = (TH1F*)h_[5][ii]->Clone();
		diboson ->Add(h_[6][ii]);
		diboson ->Add(h_[7][ii]);
TH1F* bg_sum = (TH1F*)h_[2][ii]->Clone();
bg_sum->Add(h_[3][ii]);
bg_sum->Add(h_[4][ii]);
bg_sum->Add(diboson);
TH1F*mc_bg = (TH1F*)h_[1][ii]->Clone();
mc_bg->Add(bg_sum);
//		mc_bg->Add(diboson);
//		mc_bg->Add(h_[2][ii]);
//		mc_bg->Add(h_[3][ii]);
//		mc_bg->Add(h_[4][ii]);
//		mc_bg->Add(h_[8][ii]);
}

if(integral_norm)TH1F*mc_bg = (TH1F*)h_[1][ii]->Clone();

Mont[ii]= (TH1F*)h_[0][ii]->Clone();

Mont[ii]->Divide(h_[0][ii],mc_bg,1.0,1.0);



sprintf(name1,"histo_z_plot_new/hist_%i.root",ii);
sprintf(name2,"histo_z_plot_new/hist_%i.pdf",ii);
          c[ii]->cd(1); 

mc_bg->Draw("hhist");
mc_bg->SetTitle("Z/#gamma* -> ee + BG");
mc_bg.SetLineColor(4);
mc_bg.SetLineWidth(2);

if(logx)gPad->SetLogx();
if(logy)gPad->SetLogy();
if(!integral_norm){

h_[1][ii]->Draw("hhistsames");
h_[1][ii]->SetTitle("Z/#gamma* -> ee");
h_[1][ii].SetLineStyle(3);
h_[1][ii].SetLineColor(12);
h_[1][ii].SetLineWidth(2);

bg_sum->Draw("hhistsames");
bg_sum->SetTitle("BG");
bg_sum.SetLineColor(8);
bg_sum.SetLineWidth(2);
/*
h_[2][ii]->Draw("hhistsames");
h_[2][ii]->SetTitle("Z/#gamma* -> #tau #tau");
h_[2][ii].SetLineColor(5);
h_[2][ii].SetLineWidth(2);
h_[3][ii]->Draw("hhistsames");
h_[3][ii]->SetTitle("t #bar{t}");
h_[3][ii].SetLineColor(6);
h_[3][ii].SetLineWidth(2);
h_[4][ii]->Draw("hhistsames");
h_[4][ii]->SetTitle("W + jets");
h_[4][ii].SetLineColor(7);
h_[4][ii].SetLineWidth(2);

diboson->Draw("hhistsames");
diboson->SetTitle("Diboson");
diboson.SetLineColor(8);
diboson.SetLineWidth(2);
*/
}
h_[0][ii]->Draw("e1sames");
h_[0][ii]->SetTitle("CMS Data 5.09 fb^{-1}");
h_[0][ii]->SetMarkerStyle(2);
//h_[0][ii]->SetLineWidth(3); 

if(same_sign_0){
qcd_bckg_0 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_0->Add(h_[1][ii],-1);
}

c[ii]->cd(2);
  gStyle->SetOptStat(kFALSE);
//	  Mont[ii]->Sumw2();
          Mont[ii]->SetMinimum(0.5);
	  Mont[ii]->SetMaximum(1.5);
          Mont[ii]->SetYTitle("Ratio");
	  Mont[ii]->Draw("e1");
if(logx)gPad->SetLogx();
//if(logy)gPad->SetLogy();
//TCanvas->Return(c);
c[ii]->Print(name1);
c[ii]->Print(name2);

}
/*
TCanvas *c_qcd_0 = new TCanvas();
c_qcd_0->cd();
qcd_bckg_0->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG) of MVA(e1)>0.0 MVA(e2)>0.0");
qcd_bckg_0->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_0->GetYaxis()->SetTitle("Entries");
qcd_bckg_0->Draw("E1");
c_qcd_0->Print("histo_z_plot/qcd_bckground_0.root");
c_qcd_0->Print("histo_z_plot/qcd_bckground_0.pdf");
*/


}
