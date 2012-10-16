{

 //  gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);
gROOT->SetBatch(kTRUE);

TFile * file[16];

file[0] = TFile::Open("rootfiles_new/16_10_52X_data.root");

file[1] = TFile::Open("rootfiles_new/16_10_52X_mc.root");

file[2] = TFile::Open("rootfiles_new/16_10_52X_bg_ztautau.root");

file[3] = TFile::Open("rootfiles_new/16_10_52X_bg_tt.root");

file[4] = TFile::Open("rootfiles_new/16_10_52X_bg_wjets.root");

file[5] = TFile::Open("rootfiles_new/16_10_52X_ww.root");
file[6] = TFile::Open("rootfiles_new/16_10_52X_wz.root");
file[7] = TFile::Open("rootfiles_new/16_10_52X_zz.root");

file[8] = TFile::Open("rootfiles_new/16_10_52X_mc_low.root");

TH1F *h_[500][500];
TH1F *h1_[500][500];
TH1F *h2_[500][500];
TH1F *h12_[500][500];
int numberoffiles = 9;


for(int ii=0;ii<numberoffiles;ii++){
h2_[ii][0] =(TH1F*) file[ii] ->Get("h_numberofPFjets_1_0");
h_[ii][0]= (TH1F*)h2_[ii][0]->Clone();
h_[ii][0]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity ");
h_[ii][0]->GetYaxis()->SetTitle("Entries");

h2_[ii][1] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_1_0");
h_[ii][1]= (TH1F*)h2_[ii][1]->Clone();
h_[ii][1]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE ");
h_[ii][1]->GetYaxis()->SetTitle("Entries");

h2_[ii][2] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_1_0");
h_[ii][2]= (TH1F*)h2_[ii][2]->Clone();
h_[ii][2]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF ");
h_[ii][2]->GetYaxis()->SetTitle("Entries");

h2_[ii][3] =(TH1F*) file[ii] ->Get("h_mz_1_0");
h_[ii][3]= (TH1F*)h2_[ii][3]->Clone();
h_[ii][3]->GetXaxis()->SetTitle("Z Mass ");
h_[ii][3]->GetYaxis()->SetTitle("Entries");
h_[ii][3]->GetXaxis()->SetRangeUser(0.,200.);

h2_[ii][4] =(TH1F*) file[ii] ->Get("h_mz_1j_all_1_0");
h_[ii][4]= (TH1F*)h2_[ii][4]->Clone();
h_[ii][4]->GetXaxis()->SetTitle("Z mass for excl. Z+1 jet events ");
h_[ii][4]->GetYaxis()->SetTitle("Entries");

h2_[ii][5] =(TH1F*) file[ii] ->Get("h_mz_1j_c_1_0");
h_[ii][5]= (TH1F*)h2_[ii][5]->Clone();
h_[ii][5]->GetXaxis()->SetTitle("Z mass for excl. Z+1 central jet events ");
h_[ii][5]->GetYaxis()->SetTitle("Entries");

h2_[ii][6] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_1_0");
h_[ii][6]= (TH1F*)h2_[ii][6]->Clone();
h_[ii][6]->GetXaxis()->SetTitle("Z mass for excl. Z+1 HF jet events ");
h_[ii][6]->GetYaxis()->SetTitle("Entries");

h2_[ii][7] =(TH1F*) file[ii] ->Get("h_dielec_PT_1_0");
h_[ii][7]= (TH1F*)h2_[ii][7]->Clone();
h_[ii][7]->GetXaxis()->SetTitle("Z P_{T}  ");
h_[ii][7]->GetYaxis()->SetTitle("Entries");

h2_[ii][8] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_1_0");
h_[ii][8]= (TH1F*)h2_[ii][8]->Clone();
h_[ii][8]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events ");
h_[ii][8]->GetYaxis()->SetTitle("Entries");

h2_[ii][9] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_1_0");
h_[ii][9]= (TH1F*)h2_[ii][9]->Clone();
h_[ii][9]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events ");
h_[ii][9]->GetYaxis()->SetTitle("Entries");

h2_[ii][10] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_1_0");
h_[ii][10]= (TH1F*)h2_[ii][10]->Clone();
h_[ii][10]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events ");
h_[ii][10]->GetYaxis()->SetTitle("Entries");

h2_[ii][11] =(TH1F*) file[ii] ->Get("h_dielecphi_1_0");
h_[ii][11]= (TH1F*)h2_[ii][11]->Clone();
h_[ii][11]->GetXaxis()->SetTitle("Z #phi ");
h_[ii][11]->GetYaxis()->SetTitle("Entries");
//h_[ii][11]->GetYaxis()->SetRangeUser(0.,50000.);

h2_[ii][12] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_1_0");
h_[ii][12]= (TH1F*)h2_[ii][12]->Clone();
h_[ii][12]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events ");
h_[ii][12]->GetYaxis()->SetTitle("Entries");

h2_[ii][13] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_1_0");
h_[ii][13]= (TH1F*)h2_[ii][13]->Clone();
h_[ii][13]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events ");
h_[ii][13]->GetYaxis()->SetTitle("Entries");

h2_[ii][14] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_1_0");
h_[ii][14]= (TH1F*)h2_[ii][14]->Clone();
h_[ii][14]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events ");
h_[ii][14]->GetYaxis()->SetTitle("Entries");

h2_[ii][15] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1_0");
h_[ii][15]= (TH1F*)h2_[ii][15]->Clone();
h_[ii][15]->GetXaxis()->SetTitle("Z Y ");
h_[ii][15]->GetYaxis()->SetTitle("Entries");

h2_[ii][16] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_1_0");
h_[ii][16]= (TH1F*)h2_[ii][16]->Clone();
h_[ii][16]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events ");
h_[ii][16]->GetYaxis()->SetTitle("Entries");

h2_[ii][17] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_1_0");
h_[ii][17]= (TH1F*)h2_[ii][17]->Clone();
h_[ii][17]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events ");
h_[ii][17]->GetYaxis()->SetTitle("Entries");

h2_[ii][18] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_1_0");
h_[ii][18]= (TH1F*)h2_[ii][18]->Clone();
h_[ii][18]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events ");
h_[ii][18]->GetYaxis()->SetTitle("Entries");

h2_[ii][19] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_0");
h_[ii][19]= (TH1F*)h2_[ii][19]->Clone();
h_[ii][19]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets ");
h_[ii][19]->GetYaxis()->SetTitle("Entries");

h2_[ii][20] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_0");
h_[ii][20]= (TH1F*)h2_[ii][20]->Clone();
h_[ii][20]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets ");
h_[ii][20]->GetYaxis()->SetTitle("Entries");

h2_[ii][21] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_1_0");
h_[ii][21]= (TH1F*)h2_[ii][21]->Clone();
h_[ii][21]->GetXaxis()->SetTitle("# of good vertices of Z events ");
h_[ii][21]->GetYaxis()->SetTitle("Entries");

h2_[ii][22] =(TH1F*) file[ii] ->Get("h_mz_same_sign_1_0");
h_[ii][22]= (TH1F*)h2_[ii][22]->Clone();
h_[ii][22]->GetXaxis()->SetTitle("same sign electron invariant mass with PFIso<0.1");
h_[ii][22]->GetYaxis()->SetTitle("Entries");

h2_[ii][23] =(TH1F*) file[ii] ->Get("h_jet_pt_all_1_0");
h_[ii][23]= (TH1F*)h2_[ii][23]->Clone();
h_[ii][23]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector ");
h_[ii][23]->GetYaxis()->SetTitle("Entries");

h2_[ii][24] =(TH1F*) file[ii] ->Get("h_jet_eta_all_1_0");
h_[ii][24]= (TH1F*)h2_[ii][24]->Clone();
h_[ii][24]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector ");
h_[ii][24]->GetYaxis()->SetTitle("Entries");

h2_[ii][25] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_1_0");
h_[ii][25]= (TH1F*)h2_[ii][25]->Clone();
h_[ii][25]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  ");
h_[ii][25]->GetYaxis()->SetTitle("Entries");

h2_[ii][26] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_1_0");
h_[ii][26]= (TH1F*)h2_[ii][26]->Clone();
h_[ii][26]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  ");
h_[ii][26]->GetYaxis()->SetTitle("Entries");

h2_[ii][27] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_1_0");
h_[ii][27]= (TH1F*)h2_[ii][27]->Clone();
h_[ii][27]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  ");
h_[ii][27]->GetYaxis()->SetTitle("Entries");

h2_[ii][28] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_1_0");
h_[ii][28]= (TH1F*)h2_[ii][28]->Clone();
h_[ii][28]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  ");
h_[ii][28]->GetYaxis()->SetTitle("Entries");

h2_[ii][29] =(TH1F*) file[ii] ->Get("h_eEta_mva0_do02_pfiso01_mhits0");
h_[ii][29]= (TH1F*)h2_[ii][29]->Clone();
h_[ii][29]->GetXaxis()->SetTitle("Electron Eta  ");
h_[ii][29]->GetYaxis()->SetTitle("Entries");

h2_[ii][30] =(TH1F*) file[ii] ->Get("h_ePhi_mva0_do02_pfiso01_mhits0");
h_[ii][30]= (TH1F*)h2_[ii][30]->Clone();
h_[ii][30]->GetXaxis()->SetTitle("Electron Phi  ");
h_[ii][30]->GetYaxis()->SetTitle("Entries");

h2_[ii][31] =(TH1F*) file[ii] ->Get("h_ePt_mva0_do02_pfiso01_mhits0");
h_[ii][31]= (TH1F*)h2_[ii][31]->Clone();
h_[ii][31]->GetXaxis()->SetTitle("Electron Pt  ");
h_[ii][31]->GetYaxis()->SetTitle("Entries");

h2_[ii][32] =(TH1F*) file[ii] ->Get("h_ePx_mva0_do02_pfiso01_mhits0");
h_[ii][32]= (TH1F*)h2_[ii][32]->Clone();
h_[ii][32]->GetXaxis()->SetTitle("Electron Px  ");
h_[ii][32]->GetYaxis()->SetTitle("Entries");

h2_[ii][33] =(TH1F*) file[ii] ->Get("h_ePy_mva0_do02_pfiso01_mhits0");
h_[ii][33]= (TH1F*)h2_[ii][33]->Clone();
h_[ii][33]->GetXaxis()->SetTitle("Electron Py  ");
h_[ii][33]->GetYaxis()->SetTitle("Entries");

h2_[ii][34] =(TH1F*) file[ii] ->Get("h_ePz_mva0_do02_pfiso01_mhits0");
h_[ii][34]= (TH1F*)h2_[ii][34]->Clone();
h_[ii][34]->GetXaxis()->SetTitle("Electron Pz  ");
h_[ii][34]->GetYaxis()->SetTitle("Entries");

h2_[ii][35] =(TH1F*) file[ii] ->Get("h_eM_mva0_do02_pfiso01_mhits0");
h_[ii][35]= (TH1F*)h2_[ii][35]->Clone();
h_[ii][35]->GetXaxis()->SetTitle("Electron M  ");
h_[ii][35]->GetYaxis()->SetTitle("Entries");

h2_[ii][36] =(TH1F*) file[ii] ->Get("h_eE_mva0_do02_pfiso01_mhits0");
h_[ii][36]= (TH1F*)h2_[ii][36]->Clone();
h_[ii][36]->GetXaxis()->SetTitle("Electron Energy  ");
h_[ii][36]->GetYaxis()->SetTitle("Entries");

h2_[ii][37] =(TH1F*) file[ii] ->Get("h_eMVATrigId_mva0_do02_pfiso01_mhits0");
h_[ii][37]= (TH1F*)h2_[ii][37]->Clone();
h_[ii][37]->GetXaxis()->SetTitle("Electron MVA Trig Id  ");
h_[ii][37]->GetYaxis()->SetTitle("Entries");

h2_[ii][38] =(TH1F*) file[ii] ->Get("h_escSigmaIEtaIEta_mva0_do02_pfiso01_mhits0");
h_[ii][38]= (TH1F*)h2_[ii][38]->Clone();
h_[ii][38]->GetXaxis()->SetTitle("Electron SuperCluster #sigma i #eta i #eta ");
h_[ii][38]->GetYaxis()->SetTitle("Entries");

h2_[ii][39] =(TH1F*) file[ii] ->Get("h_edeltaPhiSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0");
h_[ii][39]= (TH1F*)h2_[ii][39]->Clone();
h_[ii][39]->GetXaxis()->SetTitle("Electron SuperCluster #Delta #phi  ");
h_[ii][39]->GetYaxis()->SetTitle("Entries");

h2_[ii][40] =(TH1F*) file[ii] ->Get("h_edeltaEtaSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0");
h_[ii][40]= (TH1F*)h2_[ii][40]->Clone();
h_[ii][40]->GetXaxis()->SetTitle("Electron SuperCluster #Delta #eta  ");
h_[ii][40]->GetYaxis()->SetTitle("Entries");

h2_[ii][41] =(TH1F*) file[ii] ->Get("h_ehadronicOverEm_mva0_do02_pfiso01_mhits0");
h_[ii][41]= (TH1F*)h2_[ii][41]->Clone();
h_[ii][41]->GetXaxis()->SetTitle("Electron E/M  ");
h_[ii][41]->GetYaxis()->SetTitle("Entries");

h2_[ii][42] =(TH1F*) file[ii] ->Get("h_egsfTrack_numberOfLostHits_mva0_do02_pfiso01_mhits0");
h_[ii][42]= (TH1F*)h2_[ii][42]->Clone();
h_[ii][42]->GetXaxis()->SetTitle("Electron number of lost hits  ");
h_[ii][42]->GetYaxis()->SetTitle("Entries");

h2_[ii][43] =(TH1F*) file[ii] ->Get("h_ed0vtx_mva0_do02_pfiso01_mhits0");
h_[ii][43]= (TH1F*)h2_[ii][43]->Clone();
h_[ii][43]->GetXaxis()->SetTitle("Electron d0  ");
h_[ii][43]->GetYaxis()->SetTitle("Entries");

h2_[ii][44] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0");
h_[ii][44]= (TH1F*)h2_[ii][44]->Clone();
h_[ii][44]->GetXaxis()->SetTitle("Electron dz  ");
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
h_[ii][49]->GetXaxis()->SetTitle("Electron PF iso  ");
h_[ii][49]->GetYaxis()->SetTitle("Entries");

h2_[ii][50] =(TH1F*) file[ii] ->Get("h_erelIsodb_mva0_do02_pfiso01_mhits0");
h_[ii][50]= (TH1F*)h2_[ii][50]->Clone();
h_[ii][50]->GetXaxis()->SetTitle("db corrected Electron PF iso  ");
h_[ii][50]->GetYaxis()->SetTitle("Entries");

h2_[ii][51] =(TH1F*) file[ii] ->Get("h_erelIsorho_mva0_do02_pfiso01_mhits0");
h_[ii][51]= (TH1F*)h2_[ii][51]->Clone();
h_[ii][51]->GetXaxis()->SetTitle("EA corrected Electron PF iso  ");
h_[ii][51]->GetYaxis()->SetTitle("Entries");

h2_[ii][52] =(TH1F*) file[ii] ->Get("h_efMVAVar_fbrem_mva0_do02_pfiso01_mhits0");
h_[ii][52]= (TH1F*)h2_[ii][52]->Clone();
h_[ii][52]->GetXaxis()->SetTitle("Electron fbrem(mva variable)  ");
h_[ii][52]->GetYaxis()->SetTitle("Entries");

h2_[ii][53] =(TH1F*) file[ii] ->Get("h_efMVAVar_kfchi2_mva0_do02_pfiso01_mhits0");
h_[ii][53]= (TH1F*)h2_[ii][53]->Clone();
h_[ii][53]->GetXaxis()->SetTitle("Electron kfchi^2(mva variable)  ");
h_[ii][53]->GetYaxis()->SetTitle("Entries");

h2_[ii][54] =(TH1F*) file[ii] ->Get("h_efMVAVar_kfhits_mva0_do02_pfiso01_mhits0");
h_[ii][54]= (TH1F*)h2_[ii][54]->Clone();
h_[ii][54]->GetXaxis()->SetTitle("Electron kfhits(mva variable)  ");
h_[ii][54]->GetYaxis()->SetTitle("Entries");

h2_[ii][55] =(TH1F*) file[ii] ->Get("h_efMVAVar_gsfchi2_mva0_do02_pfiso01_mhits0");
h_[ii][55]= (TH1F*)h2_[ii][55]->Clone();
h_[ii][55]->GetXaxis()->SetTitle("Electron gsfchi^2(mva variable)  ");
h_[ii][55]->GetYaxis()->SetTitle("Entries");

h2_[ii][56] =(TH1F*) file[ii] ->Get("h_efMVAVar_detacalo_mva0_do02_pfiso01_mhits0");
h_[ii][56]= (TH1F*)h2_[ii][56]->Clone();
h_[ii][56]->GetXaxis()->SetTitle("Electron d #eta Calo(mva variable)  ");
h_[ii][56]->GetYaxis()->SetTitle("Entries");

h2_[ii][57] =(TH1F*) file[ii] ->Get("h_efMVAVar_see_mva0_do02_pfiso01_mhits0");
h_[ii][57]= (TH1F*)h2_[ii][57]->Clone();
h_[ii][57]->GetXaxis()->SetTitle("Electron  #sigma i #eta i #eta(mva variable)  ");
h_[ii][57]->GetYaxis()->SetTitle("Entries");

h2_[ii][58] =(TH1F*) file[ii] ->Get("h_efMVAVar_spp_mva0_do02_pfiso01_mhits0");
h_[ii][58]= (TH1F*)h2_[ii][58]->Clone();
h_[ii][58]->GetXaxis()->SetTitle("Electron  #sigma i #phi i #phi(mva variable)  ");
h_[ii][58]->GetYaxis()->SetTitle("Entries");

h2_[ii][59] =(TH1F*) file[ii] ->Get("h_efMVAVar_etawidth_mva0_do02_pfiso01_mhits0");
h_[ii][59]= (TH1F*)h2_[ii][59]->Clone();
h_[ii][59]->GetXaxis()->SetTitle("Electron eta width(mva variable)  ");
h_[ii][59]->GetYaxis()->SetTitle("Entries");

h2_[ii][60] =(TH1F*) file[ii] ->Get("h_efMVAVar_phiwidth_mva0_do02_pfiso01_mhits0");
h_[ii][60]= (TH1F*)h2_[ii][60]->Clone();
h_[ii][60]->GetXaxis()->SetTitle("Electron phi width(mva variable)  ");
h_[ii][60]->GetYaxis()->SetTitle("Entries");

h2_[ii][61] =(TH1F*) file[ii] ->Get("h_efMVAVar_e1x5e5x5_mva0_do02_pfiso01_mhits0");
h_[ii][61]= (TH1F*)h2_[ii][61]->Clone();
h_[ii][61]->GetXaxis()->SetTitle("Electron e1x5/e5x5(mva variable)  ");
h_[ii][61]->GetYaxis()->SetTitle("Entries");

h2_[ii][62] =(TH1F*) file[ii] ->Get("h_efMVAVar_R9_mva0_do02_pfiso01_mhits0");
h_[ii][62]= (TH1F*)h2_[ii][62]->Clone();
h_[ii][62]->GetXaxis()->SetTitle("Second Leading Electron MVA_mva0_do02_pfiso01_mhits0");
h_[ii][62]->GetYaxis()->SetTitle("Entries");

h2_[ii][63] =(TH1F*) file[ii] ->Get("h_efMVAVar_EoP_mva0_do02_pfiso01_mhits0");
h_[ii][63]= (TH1F*)h2_[ii][63]->Clone();
h_[ii][63]->GetXaxis()->SetTitle("Electron scE/p(mva variable)  ");
h_[ii][63]->GetYaxis()->SetTitle("Entries");

h2_[ii][64] =(TH1F*) file[ii] ->Get("h_efMVAVar_IoEmIoP_mva0_do02_pfiso01_mhits0");
h_[ii][64]= (TH1F*)h2_[ii][64]->Clone();
h_[ii][64]->GetXaxis()->SetTitle("Electron 1/E - 1/p(mva variable)  ");
h_[ii][64]->GetYaxis()->SetTitle("Entries");

h2_[ii][65] =(TH1F*) file[ii] ->Get("h_efMVAVar_eleEoPout_mva0_do02_pfiso01_mhits0");
h_[ii][65]= (TH1F*)h2_[ii][65]->Clone();
h_[ii][65]->GetXaxis()->SetTitle("Electron sc E/P out(mva variable)  ");
h_[ii][65]->GetYaxis()->SetTitle("Entries");

h2_[ii][66] =(TH1F*) file[ii] ->Get("h_efMVAVar_PreShowerOverRaw_mva0_do02_pfiso01_mhits0");
h_[ii][66]= (TH1F*)h2_[ii][66]->Clone();
h_[ii][66]->GetXaxis()->SetTitle("Electron PreShowerOverRaw(mva variable)  ");
h_[ii][66]->GetYaxis()->SetTitle("Entries");

h2_[ii][67] =(TH1F*) file[ii] ->Get("h_MZe_l_mva0_do02_pfiso01_mhits0");
h_[ii][67]= (TH1F*)h2_[ii][67]->Clone();
h_[ii][67]->GetXaxis()->SetTitle("Di-Electron invariant mass  ");
h_[ii][67]->GetYaxis()->SetTitle("Entries");
h_[ii][67]->GetXaxis()->SetRangeUser(50.,130.);

h2_[ii][68] =(TH1F*) file[ii] ->Get("h_mz_same_sign_anti_iso_1_0");
h_[ii][68]= (TH1F*)h2_[ii][68]->Clone();
h_[ii][68]->GetXaxis()->SetTitle("same sign electron invariant mass with PFIso>0.1");
h_[ii][68]->GetYaxis()->SetTitle("Entries");

h2_[ii][69] =(TH1F*) file[ii] ->Get("h_mz_same_sign_anti_iso_oppsign_1_0");
h_[ii][69]= (TH1F*)h2_[ii][69]->Clone();
h_[ii][69]->GetXaxis()->SetTitle("opposite sign electron invariant mass with PFIso>0.1 ");
h_[ii][69]->GetYaxis()->SetTitle("Entries");

h2_[ii][70] =(TH1F*) file[ii] ->Get("h_mz_2j_all_1_0");
h_[ii][70]= (TH1F*)h2_[ii][70]->Clone();
h_[ii][70]->GetXaxis()->SetTitle("Z mass for excl. Z+2 jet events ");
h_[ii][70]->GetYaxis()->SetTitle("Entries");

h2_[ii][71] =(TH1F*) file[ii] ->Get("h_mz_2j_c_1_0");
h_[ii][71]= (TH1F*)h2_[ii][71]->Clone();
h_[ii][71]->GetXaxis()->SetTitle("Z mass for excl. Z+2 central jet events ");
h_[ii][71]->GetYaxis()->SetTitle("Entries");

h2_[ii][72] =(TH1F*) file[ii] ->Get("h_mz_2j_hf_1_0");
h_[ii][72]= (TH1F*)h2_[ii][72]->Clone();
h_[ii][72]->GetXaxis()->SetTitle("Z mass for excl. Z+2 forward jet events ");
h_[ii][72]->GetYaxis()->SetTitle("Entries");

h2_[ii][73] =(TH1F*) file[ii] ->Get("h_dielec_PT_2j_all_1_0");
h_[ii][73]= (TH1F*)h2_[ii][73]->Clone();
h_[ii][73]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+2 jet events ");
h_[ii][73]->GetYaxis()->SetTitle("Entries");

h2_[ii][74] =(TH1F*) file[ii] ->Get("h_dielec_PT_2j_c_1_0");
h_[ii][74]= (TH1F*)h2_[ii][74]->Clone();
h_[ii][74]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+2 central jet events ");
h_[ii][74]->GetYaxis()->SetTitle("Entries");

h2_[ii][75] =(TH1F*) file[ii] ->Get("h_dielec_PT_2j_hf_1_0");
h_[ii][75]= (TH1F*)h2_[ii][75]->Clone();
h_[ii][75]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+2 forward jet events ");
h_[ii][75]->GetYaxis()->SetTitle("Entries");

h2_[ii][76] =(TH1F*) file[ii] ->Get("h_dielecphi_2j_all_1_0");
h_[ii][76]= (TH1F*)h2_[ii][76]->Clone();
h_[ii][76]->GetXaxis()->SetTitle("Z #phi for excl. Z+2 jet events ");
h_[ii][76]->GetYaxis()->SetTitle("Entries");

h2_[ii][77] =(TH1F*) file[ii] ->Get("h_dielecphi_2j_c_1_0");
h_[ii][77]= (TH1F*)h2_[ii][77]->Clone();
h_[ii][77]->GetXaxis()->SetTitle("Z #phi for excl. Z+2 central jet events ");
h_[ii][77]->GetYaxis()->SetTitle("Entries");

h2_[ii][78] =(TH1F*) file[ii] ->Get("h_dielecphi_2j_hf_1_0");
h_[ii][78]= (TH1F*)h2_[ii][78]->Clone();
h_[ii][78]->GetXaxis()->SetTitle("Z #phi for excl. Z+2 forward jet events ");
h_[ii][78]->GetYaxis()->SetTitle("Entries");

h2_[ii][79] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_2j_all_1_0");
h_[ii][79]= (TH1F*)h2_[ii][79]->Clone();
h_[ii][79]->GetXaxis()->SetTitle("Z Y for excl. Z+2 jet events ");
h_[ii][79]->GetYaxis()->SetTitle("Entries");

h2_[ii][80] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_2j_c_1_0");
h_[ii][80]= (TH1F*)h2_[ii][80]->Clone();
h_[ii][80]->GetXaxis()->SetTitle("Z Y for excl. Z+2 central jet events ");
h_[ii][80]->GetYaxis()->SetTitle("Entries");

h2_[ii][81] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_2j_hf_1_0");
h_[ii][81]= (TH1F*)h2_[ii][81]->Clone();
h_[ii][81]->GetXaxis()->SetTitle("Z Y for excl. Z+2 forward jet events ");
h_[ii][81]->GetYaxis()->SetTitle("Entries");

h2_[ii][82] =(TH1F*) file[ii] ->Get("h_jet_pt_all_2j_1_0");
h_[ii][82]= (TH1F*)h2_[ii][82]->Clone();
h_[ii][82]->GetXaxis()->SetTitle("Jet pt of leading jet in Z +2jet events");
h_[ii][82]->GetYaxis()->SetTitle("Entries");

h2_[ii][83] =(TH1F*) file[ii] ->Get("h_jet_eta_all_2j_1_0");
h_[ii][83]= (TH1F*)h2_[ii][83]->Clone();
h_[ii][83]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z +2jet events");
h_[ii][83]->GetYaxis()->SetTitle("Entries");

h2_[ii][84] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_2j_1_0");
h_[ii][84]= (TH1F*)h2_[ii][84]->Clone();
h_[ii][84]->GetXaxis()->SetTitle("Jet pt of leading jet in Z +2 central jet events");
h_[ii][84]->GetYaxis()->SetTitle("Entries");

h2_[ii][85] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_2j_1_0");
h_[ii][85]= (TH1F*)h2_[ii][85]->Clone();
h_[ii][85]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z +2 central jet events");
h_[ii][85]->GetYaxis()->SetTitle("Entries");

h2_[ii][86] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_2j_1_0");
h_[ii][86]= (TH1F*)h2_[ii][86]->Clone();
h_[ii][86]->GetXaxis()->SetTitle("Jet pt of leading jet in Z +2 forward jet events");
h_[ii][86]->GetYaxis()->SetTitle("Entries");

h2_[ii][87] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_2j_1_0");
h_[ii][87]= (TH1F*)h2_[ii][87]->Clone();
h_[ii][87]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z +2 forward jet events");
h_[ii][87]->GetYaxis()->SetTitle("Entries");

h2_[ii][88] =(TH1F*) file[ii] ->Get("h_mz_3j_all_1_0");
h_[ii][88]= (TH1F*)h2_[ii][88]->Clone();
h_[ii][88]->GetXaxis()->SetTitle("Z mass for Z + 3 jet events");
h_[ii][88]->GetYaxis()->SetTitle("Entries");

h2_[ii][89] =(TH1F*) file[ii] ->Get("h_mz_3j_c_1_0");
h_[ii][89]= (TH1F*)h2_[ii][89]->Clone();
h_[ii][89]->GetXaxis()->SetTitle("Z mass for Z + 3 central jet events");
h_[ii][89]->GetYaxis()->SetTitle("Entries");

h2_[ii][90] =(TH1F*) file[ii] ->Get("h_mz_3j_hf_1_0");
h_[ii][90]= (TH1F*)h2_[ii][90]->Clone();
h_[ii][90]->GetXaxis()->SetTitle("Z mass for Z + 3 forward jet events");
h_[ii][90]->GetYaxis()->SetTitle("Entries");

h2_[ii][91] =(TH1F*) file[ii] ->Get("h_dielec_PT_3j_all_1_0");
h_[ii][91]= (TH1F*)h2_[ii][91]->Clone();
h_[ii][91]->GetXaxis()->SetTitle("Z P_{T} for Z + 3 jet events");
h_[ii][91]->GetYaxis()->SetTitle("Entries");

h2_[ii][92] =(TH1F*) file[ii] ->Get("h_dielec_PT_3j_c_1_0");
h_[ii][92]= (TH1F*)h2_[ii][92]->Clone();
h_[ii][92]->GetXaxis()->SetTitle("Z P_{T} for Z + 3 central jet events");
h_[ii][92]->GetYaxis()->SetTitle("Entries");

h2_[ii][93] =(TH1F*) file[ii] ->Get("h_dielec_PT_3j_hf_1_0");
h_[ii][93]= (TH1F*)h2_[ii][93]->Clone();
h_[ii][93]->GetXaxis()->SetTitle("Z P_{T} for Z + 3 forward jet events");
h_[ii][93]->GetYaxis()->SetTitle("Entries");

h2_[ii][94] =(TH1F*) file[ii] ->Get("h_dielecphi_3j_all_1_0");
h_[ii][94]= (TH1F*)h2_[ii][94]->Clone();
h_[ii][94]->GetXaxis()->SetTitle("Z #phi for Z + 3 jet events");
h_[ii][94]->GetYaxis()->SetTitle("Entries");

h2_[ii][95] =(TH1F*) file[ii] ->Get("h_dielecphi_3j_c_1_0");
h_[ii][95]= (TH1F*)h2_[ii][95]->Clone();
h_[ii][95]->GetXaxis()->SetTitle("Z #phi for Z + 3 central jet events");
h_[ii][95]->GetYaxis()->SetTitle("Entries");

h2_[ii][96] =(TH1F*) file[ii] ->Get("h_dielecphi_3j_hf_1_0");
h_[ii][96]= (TH1F*)h2_[ii][96]->Clone();
h_[ii][96]->GetXaxis()->SetTitle("Z #phi for Z + 3 forward jet events");
h_[ii][96]->GetYaxis()->SetTitle("Entries");

h2_[ii][97] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_3j_all_1_0");
h_[ii][97]= (TH1F*)h2_[ii][97]->Clone();
h_[ii][97]->GetXaxis()->SetTitle("Z Y for Z + 3 jet events");
h_[ii][97]->GetYaxis()->SetTitle("Entries");

h2_[ii][98] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_3j_c_1_0");
h_[ii][98]= (TH1F*)h2_[ii][98]->Clone();
h_[ii][98]->GetXaxis()->SetTitle("Z Y for Z + 3 central jet events");
h_[ii][98]->GetYaxis()->SetTitle("Entries");

h2_[ii][99] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_3j_hf_1_0");
h_[ii][99]= (TH1F*)h2_[ii][99]->Clone();
h_[ii][99]->GetXaxis()->SetTitle("Z Y for Z + 3 forward jet events");
h_[ii][99]->GetYaxis()->SetTitle("Entries");

h2_[ii][100] =(TH1F*) file[ii] ->Get("h_jet_pt_all_3j_1_0");
h_[ii][100]= (TH1F*)h2_[ii][100]->Clone();
h_[ii][100]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 3 jet events");
h_[ii][100]->GetYaxis()->SetTitle("Entries");

h2_[ii][101] =(TH1F*) file[ii] ->Get("h_jet_eta_all_3j_1_0");
h_[ii][101]= (TH1F*)h2_[ii][101]->Clone();
h_[ii][101]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 3 jet events");
h_[ii][101]->GetYaxis()->SetTitle("Entries");

h2_[ii][102] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_3j_1_0");
h_[ii][102]= (TH1F*)h2_[ii][102]->Clone();
h_[ii][102]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 3 central jet events");
h_[ii][102]->GetYaxis()->SetTitle("Entries");

h2_[ii][103] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_3j_1_0");
h_[ii][103]= (TH1F*)h2_[ii][103]->Clone();
h_[ii][103]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 3 central jet events");
h_[ii][103]->GetYaxis()->SetTitle("Entries");

h2_[ii][104] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_3j_1_0");
h_[ii][104]= (TH1F*)h2_[ii][104]->Clone();
h_[ii][104]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 3 forward jet events");
h_[ii][104]->GetYaxis()->SetTitle("Entries");

h2_[ii][105] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_3j_1_0");
h_[ii][105]= (TH1F*)h2_[ii][105]->Clone();
h_[ii][105]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 3 forward jet events");
h_[ii][105]->GetYaxis()->SetTitle("Entries");

h2_[ii][106] =(TH1F*) file[ii] ->Get("h_mz_4j_all_1_0");
h_[ii][106]= (TH1F*)h2_[ii][106]->Clone();
h_[ii][106]->GetXaxis()->SetTitle("Z mass for Z + 4 jet events");
h_[ii][106]->GetYaxis()->SetTitle("Entries");

h2_[ii][107] =(TH1F*) file[ii] ->Get("h_mz_4j_c_1_0");
h_[ii][107]= (TH1F*)h2_[ii][107]->Clone();
h_[ii][107]->GetXaxis()->SetTitle("Z mass for Z + 4 central jet events");
h_[ii][107]->GetYaxis()->SetTitle("Entries");

h2_[ii][108] =(TH1F*) file[ii] ->Get("h_mz_4j_hf_1_0");
h_[ii][108]= (TH1F*)h2_[ii][108]->Clone();
h_[ii][108]->GetXaxis()->SetTitle("Z mass for Z + 4 forward jet events");
h_[ii][108]->GetYaxis()->SetTitle("Entries");

h2_[ii][109] =(TH1F*) file[ii] ->Get("h_dielec_PT_4j_all_1_0");
h_[ii][109]= (TH1F*)h2_[ii][109]->Clone();
h_[ii][109]->GetXaxis()->SetTitle("Z P_{T} for Z + 4 jet events");
h_[ii][109]->GetYaxis()->SetTitle("Entries");

h2_[ii][110] =(TH1F*) file[ii] ->Get("h_dielec_PT_4j_c_1_0");
h_[ii][110]= (TH1F*)h2_[ii][110]->Clone();
h_[ii][110]->GetXaxis()->SetTitle("Z P_{T} for Z + 4 central jet events");
h_[ii][110]->GetYaxis()->SetTitle("Entries");

h2_[ii][111] =(TH1F*) file[ii] ->Get("h_dielec_PT_4j_hf_1_0");
h_[ii][111]= (TH1F*)h2_[ii][111]->Clone();
h_[ii][111]->GetXaxis()->SetTitle("Z P_{T} for Z + 4 forward jet events");
h_[ii][111]->GetYaxis()->SetTitle("Entries");

h2_[ii][112] =(TH1F*) file[ii] ->Get("h_dielecphi_4j_all_1_0");
h_[ii][112]= (TH1F*)h2_[ii][112]->Clone();
h_[ii][112]->GetXaxis()->SetTitle("Z #phi for Z + 4 jet events");
h_[ii][112]->GetYaxis()->SetTitle("Entries");

h2_[ii][113] =(TH1F*) file[ii] ->Get("h_dielecphi_4j_c_1_0");
h_[ii][113]= (TH1F*)h2_[ii][113]->Clone();
h_[ii][113]->GetXaxis()->SetTitle("Z #phi for Z + 4 central jet events");
h_[ii][113]->GetYaxis()->SetTitle("Entries");

h2_[ii][114] =(TH1F*) file[ii] ->Get("h_dielecphi_4j_hf_1_0");
h_[ii][114]= (TH1F*)h2_[ii][114]->Clone();
h_[ii][114]->GetXaxis()->SetTitle("Z #phi for Z + 4 forward jet events");
h_[ii][114]->GetYaxis()->SetTitle("Entries");

h2_[ii][115] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_4j_all_1_0");
h_[ii][115]= (TH1F*)h2_[ii][115]->Clone();
h_[ii][115]->GetXaxis()->SetTitle("Z Y for Z + 4 jet events");
h_[ii][115]->GetYaxis()->SetTitle("Entries");

h2_[ii][116] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_4j_c_1_0");
h_[ii][116]= (TH1F*)h2_[ii][116]->Clone();
h_[ii][116]->GetXaxis()->SetTitle("Z Y for Z + 4 central jet events");
h_[ii][116]->GetYaxis()->SetTitle("Entries");

h2_[ii][117] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_4j_hf_1_0");
h_[ii][117]= (TH1F*)h2_[ii][117]->Clone();
h_[ii][117]->GetXaxis()->SetTitle("Z Y for Z + 4 forward jet events");
h_[ii][117]->GetYaxis()->SetTitle("Entries");

h2_[ii][118] =(TH1F*) file[ii] ->Get("h_jet_pt_all_4j_1_0");
h_[ii][118]= (TH1F*)h2_[ii][118]->Clone();
h_[ii][118]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 4 jet events");
h_[ii][118]->GetYaxis()->SetTitle("Entries");

h2_[ii][119] =(TH1F*) file[ii] ->Get("h_jet_eta_all_4j_1_0");
h_[ii][119]= (TH1F*)h2_[ii][119]->Clone();
h_[ii][119]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 4 jet events");
h_[ii][119]->GetYaxis()->SetTitle("Entries");

h2_[ii][120] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_4j_1_0");
h_[ii][120]= (TH1F*)h2_[ii][120]->Clone();
h_[ii][120]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 4 central jet events");
h_[ii][120]->GetYaxis()->SetTitle("Entries");

h2_[ii][121] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_4j_1_0");
h_[ii][121]= (TH1F*)h2_[ii][121]->Clone();
h_[ii][121]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 4 central jet events");
h_[ii][121]->GetYaxis()->SetTitle("Entries");

h2_[ii][122] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_4j_1_0");
h_[ii][122]= (TH1F*)h2_[ii][122]->Clone();
h_[ii][122]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 4 forward jet events");
h_[ii][122]->GetYaxis()->SetTitle("Entries");

h2_[ii][123] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_4j_1_0");
h_[ii][123]= (TH1F*)h2_[ii][123]->Clone();
h_[ii][123]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 4 forward jet events");
h_[ii][123]->GetYaxis()->SetTitle("Entries");

h2_[ii][124] =(TH1F*) file[ii] ->Get("h_mz_5j_all_1_0");
h_[ii][124]= (TH1F*)h2_[ii][124]->Clone();
h_[ii][124]->GetXaxis()->SetTitle("Z mass for Z + 5 jet events");
h_[ii][124]->GetYaxis()->SetTitle("Entries");

h2_[ii][125] =(TH1F*) file[ii] ->Get("h_mz_5j_c_1_0");
h_[ii][125]= (TH1F*)h2_[ii][125]->Clone();
h_[ii][125]->GetXaxis()->SetTitle("Z mass for Z + 5 central jet events");
h_[ii][125]->GetYaxis()->SetTitle("Entries");

h2_[ii][126] =(TH1F*) file[ii] ->Get("h_mz_5j_hf_1_0");
h_[ii][126]= (TH1F*)h2_[ii][126]->Clone();
h_[ii][126]->GetXaxis()->SetTitle("Z mass for Z + 5 forward jet events");
h_[ii][126]->GetYaxis()->SetTitle("Entries");

h2_[ii][127] =(TH1F*) file[ii] ->Get("h_dielec_PT_5j_all_1_0");
h_[ii][127]= (TH1F*)h2_[ii][127]->Clone();
h_[ii][127]->GetXaxis()->SetTitle("Z P_{T} for Z + 5 jet events");
h_[ii][127]->GetYaxis()->SetTitle("Entries");

h2_[ii][128] =(TH1F*) file[ii] ->Get("h_dielec_PT_5j_c_1_0");
h_[ii][128]= (TH1F*)h2_[ii][128]->Clone();
h_[ii][128]->GetXaxis()->SetTitle("Z P_{T} for Z + 5 central jet events");
h_[ii][128]->GetYaxis()->SetTitle("Entries");

h2_[ii][129] =(TH1F*) file[ii] ->Get("h_dielec_PT_5j_hf_1_0");
h_[ii][129]= (TH1F*)h2_[ii][129]->Clone();
h_[ii][129]->GetXaxis()->SetTitle("Z P_{T} for Z + 5 forward jet events");
h_[ii][129]->GetYaxis()->SetTitle("Entries");

h2_[ii][130] =(TH1F*) file[ii] ->Get("h_dielecphi_5j_all_1_0");
h_[ii][130]= (TH1F*)h2_[ii][130]->Clone();
h_[ii][130]->GetXaxis()->SetTitle("Z #phi for Z + 5 jet events");
h_[ii][130]->GetYaxis()->SetTitle("Entries");

h2_[ii][131] =(TH1F*) file[ii] ->Get("h_dielecphi_5j_c_1_0");
h_[ii][131]= (TH1F*)h2_[ii][131]->Clone();
h_[ii][131]->GetXaxis()->SetTitle("Z #phi for Z + 5 central jet events");
h_[ii][131]->GetYaxis()->SetTitle("Entries");

h2_[ii][132] =(TH1F*) file[ii] ->Get("h_dielecphi_5j_hf_1_0");
h_[ii][132]= (TH1F*)h2_[ii][132]->Clone();
h_[ii][132]->GetXaxis()->SetTitle("Z #phi for Z + 5 forward jet events");
h_[ii][132]->GetYaxis()->SetTitle("Entries");

h2_[ii][133] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_5j_all_1_0");
h_[ii][133]= (TH1F*)h2_[ii][133]->Clone();
h_[ii][133]->GetXaxis()->SetTitle("Z Y for Z + 5 jet events");
h_[ii][133]->GetYaxis()->SetTitle("Entries");

h2_[ii][134] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_5j_c_1_0");
h_[ii][134]= (TH1F*)h2_[ii][134]->Clone();
h_[ii][134]->GetXaxis()->SetTitle("Z Y for Z + 5 central jet events");
h_[ii][134]->GetYaxis()->SetTitle("Entries");

h2_[ii][135] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_5j_hf_1_0");
h_[ii][135]= (TH1F*)h2_[ii][135]->Clone();
h_[ii][135]->GetXaxis()->SetTitle("Z Y for Z + 5 forward jet events");
h_[ii][135]->GetYaxis()->SetTitle("Entries");

h2_[ii][136] =(TH1F*) file[ii] ->Get("h_jet_pt_all_5j_1_0");
h_[ii][136]= (TH1F*)h2_[ii][136]->Clone();
h_[ii][136]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 5 jet events");
h_[ii][136]->GetYaxis()->SetTitle("Entries");

h2_[ii][137] =(TH1F*) file[ii] ->Get("h_jet_eta_all_5j_1_0");
h_[ii][137]= (TH1F*)h2_[ii][137]->Clone();
h_[ii][137]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 5 jet events");
h_[ii][137]->GetYaxis()->SetTitle("Entries");

h2_[ii][138] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_5j_1_0");
h_[ii][138]= (TH1F*)h2_[ii][138]->Clone();
h_[ii][138]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 5 central jet events");
h_[ii][138]->GetYaxis()->SetTitle("Entries");


h2_[ii][139] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_5j_1_0");
h_[ii][139]= (TH1F*)h2_[ii][139]->Clone();
h_[ii][139]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 5 central jet events");
h_[ii][139]->GetYaxis()->SetTitle("Entries");

h2_[ii][140] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_5j_1_0");
h_[ii][140]= (TH1F*)h2_[ii][140]->Clone();
h_[ii][140]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 5 forward jet events");
h_[ii][140]->GetYaxis()->SetTitle("Entries");

h2_[ii][141] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_5j_1_0");
h_[ii][141]= (TH1F*)h2_[ii][141]->Clone();
h_[ii][141]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 5 forward jet events");
h_[ii][141]->GetYaxis()->SetTitle("Entries");

h2_[ii][142] =(TH1F*) file[ii] ->Get("h_dielec_PT_same_sign_1_0");
h_[ii][142]= (TH1F*)h2_[ii][142]->Clone();
h_[ii][142]->GetXaxis()->SetTitle("Z P_{T}");
h_[ii][142]->GetYaxis()->SetTitle("Entries");

h2_[ii][143] =(TH1F*) file[ii] ->Get("h_numberofPFjets_1_1");
h_[ii][143]= (TH1F*)h2_[ii][143]->Clone();
h_[ii][143]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity");
h_[ii][143]->GetYaxis()->SetTitle("Entries");

h2_[ii][144] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_1_1");
h_[ii][144]= (TH1F*)h2_[ii][144]->Clone();
h_[ii][144]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE ");
h_[ii][144]->GetYaxis()->SetTitle("Entries");

h2_[ii][145] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_1_1");
h_[ii][145]= (TH1F*)h2_[ii][145]->Clone();
h_[ii][145]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF ");
h_[ii][145]->GetYaxis()->SetTitle("Entries");

h2_[ii][146] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_1_1");
h_[ii][146]= (TH1F*)h2_[ii][146]->Clone();
h_[ii][146]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events ");
h_[ii][146]->GetYaxis()->SetTitle("Entries");

h2_[ii][147] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_1_1");
h_[ii][147]= (TH1F*)h2_[ii][147]->Clone();
h_[ii][147]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events ");
h_[ii][147]->GetYaxis()->SetTitle("Entries");

h2_[ii][148] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_1_1");
h_[ii][148]= (TH1F*)h2_[ii][148]->Clone();
h_[ii][148]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events ");
h_[ii][148]->GetYaxis()->SetTitle("Entries");

h2_[ii][149] =(TH1F*) file[ii] ->Get("h_dielecphi_1_1");
h_[ii][149]= (TH1F*)h2_[ii][149]->Clone();
h_[ii][149]->GetXaxis()->SetTitle("Z #phi ");
h_[ii][149]->GetYaxis()->SetTitle("Entries");
//h_[ii][149]->GetYaxis()->SetRangeUser(0.,50000.);

h2_[ii][150] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_1_1");
h_[ii][150]= (TH1F*)h2_[ii][150]->Clone();
h_[ii][150]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events ");
h_[ii][150]->GetYaxis()->SetTitle("Entries");

h2_[ii][151] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_1_1");
h_[ii][151]= (TH1F*)h2_[ii][151]->Clone();
h_[ii][151]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events ");
h_[ii][151]->GetYaxis()->SetTitle("Entries");

h2_[ii][152] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_1_1");
h_[ii][152]= (TH1F*)h2_[ii][152]->Clone();
h_[ii][152]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events ");
h_[ii][152]->GetYaxis()->SetTitle("Entries");

h2_[ii][153] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1_1");
h_[ii][153]= (TH1F*)h2_[ii][153]->Clone();
h_[ii][153]->GetXaxis()->SetTitle("Z Y ");
h_[ii][153]->GetYaxis()->SetTitle("Entries");

h2_[ii][154] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_1_1");
h_[ii][154]= (TH1F*)h2_[ii][154]->Clone();
h_[ii][154]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events ");
h_[ii][154]->GetYaxis()->SetTitle("Entries");

h2_[ii][155] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_1_1");
h_[ii][155]= (TH1F*)h2_[ii][155]->Clone();
h_[ii][155]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events ");
h_[ii][155]->GetYaxis()->SetTitle("Entries");

h2_[ii][156] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_1_1");
h_[ii][156]= (TH1F*)h2_[ii][156]->Clone();
h_[ii][156]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events ");
h_[ii][156]->GetYaxis()->SetTitle("Entries");

h2_[ii][157] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_1");
h_[ii][157]= (TH1F*)h2_[ii][157]->Clone();
h_[ii][157]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets ");
h_[ii][157]->GetYaxis()->SetTitle("Entries");

h2_[ii][158] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_1");
h_[ii][158]= (TH1F*)h2_[ii][158]->Clone();
h_[ii][158]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets ");
h_[ii][158]->GetYaxis()->SetTitle("Entries");

h2_[ii][159] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_1_1");
h_[ii][159]= (TH1F*)h2_[ii][159]->Clone();
h_[ii][159]->GetXaxis()->SetTitle("# of good vertices of Z events ");
h_[ii][159]->GetYaxis()->SetTitle("Entries");

h2_[ii][160] =(TH1F*) file[ii] ->Get("h_mz_same_sign_1_1");
h_[ii][160]= (TH1F*)h2_[ii][160]->Clone();
h_[ii][160]->GetXaxis()->SetTitle("same sign electron invariant mass with PFIso<0.1");
h_[ii][160]->GetYaxis()->SetTitle("Entries");

h2_[ii][161] =(TH1F*) file[ii] ->Get("h_jet_pt_all_1_1");
h_[ii][161]= (TH1F*)h2_[ii][161]->Clone();
h_[ii][161]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector ");
h_[ii][161]->GetYaxis()->SetTitle("Entries");

h2_[ii][162] =(TH1F*) file[ii] ->Get("h_jet_eta_all_1_1");
h_[ii][162]= (TH1F*)h2_[ii][162]->Clone();
h_[ii][162]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector ");
h_[ii][162]->GetYaxis()->SetTitle("Entries");

h2_[ii][163] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_1_1");
h_[ii][163]= (TH1F*)h2_[ii][163]->Clone();
h_[ii][163]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  ");
h_[ii][163]->GetYaxis()->SetTitle("Entries");

h2_[ii][164] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_1_1");
h_[ii][164]= (TH1F*)h2_[ii][164]->Clone();
h_[ii][164]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  ");
h_[ii][164]->GetYaxis()->SetTitle("Entries");

h2_[ii][165] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_1_1");
h_[ii][165]= (TH1F*)h2_[ii][165]->Clone();
h_[ii][165]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  ");
h_[ii][165]->GetYaxis()->SetTitle("Entries");

h2_[ii][166] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_1_1");
h_[ii][166]= (TH1F*)h2_[ii][166]->Clone();
h_[ii][166]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  ");
h_[ii][166]->GetYaxis()->SetTitle("Entries");

h2_[ii][167] =(TH1F*) file[ii] ->Get("h_mz_2j_all_1_1");
h_[ii][167]= (TH1F*)h2_[ii][167]->Clone();
h_[ii][167]->GetXaxis()->SetTitle("Z mass for excl. Z+2 jet events ");
h_[ii][167]->GetYaxis()->SetTitle("Entries");

h2_[ii][168] =(TH1F*) file[ii] ->Get("h_mz_2j_c_1_1");
h_[ii][168]= (TH1F*)h2_[ii][168]->Clone();
h_[ii][168]->GetXaxis()->SetTitle("Z mass for excl. Z+2 central jet events ");
h_[ii][168]->GetYaxis()->SetTitle("Entries");

h2_[ii][169] =(TH1F*) file[ii] ->Get("h_mz_2j_hf_1_1");
h_[ii][169]= (TH1F*)h2_[ii][169]->Clone();
h_[ii][169]->GetXaxis()->SetTitle("Z mass for excl. Z+2 forward jet events ");
h_[ii][169]->GetYaxis()->SetTitle("Entries");

h2_[ii][170] =(TH1F*) file[ii] ->Get("h_dielec_PT_2j_all_1_1");
h_[ii][170]= (TH1F*)h2_[ii][170]->Clone();
h_[ii][170]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+2 jet events ");
h_[ii][170]->GetYaxis()->SetTitle("Entries");

h2_[ii][171] =(TH1F*) file[ii] ->Get("h_dielec_PT_2j_c_1_1");
h_[ii][171]= (TH1F*)h2_[ii][171]->Clone();
h_[ii][171]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+2 central jet events ");
h_[ii][171]->GetYaxis()->SetTitle("Entries");

h2_[ii][172] =(TH1F*) file[ii] ->Get("h_dielec_PT_2j_hf_1_1");
h_[ii][172]= (TH1F*)h2_[ii][172]->Clone();
h_[ii][172]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+2 forward jet events ");
h_[ii][172]->GetYaxis()->SetTitle("Entries");

h2_[ii][173] =(TH1F*) file[ii] ->Get("h_dielecphi_2j_all_1_1");
h_[ii][173]= (TH1F*)h2_[ii][173]->Clone();
h_[ii][173]->GetXaxis()->SetTitle("Z #phi for excl. Z+2 jet events ");
h_[ii][173]->GetYaxis()->SetTitle("Entries");

h2_[ii][174] =(TH1F*) file[ii] ->Get("h_dielecphi_2j_c_1_1");
h_[ii][174]= (TH1F*)h2_[ii][174]->Clone();
h_[ii][174]->GetXaxis()->SetTitle("Z #phi for excl. Z+2 central jet events ");
h_[ii][174]->GetYaxis()->SetTitle("Entries");

h2_[ii][175] =(TH1F*) file[ii] ->Get("h_dielecphi_2j_hf_1_1");
h_[ii][175]= (TH1F*)h2_[ii][175]->Clone();
h_[ii][175]->GetXaxis()->SetTitle("Z #phi for excl. Z+2 forward jet events ");
h_[ii][175]->GetYaxis()->SetTitle("Entries");

h2_[ii][176] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_2j_all_1_1");
h_[ii][176]= (TH1F*)h2_[ii][176]->Clone();
h_[ii][176]->GetXaxis()->SetTitle("Z Y for excl. Z+2 jet events ");
h_[ii][176]->GetYaxis()->SetTitle("Entries");

h2_[ii][177] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_2j_c_1_1");
h_[ii][177]= (TH1F*)h2_[ii][177]->Clone();
h_[ii][177]->GetXaxis()->SetTitle("Z Y for excl. Z+2 central jet events ");
h_[ii][177]->GetYaxis()->SetTitle("Entries");

h2_[ii][178] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_2j_hf_1_1");
h_[ii][178]= (TH1F*)h2_[ii][178]->Clone();
h_[ii][178]->GetXaxis()->SetTitle("Z Y for excl. Z+2 forward jet events ");
h_[ii][178]->GetYaxis()->SetTitle("Entries");

h2_[ii][179] =(TH1F*) file[ii] ->Get("h_jet_pt_all_2j_1_1");
h_[ii][179]= (TH1F*)h2_[ii][179]->Clone();
h_[ii][179]->GetXaxis()->SetTitle("Jet pt of leading jet in Z +2jet events");
h_[ii][179]->GetYaxis()->SetTitle("Entries");

h2_[ii][180] =(TH1F*) file[ii] ->Get("h_jet_eta_all_2j_1_1");
h_[ii][180]= (TH1F*)h2_[ii][180]->Clone();
h_[ii][180]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z +2jet events");
h_[ii][180]->GetYaxis()->SetTitle("Entries");

h2_[ii][181] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_2j_1_1");
h_[ii][181]= (TH1F*)h2_[ii][181]->Clone();
h_[ii][181]->GetXaxis()->SetTitle("Jet pt of leading jet in Z +2 central jet events");
h_[ii][181]->GetYaxis()->SetTitle("Entries");

h2_[ii][182] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_2j_1_1");
h_[ii][182]= (TH1F*)h2_[ii][182]->Clone();
h_[ii][182]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z +2 central jet events");
h_[ii][182]->GetYaxis()->SetTitle("Entries");

h2_[ii][183] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_2j_1_1");
h_[ii][183]= (TH1F*)h2_[ii][183]->Clone();
h_[ii][183]->GetXaxis()->SetTitle("Jet pt of leading jet in Z +2 forward jet events");
h_[ii][183]->GetYaxis()->SetTitle("Entries");

h2_[ii][184] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_2j_1_1");
h_[ii][184]= (TH1F*)h2_[ii][184]->Clone();
h_[ii][184]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z +2 forward jet events");
h_[ii][184]->GetYaxis()->SetTitle("Entries");

h2_[ii][185] =(TH1F*) file[ii] ->Get("h_mz_3j_all_1_1");
h_[ii][185]= (TH1F*)h2_[ii][185]->Clone();
h_[ii][185]->GetXaxis()->SetTitle("Z mass for Z + 3 jet events");
h_[ii][185]->GetYaxis()->SetTitle("Entries");

h2_[ii][186] =(TH1F*) file[ii] ->Get("h_mz_3j_c_1_1");
h_[ii][186]= (TH1F*)h2_[ii][186]->Clone();
h_[ii][186]->GetXaxis()->SetTitle("Z mass for Z + 3 central jet events");
h_[ii][186]->GetYaxis()->SetTitle("Entries");

h2_[ii][187] =(TH1F*) file[ii] ->Get("h_mz_3j_hf_1_1");
h_[ii][187]= (TH1F*)h2_[ii][187]->Clone();
h_[ii][187]->GetXaxis()->SetTitle("Z mass for Z + 3 forward jet events");
h_[ii][187]->GetYaxis()->SetTitle("Entries");

h2_[ii][188] =(TH1F*) file[ii] ->Get("h_dielec_PT_3j_all_1_1");
h_[ii][188]= (TH1F*)h2_[ii][188]->Clone();
h_[ii][188]->GetXaxis()->SetTitle("Z P_{T} for Z + 3 jet events");
h_[ii][188]->GetYaxis()->SetTitle("Entries");

h2_[ii][189] =(TH1F*) file[ii] ->Get("h_dielec_PT_3j_c_1_1");
h_[ii][189]= (TH1F*)h2_[ii][189]->Clone();
h_[ii][189]->GetXaxis()->SetTitle("Z P_{T} for Z + 3 central jet events");
h_[ii][189]->GetYaxis()->SetTitle("Entries");

h2_[ii][190] =(TH1F*) file[ii] ->Get("h_dielec_PT_3j_hf_1_1");
h_[ii][190]= (TH1F*)h2_[ii][190]->Clone();
h_[ii][190]->GetXaxis()->SetTitle("Z P_{T} for Z + 3 forward jet events");
h_[ii][190]->GetYaxis()->SetTitle("Entries");

h2_[ii][191] =(TH1F*) file[ii] ->Get("h_dielecphi_3j_all_1_1");
h_[ii][191]= (TH1F*)h2_[ii][191]->Clone();
h_[ii][191]->GetXaxis()->SetTitle("Z #phi for Z + 3 jet events");
h_[ii][191]->GetYaxis()->SetTitle("Entries");

h2_[ii][192] =(TH1F*) file[ii] ->Get("h_dielecphi_3j_c_1_1");
h_[ii][192]= (TH1F*)h2_[ii][192]->Clone();
h_[ii][192]->GetXaxis()->SetTitle("Z #phi for Z + 3 central jet events");
h_[ii][192]->GetYaxis()->SetTitle("Entries");

h2_[ii][193] =(TH1F*) file[ii] ->Get("h_dielecphi_3j_hf_1_1");
h_[ii][193]= (TH1F*)h2_[ii][193]->Clone();
h_[ii][193]->GetXaxis()->SetTitle("Z #phi for Z + 3 forward jet events");
h_[ii][193]->GetYaxis()->SetTitle("Entries");

h2_[ii][194] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_3j_all_1_1");
h_[ii][194]= (TH1F*)h2_[ii][194]->Clone();
h_[ii][194]->GetXaxis()->SetTitle("Z Y for Z + 3 jet events");
h_[ii][194]->GetYaxis()->SetTitle("Entries");

h2_[ii][195] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_3j_c_1_1");
h_[ii][195]= (TH1F*)h2_[ii][195]->Clone();
h_[ii][195]->GetXaxis()->SetTitle("Z Y for Z + 3 central jet events");
h_[ii][195]->GetYaxis()->SetTitle("Entries");

h2_[ii][196] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_3j_hf_1_1");
h_[ii][196]= (TH1F*)h2_[ii][196]->Clone();
h_[ii][196]->GetXaxis()->SetTitle("Z Y for Z + 3 forward jet events");
h_[ii][196]->GetYaxis()->SetTitle("Entries");

h2_[ii][197] =(TH1F*) file[ii] ->Get("h_jet_pt_all_3j_1_1");
h_[ii][197]= (TH1F*)h2_[ii][197]->Clone();
h_[ii][197]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 3 jet events");
h_[ii][197]->GetYaxis()->SetTitle("Entries");

h2_[ii][198] =(TH1F*) file[ii] ->Get("h_jet_eta_all_3j_1_1");
h_[ii][198]= (TH1F*)h2_[ii][198]->Clone();
h_[ii][198]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 3 jet events");
h_[ii][198]->GetYaxis()->SetTitle("Entries");

h2_[ii][199] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_3j_1_1");
h_[ii][199]= (TH1F*)h2_[ii][199]->Clone();
h_[ii][199]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 3 central jet events");
h_[ii][199]->GetYaxis()->SetTitle("Entries");

h2_[ii][200] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_3j_1_1");
h_[ii][200]= (TH1F*)h2_[ii][200]->Clone();
h_[ii][200]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 3 central jet events");
h_[ii][200]->GetYaxis()->SetTitle("Entries");

h2_[ii][201] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_3j_1_1");
h_[ii][201]= (TH1F*)h2_[ii][201]->Clone();
h_[ii][201]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 3 forward jet events");
h_[ii][201]->GetYaxis()->SetTitle("Entries");

h2_[ii][202] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_3j_1_1");
h_[ii][202]= (TH1F*)h2_[ii][202]->Clone();
h_[ii][202]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 3 forward jet events");
h_[ii][202]->GetYaxis()->SetTitle("Entries");

h2_[ii][203] =(TH1F*) file[ii] ->Get("h_mz_4j_all_1_1");
h_[ii][203]= (TH1F*)h2_[ii][203]->Clone();
h_[ii][203]->GetXaxis()->SetTitle("Z mass for Z + 4 jet events");
h_[ii][203]->GetYaxis()->SetTitle("Entries");

h2_[ii][204] =(TH1F*) file[ii] ->Get("h_mz_4j_c_1_1");
h_[ii][204]= (TH1F*)h2_[ii][204]->Clone();
h_[ii][204]->GetXaxis()->SetTitle("Z mass for Z + 4 central jet events");
h_[ii][204]->GetYaxis()->SetTitle("Entries");

h2_[ii][205] =(TH1F*) file[ii] ->Get("h_mz_4j_hf_1_1");
h_[ii][205]= (TH1F*)h2_[ii][205]->Clone();
h_[ii][205]->GetXaxis()->SetTitle("Z mass for Z + 4 forward jet events");
h_[ii][205]->GetYaxis()->SetTitle("Entries");

h2_[ii][206] =(TH1F*) file[ii] ->Get("h_dielec_PT_4j_all_1_1");
h_[ii][206]= (TH1F*)h2_[ii][206]->Clone();
h_[ii][206]->GetXaxis()->SetTitle("Z P_{T} for Z + 4 jet events");
h_[ii][206]->GetYaxis()->SetTitle("Entries");

h2_[ii][207] =(TH1F*) file[ii] ->Get("h_dielec_PT_4j_c_1_1");
h_[ii][207]= (TH1F*)h2_[ii][207]->Clone();
h_[ii][207]->GetXaxis()->SetTitle("Z P_{T} for Z + 4 central jet events");
h_[ii][207]->GetYaxis()->SetTitle("Entries");

h2_[ii][208] =(TH1F*) file[ii] ->Get("h_dielec_PT_4j_hf_1_1");
h_[ii][208]= (TH1F*)h2_[ii][208]->Clone();
h_[ii][208]->GetXaxis()->SetTitle("Z P_{T} for Z + 4 forward jet events");
h_[ii][208]->GetYaxis()->SetTitle("Entries");

h2_[ii][209] =(TH1F*) file[ii] ->Get("h_dielecphi_4j_all_1_1");
h_[ii][209]= (TH1F*)h2_[ii][209]->Clone();
h_[ii][209]->GetXaxis()->SetTitle("Z #phi for Z + 4 jet events");
h_[ii][209]->GetYaxis()->SetTitle("Entries");

h2_[ii][210] =(TH1F*) file[ii] ->Get("h_dielecphi_4j_c_1_1");
h_[ii][210]= (TH1F*)h2_[ii][210]->Clone();
h_[ii][210]->GetXaxis()->SetTitle("Z #phi for Z + 4 central jet events");
h_[ii][210]->GetYaxis()->SetTitle("Entries");

h2_[ii][211] =(TH1F*) file[ii] ->Get("h_dielecphi_4j_hf_1_1");
h_[ii][211]= (TH1F*)h2_[ii][211]->Clone();
h_[ii][211]->GetXaxis()->SetTitle("Z #phi for Z + 4 forward jet events");
h_[ii][211]->GetYaxis()->SetTitle("Entries");

h2_[ii][212] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_4j_all_1_1");
h_[ii][212]= (TH1F*)h2_[ii][212]->Clone();
h_[ii][212]->GetXaxis()->SetTitle("Z Y for Z + 4 jet events");
h_[ii][212]->GetYaxis()->SetTitle("Entries");

h2_[ii][213] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_4j_c_1_1");
h_[ii][213]= (TH1F*)h2_[ii][213]->Clone();
h_[ii][213]->GetXaxis()->SetTitle("Z Y for Z + 4 central jet events");
h_[ii][213]->GetYaxis()->SetTitle("Entries");

h2_[ii][214] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_4j_hf_1_1");
h_[ii][214]= (TH1F*)h2_[ii][214]->Clone();
h_[ii][214]->GetXaxis()->SetTitle("Z Y for Z + 4 forward jet events");
h_[ii][214]->GetYaxis()->SetTitle("Entries");

h2_[ii][215] =(TH1F*) file[ii] ->Get("h_jet_pt_all_4j_1_1");
h_[ii][215]= (TH1F*)h2_[ii][215]->Clone();
h_[ii][215]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 4 jet events");
h_[ii][215]->GetYaxis()->SetTitle("Entries");

h2_[ii][216] =(TH1F*) file[ii] ->Get("h_jet_eta_all_4j_1_1");
h_[ii][216]= (TH1F*)h2_[ii][216]->Clone();
h_[ii][216]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 4 jet events");
h_[ii][216]->GetYaxis()->SetTitle("Entries");

h2_[ii][217] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_4j_1_1");
h_[ii][217]= (TH1F*)h2_[ii][217]->Clone();
h_[ii][217]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 4 central jet events");
h_[ii][217]->GetYaxis()->SetTitle("Entries");

h2_[ii][218] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_4j_1_1");
h_[ii][218]= (TH1F*)h2_[ii][218]->Clone();
h_[ii][218]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 4 central jet events");
h_[ii][218]->GetYaxis()->SetTitle("Entries");

h2_[ii][219] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_4j_1_1");
h_[ii][219]= (TH1F*)h2_[ii][219]->Clone();
h_[ii][219]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 4 forward jet events");
h_[ii][219]->GetYaxis()->SetTitle("Entries");

h2_[ii][220] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_4j_1_1");
h_[ii][220]= (TH1F*)h2_[ii][220]->Clone();
h_[ii][220]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 4 forward jet events");
h_[ii][220]->GetYaxis()->SetTitle("Entries");

h2_[ii][221] =(TH1F*) file[ii] ->Get("h_mz_5j_all_1_1");
h_[ii][221]= (TH1F*)h2_[ii][221]->Clone();
h_[ii][221]->GetXaxis()->SetTitle("Z mass for Z + 5 jet events");
h_[ii][221]->GetYaxis()->SetTitle("Entries");

h2_[ii][222] =(TH1F*) file[ii] ->Get("h_mz_5j_c_1_1");
h_[ii][222]= (TH1F*)h2_[ii][222]->Clone();
h_[ii][222]->GetXaxis()->SetTitle("Z mass for Z + 5 central jet events");
h_[ii][222]->GetYaxis()->SetTitle("Entries");

h2_[ii][223] =(TH1F*) file[ii] ->Get("h_mz_5j_hf_1_1");
h_[ii][223]= (TH1F*)h2_[ii][223]->Clone();
h_[ii][223]->GetXaxis()->SetTitle("Z mass for Z + 5 forward jet events");
h_[ii][223]->GetYaxis()->SetTitle("Entries");

h2_[ii][224] =(TH1F*) file[ii] ->Get("h_dielec_PT_5j_all_1_1");
h_[ii][224]= (TH1F*)h2_[ii][224]->Clone();
h_[ii][224]->GetXaxis()->SetTitle("Z P_{T} for Z + 5 jet events");
h_[ii][224]->GetYaxis()->SetTitle("Entries");

h2_[ii][225] =(TH1F*) file[ii] ->Get("h_dielec_PT_5j_c_1_1");
h_[ii][225]= (TH1F*)h2_[ii][225]->Clone();
h_[ii][225]->GetXaxis()->SetTitle("Z P_{T} for Z + 5 central jet events");
h_[ii][225]->GetYaxis()->SetTitle("Entries");

h2_[ii][226] =(TH1F*) file[ii] ->Get("h_dielec_PT_5j_hf_1_1");
h_[ii][226]= (TH1F*)h2_[ii][226]->Clone();
h_[ii][226]->GetXaxis()->SetTitle("Z P_{T} for Z + 5 forward jet events");
h_[ii][226]->GetYaxis()->SetTitle("Entries");

h2_[ii][227] =(TH1F*) file[ii] ->Get("h_dielecphi_5j_all_1_1");
h_[ii][227]= (TH1F*)h2_[ii][227]->Clone();
h_[ii][227]->GetXaxis()->SetTitle("Z #phi for Z + 5 jet events");
h_[ii][227]->GetYaxis()->SetTitle("Entries");

h2_[ii][228] =(TH1F*) file[ii] ->Get("h_dielecphi_5j_c_1_1");
h_[ii][228]= (TH1F*)h2_[ii][228]->Clone();
h_[ii][228]->GetXaxis()->SetTitle("Z #phi for Z + 5 central jet events");
h_[ii][228]->GetYaxis()->SetTitle("Entries");

h2_[ii][229] =(TH1F*) file[ii] ->Get("h_dielecphi_5j_hf_1_1");
h_[ii][229]= (TH1F*)h2_[ii][229]->Clone();
h_[ii][229]->GetXaxis()->SetTitle("Z #phi for Z + 5 forward jet events");
h_[ii][229]->GetYaxis()->SetTitle("Entries");

h2_[ii][230] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_5j_all_1_1");
h_[ii][230]= (TH1F*)h2_[ii][230]->Clone();
h_[ii][230]->GetXaxis()->SetTitle("Z Y for Z + 5 jet events");
h_[ii][230]->GetYaxis()->SetTitle("Entries");

h2_[ii][231] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_5j_c_1_1");
h_[ii][231]= (TH1F*)h2_[ii][231]->Clone();
h_[ii][231]->GetXaxis()->SetTitle("Z Y for Z + 5 central jet events");
h_[ii][231]->GetYaxis()->SetTitle("Entries");

h2_[ii][232] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_5j_hf_1_1");
h_[ii][232]= (TH1F*)h2_[ii][232]->Clone();
h_[ii][232]->GetXaxis()->SetTitle("Z Y for Z + 5 forward jet events");
h_[ii][232]->GetYaxis()->SetTitle("Entries");

h2_[ii][233] =(TH1F*) file[ii] ->Get("h_jet_pt_all_5j_1_1");
h_[ii][233]= (TH1F*)h2_[ii][233]->Clone();
h_[ii][233]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 5 jet events");
h_[ii][233]->GetYaxis()->SetTitle("Entries");

h2_[ii][234] =(TH1F*) file[ii] ->Get("h_jet_eta_all_5j_1_1");
h_[ii][234]= (TH1F*)h2_[ii][234]->Clone();
h_[ii][234]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 5 jet events");
h_[ii][234]->GetYaxis()->SetTitle("Entries");

h2_[ii][235] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_5j_1_1");
h_[ii][235]= (TH1F*)h2_[ii][235]->Clone();
h_[ii][235]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 5 central jet events");
h_[ii][235]->GetYaxis()->SetTitle("Entries");

h2_[ii][236] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_5j_1_1");
h_[ii][236]= (TH1F*)h2_[ii][236]->Clone();
h_[ii][236]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 5 central jet events");
h_[ii][236]->GetYaxis()->SetTitle("Entries");

h2_[ii][237] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_5j_1_1");
h_[ii][237]= (TH1F*)h2_[ii][237]->Clone();
h_[ii][237]->GetXaxis()->SetTitle("Jet pt of leading jet in Z + 5 forward jet events");
h_[ii][237]->GetYaxis()->SetTitle("Entries");

h2_[ii][238] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_5j_1_1");
h_[ii][238]= (TH1F*)h2_[ii][238]->Clone();
h_[ii][238]->GetXaxis()->SetTitle("Jet #eta of leading jet in Z + 5 forward jet events");
h_[ii][238]->GetYaxis()->SetTitle("Entries");

h2_[ii][239] =(TH1F*) file[ii] ->Get("h_eEta_mva0_do02_pfiso01_mhits0_ss");
h_[ii][239]= (TH1F*)h2_[ii][239]->Clone();
h_[ii][239]->GetXaxis()->SetTitle("Electron Eta  ");
h_[ii][239]->GetYaxis()->SetTitle("Entries");

h2_[ii][240] =(TH1F*) file[ii] ->Get("h_ePhi_mva0_do02_pfiso01_mhits0_ss");
h_[ii][240]= (TH1F*)h2_[ii][240]->Clone();
h_[ii][240]->GetXaxis()->SetTitle("Electron Phi  ");
h_[ii][240]->GetYaxis()->SetTitle("Entries");

h2_[ii][241] =(TH1F*) file[ii] ->Get("h_ePt_mva0_do02_pfiso01_mhits0_ss");
h_[ii][241]= (TH1F*)h2_[ii][241]->Clone();
h_[ii][241]->GetXaxis()->SetTitle("Electron Pt  ");
h_[ii][241]->GetYaxis()->SetTitle("Entries");

h2_[ii][242] =(TH1F*) file[ii] ->Get("h_ePx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][242]= (TH1F*)h2_[ii][242]->Clone();
h_[ii][242]->GetXaxis()->SetTitle("Electron Px  ");
h_[ii][242]->GetYaxis()->SetTitle("Entries");

h2_[ii][243] =(TH1F*) file[ii] ->Get("h_ePy_mva0_do02_pfiso01_mhits0_ss");
h_[ii][243]= (TH1F*)h2_[ii][243]->Clone();
h_[ii][243]->GetXaxis()->SetTitle("Electron Py  ");
h_[ii][243]->GetYaxis()->SetTitle("Entries");

h2_[ii][244] =(TH1F*) file[ii] ->Get("h_ePz_mva0_do02_pfiso01_mhits0_ss");
h_[ii][244]= (TH1F*)h2_[ii][244]->Clone();
h_[ii][244]->GetXaxis()->SetTitle("Electron Pz  ");
h_[ii][244]->GetYaxis()->SetTitle("Entries");

h2_[ii][245] =(TH1F*) file[ii] ->Get("h_eM_mva0_do02_pfiso01_mhits0_ss");
h_[ii][245]= (TH1F*)h2_[ii][245]->Clone();
h_[ii][245]->GetXaxis()->SetTitle("Electron M  ");
h_[ii][245]->GetYaxis()->SetTitle("Entries");

h2_[ii][246] =(TH1F*) file[ii] ->Get("h_eE_mva0_do02_pfiso01_mhits0_ss");
h_[ii][246]= (TH1F*)h2_[ii][246]->Clone();
h_[ii][246]->GetXaxis()->SetTitle("Electron Energy  ");
h_[ii][246]->GetYaxis()->SetTitle("Entries");

h2_[ii][247] =(TH1F*) file[ii] ->Get("h_eMVATrigId_mva0_do02_pfiso01_mhits0_ss");
h_[ii][247]= (TH1F*)h2_[ii][247]->Clone();
h_[ii][247]->GetXaxis()->SetTitle("Electron MVA Trig Id  ");
h_[ii][247]->GetYaxis()->SetTitle("Entries");

h2_[ii][248] =(TH1F*) file[ii] ->Get("h_escSigmaIEtaIEta_mva0_do02_pfiso01_mhits0_ss");
h_[ii][248]= (TH1F*)h2_[ii][248]->Clone();
h_[ii][248]->GetXaxis()->SetTitle("Electron SuperCluster #sigma i #eta i #eta ");
h_[ii][248]->GetYaxis()->SetTitle("Entries");

h2_[ii][249] =(TH1F*) file[ii] ->Get("h_edeltaPhiSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][249]= (TH1F*)h2_[ii][249]->Clone();
h_[ii][249]->GetXaxis()->SetTitle("Electron SuperCluster #Delta #phi  ");
h_[ii][249]->GetYaxis()->SetTitle("Entries");

h2_[ii][250] =(TH1F*) file[ii] ->Get("h_edeltaEtaSuperClusterTrackAtVtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][250]= (TH1F*)h2_[ii][250]->Clone();
h_[ii][250]->GetXaxis()->SetTitle("Electron SuperCluster #Delta #eta  ");
h_[ii][250]->GetYaxis()->SetTitle("Entries");

h2_[ii][251] =(TH1F*) file[ii] ->Get("h_ehadronicOverEm_mva0_do02_pfiso01_mhits0_ss");
h_[ii][251]= (TH1F*)h2_[ii][251]->Clone();
h_[ii][251]->GetXaxis()->SetTitle("Electron E/M  ");
h_[ii][251]->GetYaxis()->SetTitle("Entries");

h2_[ii][252] =(TH1F*) file[ii] ->Get("h_egsfTrack_numberOfLostHits_mva0_do02_pfiso01_mhits0_ss");
h_[ii][252]= (TH1F*)h2_[ii][252]->Clone();
h_[ii][252]->GetXaxis()->SetTitle("Electron number of lost hits  ");
h_[ii][252]->GetYaxis()->SetTitle("Entries");

h2_[ii][253] =(TH1F*) file[ii] ->Get("h_ed0vtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][253]= (TH1F*)h2_[ii][253]->Clone();
h_[ii][253]->GetXaxis()->SetTitle("Electron d0  ");
h_[ii][253]->GetYaxis()->SetTitle("Entries");

h2_[ii][254] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][254]= (TH1F*)h2_[ii][254]->Clone();
h_[ii][254]->GetXaxis()->SetTitle("Electron dz  ");
h_[ii][254]->GetYaxis()->SetTitle("Entries");

h2_[ii][255] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][255]= (TH1F*)h2_[ii][255]->Clone();
h_[ii][255]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0_ss");
h_[ii][255]->GetYaxis()->SetTitle("Entries");

h2_[ii][256] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][256]= (TH1F*)h2_[ii][256]->Clone();
h_[ii][256]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0_ss");
h_[ii][256]->GetYaxis()->SetTitle("Entries");

h2_[ii][257] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][257]= (TH1F*)h2_[ii][257]->Clone();
h_[ii][257]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0_ss");
h_[ii][257]->GetYaxis()->SetTitle("Entries");

h2_[ii][258] =(TH1F*) file[ii] ->Get("h_edzvtx_mva0_do02_pfiso01_mhits0_ss");
h_[ii][258]= (TH1F*)h2_[ii][258]->Clone();
h_[ii][258]->GetXaxis()->SetTitle("OBSOLETE_mva0_do02_pfiso01_mhits0_ss");
h_[ii][258]->GetYaxis()->SetTitle("Entries");

h2_[ii][259] =(TH1F*) file[ii] ->Get("h_erelIso_mva0_do02_pfiso01_mhits0_ss");
h_[ii][259]= (TH1F*)h2_[ii][259]->Clone();
h_[ii][259]->GetXaxis()->SetTitle("Electron PF iso  ");
h_[ii][259]->GetYaxis()->SetTitle("Entries");

h2_[ii][260] =(TH1F*) file[ii] ->Get("h_erelIsodb_mva0_do02_pfiso01_mhits0_ss");
h_[ii][260]= (TH1F*)h2_[ii][260]->Clone();
h_[ii][260]->GetXaxis()->SetTitle("db corrected Electron PF iso  ");
h_[ii][260]->GetYaxis()->SetTitle("Entries");

h2_[ii][261] =(TH1F*) file[ii] ->Get("h_erelIsorho_mva0_do02_pfiso01_mhits0_ss");
h_[ii][261]= (TH1F*)h2_[ii][261]->Clone();
h_[ii][261]->GetXaxis()->SetTitle("EA corrected Electron PF iso  ");
h_[ii][261]->GetYaxis()->SetTitle("Entries");

h2_[ii][262] =(TH1F*) file[ii] ->Get("h_efMVAVar_fbrem_mva0_do02_pfiso01_mhits0_ss");
h_[ii][262]= (TH1F*)h2_[ii][262]->Clone();
h_[ii][262]->GetXaxis()->SetTitle("Electron fbrem(mva variable)  ");
h_[ii][262]->GetYaxis()->SetTitle("Entries");

h2_[ii][263] =(TH1F*) file[ii] ->Get("h_efMVAVar_kfchi2_mva0_do02_pfiso01_mhits0_ss");
h_[ii][263]= (TH1F*)h2_[ii][263]->Clone();
h_[ii][263]->GetXaxis()->SetTitle("Electron kfchi^2(mva variable)  ");
h_[ii][263]->GetYaxis()->SetTitle("Entries");

h2_[ii][264] =(TH1F*) file[ii] ->Get("h_efMVAVar_kfhits_mva0_do02_pfiso01_mhits0_ss");
h_[ii][264]= (TH1F*)h2_[ii][264]->Clone();
h_[ii][264]->GetXaxis()->SetTitle("Electron kfhits(mva variable)  ");
h_[ii][264]->GetYaxis()->SetTitle("Entries");

h2_[ii][265] =(TH1F*) file[ii] ->Get("h_efMVAVar_gsfchi2_mva0_do02_pfiso01_mhits0_ss");
h_[ii][265]= (TH1F*)h2_[ii][265]->Clone();
h_[ii][265]->GetXaxis()->SetTitle("Electron gsfchi^2(mva variable)  ");
h_[ii][265]->GetYaxis()->SetTitle("Entries");

h2_[ii][266] =(TH1F*) file[ii] ->Get("h_efMVAVar_detacalo_mva0_do02_pfiso01_mhits0_ss");
h_[ii][266]= (TH1F*)h2_[ii][266]->Clone();
h_[ii][266]->GetXaxis()->SetTitle("Electron d #eta Calo(mva variable)  ");
h_[ii][266]->GetYaxis()->SetTitle("Entries");

h2_[ii][267] =(TH1F*) file[ii] ->Get("h_efMVAVar_see_mva0_do02_pfiso01_mhits0_ss");
h_[ii][267]= (TH1F*)h2_[ii][267]->Clone();
h_[ii][267]->GetXaxis()->SetTitle("Electron  #sigma i #eta i #eta(mva variable)  ");
h_[ii][267]->GetYaxis()->SetTitle("Entries");

h2_[ii][268] =(TH1F*) file[ii] ->Get("h_efMVAVar_spp_mva0_do02_pfiso01_mhits0_ss");
h_[ii][268]= (TH1F*)h2_[ii][268]->Clone();
h_[ii][268]->GetXaxis()->SetTitle("Electron  #sigma i #phi i #phi(mva variable)  ");
h_[ii][268]->GetYaxis()->SetTitle("Entries");

h2_[ii][269] =(TH1F*) file[ii] ->Get("h_efMVAVar_etawidth_mva0_do02_pfiso01_mhits0_ss");
h_[ii][269]= (TH1F*)h2_[ii][269]->Clone();
h_[ii][269]->GetXaxis()->SetTitle("Electron eta width(mva variable)  ");
h_[ii][269]->GetYaxis()->SetTitle("Entries");

h2_[ii][270] =(TH1F*) file[ii] ->Get("h_efMVAVar_phiwidth_mva0_do02_pfiso01_mhits0_ss");
h_[ii][270]= (TH1F*)h2_[ii][270]->Clone();
h_[ii][270]->GetXaxis()->SetTitle("Electron phi width(mva variable)  ");
h_[ii][270]->GetYaxis()->SetTitle("Entries");

h2_[ii][271] =(TH1F*) file[ii] ->Get("h_efMVAVar_e1x5e5x5_mva0_do02_pfiso01_mhits0_ss");
h_[ii][271]= (TH1F*)h2_[ii][271]->Clone();
h_[ii][271]->GetXaxis()->SetTitle("Electron e1x5/e5x5(mva variable)  ");
h_[ii][271]->GetYaxis()->SetTitle("Entries");

h2_[ii][272] =(TH1F*) file[ii] ->Get("h_efMVAVar_R9_mva0_do02_pfiso01_mhits0_ss");
h_[ii][272]= (TH1F*)h2_[ii][272]->Clone();
h_[ii][272]->GetXaxis()->SetTitle("Second Leading Electron MVA_mva0_do02_pfiso01_mhits0_ss");
h_[ii][272]->GetYaxis()->SetTitle("Entries");

h2_[ii][273] =(TH1F*) file[ii] ->Get("h_efMVAVar_EoP_mva0_do02_pfiso01_mhits0_ss");
h_[ii][273]= (TH1F*)h2_[ii][273]->Clone();
h_[ii][273]->GetXaxis()->SetTitle("Electron scE/p(mva variable)  ");
h_[ii][273]->GetYaxis()->SetTitle("Entries");

h2_[ii][274] =(TH1F*) file[ii] ->Get("h_efMVAVar_IoEmIoP_mva0_do02_pfiso01_mhits0_ss");
h_[ii][274]= (TH1F*)h2_[ii][274]->Clone();
h_[ii][274]->GetXaxis()->SetTitle("Electron 1/E - 1/p(mva variable)  ");
h_[ii][274]->GetYaxis()->SetTitle("Entries");

h2_[ii][275] =(TH1F*) file[ii] ->Get("h_efMVAVar_eleEoPout_mva0_do02_pfiso01_mhits0_ss");
h_[ii][275]= (TH1F*)h2_[ii][275]->Clone();
h_[ii][275]->GetXaxis()->SetTitle("Electron sc E/P out(mva variable)  ");
h_[ii][275]->GetYaxis()->SetTitle("Entries");

h2_[ii][276] =(TH1F*) file[ii] ->Get("h_efMVAVar_PreShowerOverRaw_mva0_do02_pfiso01_mhits0_ss");
h_[ii][276]= (TH1F*)h2_[ii][276]->Clone();
h_[ii][276]->GetXaxis()->SetTitle("Electron PreShowerOverRaw(mva variable)  ");
h_[ii][276]->GetYaxis()->SetTitle("Entries");



}

int number_of_histos = 277;
char name1[500];
char name2[500];
// TH1F *hsum[100];
TCanvas *c[500];
TH1F*qcd_bckg_0;
TH1F*qcd_bckg_1;
TH1F*qcd_bckg_2;
TH1F*qcd_bckg_3;
TH1F *Mont[500];

int ss_index[500];

ss_index[0]=143;

ss_index[1]=144;

ss_index[2]=145;

ss_index[3]=22;

ss_index[7]=142;
for(int j==0;j<numberoffiles;j++){
for(int ii = 0; ii<number_of_histos;ii++){
if(ii>7 && ii<29)ss_index[ii]=ii+138;
if(ii>28&&ii<67)ss_index[ii]=ii + 210;
if(ii>69 && ii<142)ss_index[ii]=ii+97;
int jj= ss_index[ii];
cout<<ii<<"  "<<ss_index[ii]<<endl;
if(ii==138||ii==136||ii==128||ii==127||ii==120||ii==102||ii==118||ii==110||ii==109||ii==92||ii==91||ii==75){
h_[j][ii]->Rebin(2);
h_[j][jj]->Rebin(2);
}}}

for(int ii = 0; ii<number_of_histos;ii++){
int jj= ss_index[ii];
if(ii>142)continue;
//if(ii>7)continue;
//if(ii>37&&ii<68)continue;
//if(ii==20) continue;
//if(ii!=22 &&ii!=68)continue;
//if(!(ii==58||ii==57||(ii<49 && ii>44)||ii==40||(ii<29&&ii>18 &&ii!=22)||ii==11))continue;
//if(ii!=11)continue;
bool logx=false;
bool logy=false;
bool integral_norm = false;
bool same_sign_0 = false;
bool same_sign_1 = false;
bool same_sign_2 = false;
bool same_sign_add = false;
if(ii==0 ||ii==1||ii==2||ii==3||(ii>6&&ii<29)||(ii>69&&ii<142))same_sign_add=true;
//cout<<jj<<endl;
if((ii>58&&ii<68&& ii!=62)||ii==56||ii==55||ii==54||ii==53||ii==51||ii==50||ii==49 ||(ii>30&&ii<45&& ii!=35)||ii==22||ii==29||(ii<19&&ii>=0&&(ii!=7||ii!=11||ii!=15)) ||(ii<29&&ii>22) ||ii==125||ii==124||ii==107||ii==106||ii==89||ii==88||ii==72||ii==71||ii==70)logy=true;

if(ii==66||ii==65||ii==63||ii==60||ii==59||ii==55||ii==51||ii==50||ii==49||ii==41)logx=true;

//if(ii==68) same_sign_0=true;
//if(ii==69) same_sign_2=true;
//if(ii==22) same_sign_1=true;


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
if(same_sign_add)mc_bg->Add(h_[0][jj]);
//		mc_bg->Add(diboson);
//		mc_bg->Add(h_[2][ii]);
//		mc_bg->Add(h_[3][ii]);
//		mc_bg->Add(h_[4][ii]);
//		mc_bg->Add(h_[8][ii]);
}

if(integral_norm)TH1F*mc_bg = (TH1F*)h_[1][ii]->Clone();

Mont[ii]= (TH1F*)h_[0][ii]->Clone();

Mont[ii]->Divide(h_[0][ii],mc_bg,1.0,1.0);



sprintf(name1,"histo_z_plot_new_52X/hist_%i.root",ii);
sprintf(name2,"histo_z_plot_new_52X/hist_%i.pdf",ii);
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
if(same_sign_add){
h_[0][jj]->Draw("hhistsames");
h_[0][jj]->SetTitle("Same sign");
h_[0][jj].SetLineStyle(4);
h_[0][jj].SetLineColor(42);
h_[0][jj].SetLineWidth(2);
}
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
if(same_sign_0){
qcd_bckg_0 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_0->Add(h_[1][ii],-1);
}

if(same_sign_1){
qcd_bckg_1 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_1->Add(h_[1][ii],-1);
}
if(same_sign_2){
qcd_bckg_2 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_2->Add(h_[1][ii],-1);
}

c[ii]->cd(2);
  gStyle->SetOptStat(kFALSE);
          Mont[ii]->SetMinimum(0.5);
	  Mont[ii]->SetMaximum(1.5);
          Mont[ii]->SetYTitle("Ratio");
	  Mont[ii]->Draw("e1");
if(logx)gPad->SetLogx();
c[ii]->Print(name1);
c[ii]->Print(name2);

}
/*
TCanvas *c_qcd_0 = new TCanvas();
c_qcd_0->cd();
qcd_bckg_0->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG) with electrons PFIso>0.1");
qcd_bckg_0->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_0->GetYaxis()->SetTitle("Entries");
qcd_bckg_0->Draw("E1");
c_qcd_0->Print("histo_z_plot_new/qcd_bckground_0.root");
c_qcd_0->Print("histo_z_plot_new/qcd_bckground_0.pdf");

TCanvas *c_qcd_1 = new TCanvas();
c_qcd_1->cd();
qcd_bckg_1->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG)with electrons PFIso<0.1");
qcd_bckg_1->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_1->GetYaxis()->SetTitle("Entries");
qcd_bckg_1->Draw("E1");
c_qcd_1->Print("histo_z_plot_new/qcd_bckground_1.root");
c_qcd_1->Print("histo_z_plot_new/qcd_bckground_1.pdf");

TCanvas *c_qcd_2 = new TCanvas();
c_qcd_2->cd();
qcd_bckg_2->GetXaxis()->SetTitle("opposite-sign ee inv.M. Data-(MC+BG)with electrons PFIso>0.1");
qcd_bckg_2->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_2->GetYaxis()->SetTitle("Entries");
qcd_bckg_2->Draw("E1");
c_qcd_2->Print("histo_z_plot_new/qcd_bckground_2.root");
c_qcd_2->Print("histo_z_plot_new/qcd_bckground_2.pdf");
*/
}
