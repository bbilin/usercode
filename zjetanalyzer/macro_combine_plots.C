{

 //  gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

TFile * file[500];

file[0] = TFile::Open("rootfiles/29_08_data.root");

file[1] = TFile::Open("rootfiles/29_08_mc.root");

file[2] = TFile::Open("rootfiles/29_08_bg_ztautau.root");

file[3] = TFile::Open("rootfiles/29_08_bg_tt.root");

file[4] = TFile::Open("rootfiles/29_08_bg_wjets.root");

file[5] = TFile::Open("rootfiles/29_08_ww.root");
file[6] = TFile::Open("rootfiles/29_08_wz.root");
file[7] = TFile::Open("rootfiles/29_08_zz.root");

file[8] = TFile::Open("rootfiles/29_08_mc_low.root");

TH1F *h_[500][500];
TH1F *h1_[500][500];
int numberoffiles = 9;
int number_of_histos = 79;

for(int ii=0;ii<numberoffiles;ii++){
h_[ii][0] =(TH1F*) file[ii] ->Get("h_elec_mva_l_mva_with_transition_cut");
h_[ii][0]->GetXaxis()->SetTitle("Leading Electron MVA + transition cut");
h_[ii][0]->GetYaxis()->SetTitle("Entries");
h_[ii][1] =(TH1F*) file[ii] ->Get("h_elec_mva_nl_mva_with_transition_cut");
h_[ii][1]->GetXaxis()->SetTitle("Second Leading Electron MVA + transition cut");
h_[ii][1]->GetYaxis()->SetTitle("Entries");
h_[ii][2] =(TH1F*) file[ii] ->Get("h_hfeL1");
h_[ii][2]->GetXaxis()->SetTitle("HF Electron L1");
h_[ii][2]->GetYaxis()->SetTitle("Entries");
h_[ii][3] =(TH1F*) file[ii] ->Get("h_hfeL9");
h_[ii][3]->GetXaxis()->SetTitle("HF Electron L9");
h_[ii][3]->GetYaxis()->SetTitle("Entries");
h_[ii][4] =(TH1F*) file[ii] ->Get("h_hfeS9");
h_[ii][4]->GetXaxis()->SetTitle("HF Electron S9");
h_[ii][4]->GetYaxis()->SetTitle("Entries");
h_[ii][5] =(TH1F*) file[ii] ->Get("h_hfeL25");
h_[ii][5]->GetXaxis()->SetTitle("HF Electron L25");
h_[ii][5]->GetYaxis()->SetTitle("Entries");
h_[ii][6] =(TH1F*) file[ii] ->Get("h_hfecut1");
h_[ii][6]->GetXaxis()->SetTitle("HF Electron L9/L25");
h_[ii][6]->GetYaxis()->SetTitle("Entries");
h_[ii][7] =(TH1F*) file[ii] ->Get("h_hfecut2");
h_[ii][7]->GetXaxis()->SetTitle("HF Electron L1/L9 - 1.125#times S9/L9");
h_[ii][7]->GetYaxis()->SetTitle("Entries");
h_[ii][8] =(TH1F*) file[ii] ->Get("h_hfe_Pt");
h_[ii][8]->GetXaxis()->SetTitle("Leading HF Electron P_{T}");
h_[ii][8]->GetYaxis()->SetTitle("Entries");
h_[ii][9] =(TH1F*) file[ii] ->Get("h_hfe_Eta");
h_[ii][9]->GetXaxis()->SetTitle("Leading HF Electron #eta");
h_[ii][9]->GetYaxis()->SetTitle("Entries");
h_[ii][10] =(TH1F*) file[ii] ->Get("h_no_good_vtx");
h_[ii][10]->GetXaxis()->SetTitle("# of good vertices");
h_[ii][10]->GetYaxis()->SetTitle("Entries");


h_[ii][11] =(TH1F*) file[ii] ->Get("h_elec_pt_0_1");
h_[ii][11]->GetXaxis()->SetTitle("Leading Electron P_{T} of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][11]->GetYaxis()->SetTitle("Entries");
h_[ii][12] =(TH1F*) file[ii] ->Get("h_elec_eta_0_1");
h_[ii][12]->GetXaxis()->SetTitle("Leading Electron #eta of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][12]->GetYaxis()->SetTitle("Entries");
h_[ii][13] =(TH1F*) file[ii] ->Get("h_numberofPFjets_0_1");
h_[ii][13]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][13]->GetYaxis()->SetTitle("Entries");
h_[ii][14] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_0_1");
h_[ii][14]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][14]->GetYaxis()->SetTitle("Entries");
h_[ii][15] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_0_1");
h_[ii][15]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][15]->GetYaxis()->SetTitle("Entries");
h_[ii][16] =(TH1F*) file[ii] ->Get("h_mz_0_1");
h_[ii][16]->GetXaxis()->SetTitle("Z Mass of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][16]->GetYaxis()->SetTitle("Entries");
h_[ii][17] =(TH1F*) file[ii] ->Get("h_mz_1j_all_0_1");
h_[ii][17]->GetXaxis()->SetTitle("Z mass for excl. Z+1 jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][17]->GetYaxis()->SetTitle("Entries");
h_[ii][18] =(TH1F*) file[ii] ->Get("h_mz_1j_c_0_1");
h_[ii][18]->GetXaxis()->SetTitle("Z mass for excl. Z+1 central jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][18]->GetYaxis()->SetTitle("Entries");
h_[ii][19] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_0_1");
h_[ii][19]->GetXaxis()->SetTitle("Z mass for excl. Z+1 HF jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][19]->GetYaxis()->SetTitle("Entries");
h_[ii][20] =(TH1F*) file[ii] ->Get("h_dielec_PT_0_1");
h_[ii][20]->GetXaxis()->SetTitle("Z P_{T}  of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][20]->GetYaxis()->SetTitle("Entries");
h_[ii][21] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_0_1");
h_[ii][21]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][21]->GetYaxis()->SetTitle("Entries");
h_[ii][22] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_0_1");
h_[ii][22]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][22]->GetYaxis()->SetTitle("Entries");
h_[ii][23] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_0_1");
h_[ii][23]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][23]->GetYaxis()->SetTitle("Entries");
h_[ii][24] =(TH1F*) file[ii] ->Get("h_dielecphi_0_1");
h_[ii][24]->GetXaxis()->SetTitle("Z #phi of MVA(e1)>0.0 MVA(e2)>0.0 ");
h_[ii][24]->GetYaxis()->SetTitle("Entries");
h_[ii][25] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_0_1");
h_[ii][25]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][25]->GetYaxis()->SetTitle("Entries");
h_[ii][26] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_0_1");
h_[ii][26]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][26]->GetYaxis()->SetTitle("Entries");
h_[ii][27] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_0_1");
h_[ii][27]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][27]->GetYaxis()->SetTitle("Entries");
h_[ii][28] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_0_1");
h_[ii][28]->GetXaxis()->SetTitle("Z Y of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][28]->GetYaxis()->SetTitle("Entries");
h_[ii][29] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_0_1");
h_[ii][29]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][29]->GetYaxis()->SetTitle("Entries");
h_[ii][30] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_0_1");
h_[ii][30]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][30]->GetYaxis()->SetTitle("Entries");

h_[ii][31] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_0_1");
h_[ii][31]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][31]->GetYaxis()->SetTitle("Entries");

h_[ii][32] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_0_1");
h_[ii][32]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][32]->GetYaxis()->SetTitle("Entries");

h_[ii][33] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_hf_0_1");
h_[ii][33]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][33]->GetYaxis()->SetTitle("Entries");

h_[ii][34] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_0_1");
h_[ii][34]->GetXaxis()->SetTitle("# of good vertices of Z events of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][34]->GetYaxis()->SetTitle("Entries");

h_[ii][35] =(TH1F*) file[ii] ->Get("h_mz_same_sign_0_1");
h_[ii][35]->GetXaxis()->SetTitle("same sign electron invariant mass of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][35]->GetYaxis()->SetTitle("Entries");

h_[ii][36] =(TH1F*) file[ii] ->Get("h_jet_pt_all_0_1");
h_[ii][36]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][36]->GetYaxis()->SetTitle("Entries");

h_[ii][37] =(TH1F*) file[ii] ->Get("h_jet_eta_all_0_1");
h_[ii][37]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][37]->GetYaxis()->SetTitle("Entries");

h_[ii][38] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_0_1");
h_[ii][38]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][38]->GetYaxis()->SetTitle("Entries");

h_[ii][39] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_0_1");
h_[ii][39]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][39]->GetYaxis()->SetTitle("Entries");

h_[ii][40] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_0_1");
h_[ii][40]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][40]->GetYaxis()->SetTitle("Entries");

h_[ii][41] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_0_1");
h_[ii][41]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  of MVA(e1)>0.0 MVA(e2)>0.0");
h_[ii][41]->GetYaxis()->SetTitle("Entries");

h_[ii][42] =(TH1F*) file[ii] ->Get("h_elec_pt_1_0");
h_[ii][42]->GetXaxis()->SetTitle("Leading Electron P_{T} of cut based electron pairs");
h_[ii][42]->GetYaxis()->SetTitle("Entries");

h_[ii][43] =(TH1F*) file[ii] ->Get("h_elec_eta_1_0");
h_[ii][43]->GetXaxis()->SetTitle("Leading Electron #eta of cut based electron pairs");
h_[ii][43]->GetYaxis()->SetTitle("Entries");

h_[ii][44] =(TH1F*) file[ii] ->Get("h_numberofPFjets_1_0");
h_[ii][44]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity of cut based electron pairs");
h_[ii][44]->GetYaxis()->SetTitle("Entries");

h_[ii][45] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_1_0");
h_[ii][45]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE of cut based electron pairs");
h_[ii][45]->GetYaxis()->SetTitle("Entries");

h_[ii][46] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_1_0");
h_[ii][46]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF of cut based electron pairs");
h_[ii][46]->GetYaxis()->SetTitle("Entries");

h_[ii][47] =(TH1F*) file[ii] ->Get("h_mz_1_0");
h_[ii][47]->GetXaxis()->SetTitle("Z Mass of cut based electron pairs");
h_[ii][47]->GetYaxis()->SetTitle("Entries");

h_[ii][48] =(TH1F*) file[ii] ->Get("h_mz_1j_all_1_0");
h_[ii][48]->GetXaxis()->SetTitle("Z mass for excl. Z+1 jet events of cut based electron pairs");
h_[ii][48]->GetYaxis()->SetTitle("Entries");

h_[ii][49] =(TH1F*) file[ii] ->Get("h_mz_1j_c_1_0");
h_[ii][49]->GetXaxis()->SetTitle("Z mass for excl. Z+1 central jet events of cut based electron pairs");
h_[ii][49]->GetYaxis()->SetTitle("Entries");

h_[ii][50] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_1_0");
h_[ii][50]->GetXaxis()->SetTitle("Z mass for excl. Z+1 HF jet events of cut based electron pairs");
h_[ii][50]->GetYaxis()->SetTitle("Entries");

h_[ii][51] =(TH1F*) file[ii] ->Get("h_dielec_PT_1_0");
h_[ii][51]->GetXaxis()->SetTitle("Z P_{T}  of cut based electron pairs");
h_[ii][51]->GetYaxis()->SetTitle("Entries");

h_[ii][52] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_1_0");
h_[ii][52]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events of cut based electron pairs");
h_[ii][52]->GetYaxis()->SetTitle("Entries");

h_[ii][53] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_1_0");
h_[ii][53]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events of cut based electron pairs");
h_[ii][53]->GetYaxis()->SetTitle("Entries");

h_[ii][54] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_1_0");
h_[ii][54]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events of cut based electron pairs");
h_[ii][54]->GetYaxis()->SetTitle("Entries");

h_[ii][55] =(TH1F*) file[ii] ->Get("h_dielecphi_1_0");
h_[ii][55]->GetXaxis()->SetTitle("Z #phi of cut based electron pairs ");
h_[ii][55]->GetYaxis()->SetTitle("Entries");

h_[ii][56] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_1_0");
h_[ii][56]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events of cut based electron pairs");
h_[ii][56]->GetYaxis()->SetTitle("Entries");

h_[ii][57] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_1_0");
h_[ii][57]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events of cut based electron pairs");
h_[ii][57]->GetYaxis()->SetTitle("Entries");

h_[ii][58] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_1_0");
h_[ii][58]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events of cut based electron pairs");
h_[ii][58]->GetYaxis()->SetTitle("Entries");

h_[ii][59] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1_0");
h_[ii][59]->GetXaxis()->SetTitle("Z Y of cut based electron pairs");
h_[ii][59]->GetYaxis()->SetTitle("Entries");

h_[ii][60] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_1_0");
h_[ii][60]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events of cut based electron pairs");
h_[ii][60]->GetYaxis()->SetTitle("Entries");

h_[ii][61] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_1_0");
h_[ii][61]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events of cut based electron pairs");
h_[ii][61]->GetYaxis()->SetTitle("Entries");

h_[ii][62] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_1_0");
h_[ii][62]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events of cut based electron pairs");
h_[ii][62]->GetYaxis()->SetTitle("Entries");

h_[ii][63] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_0");
h_[ii][63]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets of cut based electron pairs");
h_[ii][63]->GetYaxis()->SetTitle("Entries");

h_[ii][64] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_hf_1_0");
h_[ii][64]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets of cut based electron pairs");
h_[ii][64]->GetYaxis()->SetTitle("Entries");

h_[ii][65] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_1_0");
h_[ii][65]->GetXaxis()->SetTitle("# of good vertices of Z events of cut based electron pairs");
h_[ii][65]->GetYaxis()->SetTitle("Entries");

h_[ii][66] =(TH1F*) file[ii] ->Get("h_mz_same_sign_1_0");
h_[ii][66]->GetXaxis()->SetTitle("same sign electron invariant mass of cut based electron pairs");
h_[ii][66]->GetYaxis()->SetTitle("Entries");

h_[ii][67] =(TH1F*) file[ii] ->Get("h_jet_pt_all_1_0");
h_[ii][67]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector of cut based electron pairs");
h_[ii][67]->GetYaxis()->SetTitle("Entries");

h_[ii][68] =(TH1F*) file[ii] ->Get("h_jet_eta_all_1_0");
h_[ii][68]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector of cut based electron pairs");
h_[ii][68]->GetYaxis()->SetTitle("Entries");

h_[ii][69] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_1_0");
h_[ii][69]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  of cut based electron pairs");
h_[ii][69]->GetYaxis()->SetTitle("Entries");

h_[ii][70] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_1_0");
h_[ii][70]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  of cut based electron pairs");
h_[ii][70]->GetYaxis()->SetTitle("Entries");

h_[ii][71] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_1_0");
h_[ii][71]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  of cut based electron pairs");
h_[ii][71]->GetYaxis()->SetTitle("Entries");

h_[ii][72] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_1_0");
h_[ii][72]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  of cut based electron pairs");
h_[ii][72]->GetYaxis()->SetTitle("Entries");

h_[ii][73] =(TH1F*) file[ii] ->Get("h_elec_mva_l_cut");
h_[ii][73]->GetXaxis()->SetTitle("Leading Electron MVA of cut based electron");
h_[ii][73]->GetYaxis()->SetTitle("Entries");

h_[ii][74] =(TH1F*) file[ii] ->Get("h_elec_mva_nl_cut");
h_[ii][74]->GetXaxis()->SetTitle("Second Leading Electron MVA of cut based electron");
h_[ii][74]->GetYaxis()->SetTitle("Entries");

h_[ii][75] =(TH1F*) file[ii] ->Get("h_elec_mva");
h_[ii][75]->GetXaxis()->SetTitle("Electron MVA (normalized to area)");
h_[ii][75]->GetYaxis()->SetTitle("Entries");

h_[ii][76] =(TH1F*) file[ii] ->Get("h_elec_mva_pt25");
h_[ii][76]->GetXaxis()->SetTitle("Electron MVA for P_{T}(e)>25 GeV. (normalized to area)");
h_[ii][76]->GetYaxis()->SetTitle("Entries");

h_[ii][77] =(TH1F*) file[ii] ->Get("h_elec_mva_l_mva_no_transition_cut");
h_[ii][77]->GetXaxis()->SetTitle("Leading Electron MVA ");
h_[ii][77]->GetYaxis()->SetTitle("Entries");

h_[ii][78] =(TH1F*) file[ii] ->Get("h_elec_mva_nl_mva_no_transition_cut");
h_[ii][78]->GetXaxis()->SetTitle("Second Leading Electron MVA ");
h_[ii][78]->GetYaxis()->SetTitle("Entries");

h1_[ii][11] =(TH1F*) file[ii] ->Get("h_elec_pt_1_1");
h1_[ii][12] =(TH1F*) file[ii] ->Get("h_elec_eta_1_1");
h1_[ii][13] =(TH1F*) file[ii] ->Get("h_numberofPFjets_1_1");
h1_[ii][14] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_1_1");
h1_[ii][15] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_1_1");
h1_[ii][16] =(TH1F*) file[ii] ->Get("h_mz_1_1");
h1_[ii][17] =(TH1F*) file[ii] ->Get("h_mz_1j_all_1_1");
h1_[ii][18] =(TH1F*) file[ii] ->Get("h_mz_1j_c_1_1");
h1_[ii][19] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_1_1");
h1_[ii][20] =(TH1F*) file[ii] ->Get("h_dielec_PT_1_1");
h1_[ii][21] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_1_1");
h1_[ii][22] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_1_1");
h1_[ii][23] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_1_1");
h1_[ii][24] =(TH1F*) file[ii] ->Get("h_dielecphi_1_1");
h1_[ii][25] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_1_1");
h1_[ii][26] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_1_1");
h1_[ii][27] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_1_1");
h1_[ii][28] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1_1");
h1_[ii][29] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_1_1");
h1_[ii][30] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_1_1");
h1_[ii][31] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_1_1");
h1_[ii][32] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_1_1");
h1_[ii][33] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_hf_1_1");
h1_[ii][34] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_1_1");
h1_[ii][35] =(TH1F*) file[ii] ->Get("h_mz_same_sign_1_1");
h1_[ii][36] =(TH1F*) file[ii] ->Get("h_jet_pt_all_1_1");
h1_[ii][37] =(TH1F*) file[ii] ->Get("h_jet_eta_all_1_1");
h1_[ii][38] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_1_1");
h1_[ii][39] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe_1_1");
h1_[ii][40] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_1_1");
h1_[ii][41] =(TH1F*) file[ii] ->Get("h_jet_eta_hf_1_1");


/*
for(int jj=0;jj<number_of_histos;jj++){
h_[ii][jj]->GetXaxis()->SetTitleSize(0.028);
if(jj>10 && jj<42)h_[ii][jj]->Add(h1_[ii][jj]);
if(jj>41 && jj<73)h_[ii][jj]->Add(h1_[ii][jj-31]);
}*/
}

char name1[500];
char name2[500];
// TH1F *hsum[100];
TCanvas *c[500];
TH1F*qcd_bckg_0;
TH1F*qcd_bckg_1;
TH1F*qcd_bckg_4;
TH1F*qcd_bckg_5;
for(int ii = 0; ii<number_of_histos;ii++){

//if(!(ii==76 ||ii==75)) continue;

bool integral_norm = false;
bool same_sign_0 = false;
bool same_sign_1 = false;
bool same_sign_4 = false;
bool same_sign_5 = false;
if((ii>1 && ii<11)|| ii==106 || ii==107 )integral_norm = true;
if(ii==35) same_sign_0=true; if(ii==66) same_sign_1=true;  

cout<<ii<<same_sign_0<<endl;
c[ii] = new TCanvas();
//c[ii]->Divide(1,2,0.01,0);


h_[0][ii]->Sumw2();
if(!integral_norm){
h_[1][ii]->Sumw2();
h_[1][ii]->Scale((3503.75*5097.)/27759468.);
h_[2][ii]->Sumw2();
h_[2][ii]->Scale((3503.75*5097.)/27759468.);
h_[3][ii]->Sumw2();
h_[3][ii]->Scale((225.197*5097.)/6642654.);
h_[4][ii]->Sumw2();
h_[4][ii]->Scale((36257.2*5097.)/16487626.);
h_[5][ii]->Sumw2();
h_[5][ii]->Scale((54.838*5097.)/10000431.);
h_[6][ii]->Sumw2();
h_[6][ii]->Scale((32.31*5097.)/9996622.);
h_[7][ii]->Sumw2();
h_[7][ii]->Scale((17.654*5097.)/9799908.);

h_[8][ii]->Scale((1000.25*5097.)/6665597.);
}
if(integral_norm){
h_[1][ii]->Scale(h_[0][ii]->Integral()/h_[1][ii]->Integral());
}
if(!integral_norm){
TH1F *diboson = (TH1F*)h_[5][ii]->Clone("diboson");
		diboson ->Add(h_[6][ii]);
		diboson ->Add(h_[7][ii]);

TH1F*mc_bg = (TH1F*)diboson->Clone("mc_bg");
		mc_bg->Add(h_[1][ii]);
		mc_bg->Add(h_[2][ii]);
		mc_bg->Add(h_[3][ii]);
		mc_bg->Add(h_[4][ii]);
		mc_bg->Add(h_[8][ii]);
}

if(integral_norm)TH1F*mc_bg = (TH1F*)h_[1][ii]->Clone("mc_bg");

TH1F *Mont= (TH1F*)h_[0][ii]->Clone("Mont");

Mont->Divide(h_[0][ii],mc_bg,1.0,1.0);



sprintf(name1,"histo_z_plot/hist_%i.root",ii);
sprintf(name2,"histo_z_plot/hist_%i.pdf",ii);
          c[ii]->cd(1); 

mc_bg->Draw("hhist");
mc_bg->SetTitle("Z/#gamma* -> ee + BG");
mc_bg.SetLineColor(4);
mc_bg.SetLineWidth(2);

if(!integral_norm){

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

}
h_[0][ii]->Draw("e1sames");
h_[0][ii]->SetTitle("CMS Data 5.09 fb^{-1}");
h_[0][ii]->SetMarkerStyle(2);
//h_[0][ii]->SetLineWidth(3); 

if(same_sign_0){
qcd_bckg_0 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_0->Add(h_[1][ii],-1);
}
if(same_sign_1){
qcd_bckg_1 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_1->Add(h_[1][ii],-1);
}
//c[ii]->cd(2);
//  gStyle->SetOptStat(kFALSE);
//          Mont->SetMinimum(0.5);
//	  Mont->SetMaximum(1.5);
//          Mont->SetYTitle("Ratio");
//	  Mont->Draw("e1");


//TCanvas->Return(c);
c[ii]->Print(name1);
c[ii]->Print(name2);

}
/*
TCanvas *c_qcd_0 = new TCanvas();
TCanvas *c_qcd_1 = new TCanvas();
c_qcd_0->cd();
qcd_bckg_0->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG) of MVA(e1)>0.0 MVA(e2)>0.0");
qcd_bckg_0->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_0->GetYaxis()->SetTitle("Entries");
qcd_bckg_0->Draw("E1");
c_qcd_0->Print("histo_z_plot/qcd_bckground_0.root");
c_qcd_0->Print("histo_z_plot/qcd_bckground_0.pdf");
c_qcd_1->cd();
qcd_bckg_1->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG) of cut based electron pairs");
qcd_bckg_1->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_1->GetYaxis()->SetTitle("Entries");
qcd_bckg_1->Draw("E1");
c_qcd_1->Print("histo_z_plot/qcd_bckground_1.root");
c_qcd_1->Print("histo_z_plot/qcd_bckground_1.pdf");
*/


}
