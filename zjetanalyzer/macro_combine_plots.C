{
//gROOT->SetStyle("tdrStyle");



//  gROOT->SetBatch(1);
//gStyle->SetHistFillColor(kBlue);

//TFile * file[0] = TFile::Open("rootfiles/elec_data_bckup_nozee_test.root");



TFile * file[500];

file[0] = TFile::Open("rootfiles/16_08_data.root");

file[1] = TFile::Open("rootfiles/16_08_mc.root");

file[2] = TFile::Open("rootfiles/16_08_bg_ztautau.root");

file[3] = TFile::Open("rootfiles/16_08_bg_tt.root");

file[4] = TFile::Open("rootfiles/16_08_bg_wjets.root");

file[5] = TFile::Open("rootfiles/16_08_ww.root");
file[6] = TFile::Open("rootfiles/16_08_wz.root");
file[7] = TFile::Open("rootfiles/16_08_zz.root");

file[8] = TFile::Open("rootfiles/16_08_mc_low.root");

TH1F *h_[500][500];
int numberoffiles = 9;
int number_of_histos = 108;

for(int ii=0;ii<numberoffiles;ii++){
h_[ii][0] =(TH1F*) file[ii] ->Get("h_elec_mva_l");
h_[ii][0]->GetXaxis()->SetTitle("Leading Electron MVA");
h_[ii][0]->GetYaxis()->SetTitle("Entries");
h_[ii][1] =(TH1F*) file[ii] ->Get("h_elec_mva_nl_mva");
h_[ii][1]->GetXaxis()->SetTitle("Second Leading Electron MVA");
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


h_[ii][11] =(TH1F*) file[ii] ->Get("h_elec_pt_0123");
h_[ii][11]->GetXaxis()->SetTitle("Leading Electron P_{T} of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][11]->GetYaxis()->SetTitle("Entries");
h_[ii][12] =(TH1F*) file[ii] ->Get("h_elec_eta_0123");
h_[ii][12]->GetXaxis()->SetTitle("Leading Electron #eta of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][12]->GetYaxis()->SetTitle("Entries");
h_[ii][13] =(TH1F*) file[ii] ->Get("h_numberofPFjets_0123");
h_[ii][13]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][13]->GetYaxis()->SetTitle("Entries");
h_[ii][14] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_0123");
h_[ii][14]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][14]->GetYaxis()->SetTitle("Entries");
h_[ii][15] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_0123");
h_[ii][15]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][15]->GetYaxis()->SetTitle("Entries");
h_[ii][16] =(TH1F*) file[ii] ->Get("h_mz_0123");
h_[ii][16]->GetXaxis()->SetTitle("Z Mass of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][16]->GetYaxis()->SetTitle("Entries");
h_[ii][17] =(TH1F*) file[ii] ->Get("h_mz_1j_all_0123");
h_[ii][17]->GetXaxis()->SetTitle("Z mass for excl. Z+1 jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][17]->GetYaxis()->SetTitle("Entries");
h_[ii][18] =(TH1F*) file[ii] ->Get("h_mz_1j_c_0123");
h_[ii][18]->GetXaxis()->SetTitle("Z mass for excl. Z+1 central jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][18]->GetYaxis()->SetTitle("Entries");
h_[ii][19] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_0123");
h_[ii][19]->GetXaxis()->SetTitle("Z mass for excl. Z+1 HF jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][19]->GetYaxis()->SetTitle("Entries");
h_[ii][20] =(TH1F*) file[ii] ->Get("h_dielec_PT_0123");
h_[ii][20]->GetXaxis()->SetTitle("Z P_{T}  of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][20]->GetYaxis()->SetTitle("Entries");
h_[ii][21] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_0123");
h_[ii][21]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][21]->GetYaxis()->SetTitle("Entries");
h_[ii][22] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_0123");
h_[ii][22]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][22]->GetYaxis()->SetTitle("Entries");
h_[ii][23] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_0123");
h_[ii][23]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][23]->GetYaxis()->SetTitle("Entries");
h_[ii][24] =(TH1F*) file[ii] ->Get("h_dielecphi_0123");
h_[ii][24]->GetXaxis()->SetTitle("Z #phi of MVA(e1)<0.5 MVA(e2)<0.5 ");
h_[ii][24]->GetYaxis()->SetTitle("Entries");
h_[ii][25] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_0123");
h_[ii][25]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][25]->GetYaxis()->SetTitle("Entries");
h_[ii][26] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_0123");
h_[ii][26]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][26]->GetYaxis()->SetTitle("Entries");
h_[ii][27] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_0123");
h_[ii][27]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][27]->GetYaxis()->SetTitle("Entries");
h_[ii][28] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_0123");
h_[ii][28]->GetXaxis()->SetTitle("Z Y of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][28]->GetYaxis()->SetTitle("Entries");
h_[ii][29] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_0123");
h_[ii][29]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][29]->GetYaxis()->SetTitle("Entries");
h_[ii][30] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_0123");
h_[ii][30]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][30]->GetYaxis()->SetTitle("Entries");
h_[ii][31] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_0123");
h_[ii][31]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][31]->GetYaxis()->SetTitle("Entries");
h_[ii][32] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_0123");
h_[ii][32]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][32]->GetYaxis()->SetTitle("Entries");
h_[ii][33] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_hf_0123");
h_[ii][33]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][33]->GetYaxis()->SetTitle("Entries");
h_[ii][34] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_0123");
h_[ii][34]->GetXaxis()->SetTitle("# of good vertices of Z events of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][34]->GetYaxis()->SetTitle("Entries");
h_[ii][35] =(TH1F*) file[ii] ->Get("h_mz_same_sign_0123");
h_[ii][35]->GetXaxis()->SetTitle("same sign electron invariant mass of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][35]->GetYaxis()->SetTitle("Entries");
h_[ii][36] =(TH1F*) file[ii] ->Get("h_jet_pt_all_0123");
h_[ii][36]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][36]->GetYaxis()->SetTitle("Entries");
h_[ii][37] =(TH1F*) file[ii] ->Get("h_jet_eta_all0123");
h_[ii][37]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][37]->GetYaxis()->SetTitle("Entries");
h_[ii][38] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_0123");
h_[ii][38]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][38]->GetYaxis()->SetTitle("Entries");
h_[ii][39] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe0123");
h_[ii][39]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][39]->GetYaxis()->SetTitle("Entries");
h_[ii][40] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_0123");
h_[ii][40]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][40]->GetYaxis()->SetTitle("Entries");
h_[ii][41] =(TH1F*) file[ii] ->Get("h_jet_eta_hf0123");
h_[ii][41]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  of MVA(e1)<0.5 MVA(e2)<0.5");
h_[ii][41]->GetYaxis()->SetTitle("Entries");

h_[ii][42] =(TH1F*) file[ii] ->Get("h_elec_pt_0123456");
h_[ii][42]->GetXaxis()->SetTitle("Leading Electron P_{T} of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][42]->GetYaxis()->SetTitle("Entries");
h_[ii][43] =(TH1F*) file[ii] ->Get("h_elec_eta_0123456");
h_[ii][43]->GetXaxis()->SetTitle("Leading Electron #eta of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][43]->GetYaxis()->SetTitle("Entries");
h_[ii][44] =(TH1F*) file[ii] ->Get("h_numberofPFjets_0123456");
h_[ii][44]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][44]->GetYaxis()->SetTitle("Entries");
h_[ii][45] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_0123456");
h_[ii][45]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][45]->GetYaxis()->SetTitle("Entries");
h_[ii][46] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_0123456");
h_[ii][46]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][46]->GetYaxis()->SetTitle("Entries");
h_[ii][47] =(TH1F*) file[ii] ->Get("h_mz_0123456");
h_[ii][47]->GetXaxis()->SetTitle("Z Mass of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][47]->GetYaxis()->SetTitle("Entries");
h_[ii][48] =(TH1F*) file[ii] ->Get("h_mz_1j_all_0123456");
h_[ii][48]->GetXaxis()->SetTitle("Z mass for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][48]->GetYaxis()->SetTitle("Entries");
h_[ii][49] =(TH1F*) file[ii] ->Get("h_mz_1j_c_0123456");
h_[ii][49]->GetXaxis()->SetTitle("Z mass for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][49]->GetYaxis()->SetTitle("Entries");
h_[ii][50] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_0123456");
h_[ii][50]->GetXaxis()->SetTitle("Z mass for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][50]->GetYaxis()->SetTitle("Entries");
h_[ii][51] =(TH1F*) file[ii] ->Get("h_dielec_PT_0123456");
h_[ii][51]->GetXaxis()->SetTitle("Z P_{T}  of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][51]->GetYaxis()->SetTitle("Entries");
h_[ii][52] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_0123456");
h_[ii][52]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][52]->GetYaxis()->SetTitle("Entries");
h_[ii][53] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_0123456");
h_[ii][53]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][53]->GetYaxis()->SetTitle("Entries");
h_[ii][54] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_0123456");
h_[ii][54]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][54]->GetYaxis()->SetTitle("Entries");
h_[ii][55] =(TH1F*) file[ii] ->Get("h_dielecphi_0123456");
h_[ii][55]->GetXaxis()->SetTitle("Z #phi of MVA(e1)>0.5 MVA(e2)>0.5 ");
h_[ii][55]->GetYaxis()->SetTitle("Entries");
h_[ii][56] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_0123456");
h_[ii][56]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][56]->GetYaxis()->SetTitle("Entries");
h_[ii][57] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_0123456");
h_[ii][57]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][57]->GetYaxis()->SetTitle("Entries");
h_[ii][58] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_0123456");
h_[ii][58]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][58]->GetYaxis()->SetTitle("Entries");
h_[ii][59] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_0123456");
h_[ii][59]->GetXaxis()->SetTitle("Z Y of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][59]->GetYaxis()->SetTitle("Entries");
h_[ii][60] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_0123456");
h_[ii][60]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][60]->GetYaxis()->SetTitle("Entries");
h_[ii][61] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_0123456");
h_[ii][61]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][61]->GetYaxis()->SetTitle("Entries");
h_[ii][62] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_0123456");
h_[ii][62]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][62]->GetYaxis()->SetTitle("Entries");
h_[ii][63] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_0123456");
h_[ii][63]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][63]->GetYaxis()->SetTitle("Entries");
h_[ii][64] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_hf_0123456");
h_[ii][64]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][64]->GetYaxis()->SetTitle("Entries");
h_[ii][65] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_0123456");
h_[ii][65]->GetXaxis()->SetTitle("# of good vertices of Z events of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][65]->GetYaxis()->SetTitle("Entries");
h_[ii][66] =(TH1F*) file[ii] ->Get("h_mz_same_sign_0123456");
h_[ii][66]->GetXaxis()->SetTitle("same sign electron invariant mass of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][66]->GetYaxis()->SetTitle("Entries");
h_[ii][67] =(TH1F*) file[ii] ->Get("h_jet_pt_all_0123456");
h_[ii][67]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][67]->GetYaxis()->SetTitle("Entries");
h_[ii][68] =(TH1F*) file[ii] ->Get("h_jet_eta_all0123456");
h_[ii][68]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][68]->GetYaxis()->SetTitle("Entries");
h_[ii][69] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_0123456");
h_[ii][69]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][69]->GetYaxis()->SetTitle("Entries");
h_[ii][70] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe0123456");
h_[ii][70]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][70]->GetYaxis()->SetTitle("Entries");
h_[ii][71] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_0123456");
h_[ii][71]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][71]->GetYaxis()->SetTitle("Entries");
h_[ii][72] =(TH1F*) file[ii] ->Get("h_jet_eta_hf0123456");
h_[ii][72]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  of MVA(e1)>0.5 MVA(e2)>0.5");
h_[ii][72]->GetYaxis()->SetTitle("Entries");

h_[ii][73] =(TH1F*) file[ii] ->Get("h_elec_pt_01234567");
h_[ii][73]->GetXaxis()->SetTitle("Leading Electron P_{T} of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][73]->GetYaxis()->SetTitle("Entries");
h_[ii][74] =(TH1F*) file[ii] ->Get("h_elec_eta_01234567");
h_[ii][74]->GetXaxis()->SetTitle("Leading Electron #eta of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][74]->GetYaxis()->SetTitle("Entries");
h_[ii][75] =(TH1F*) file[ii] ->Get("h_numberofPFjets_01234567");
h_[ii][75]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][75]->GetYaxis()->SetTitle("Entries");
h_[ii][76] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hbhe_01234567");
h_[ii][76]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HBHE of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][76]->GetYaxis()->SetTitle("Entries");
h_[ii][77] =(TH1F*) file[ii] ->Get("h_numberofPFjets_hf_01234567");
h_[ii][77]->GetXaxis()->SetTitle("Exclusive Jet Multiplicity in HF of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][77]->GetYaxis()->SetTitle("Entries");
h_[ii][78] =(TH1F*) file[ii] ->Get("h_mz_01234567");
h_[ii][78]->GetXaxis()->SetTitle("Z Mass of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][78]->GetYaxis()->SetTitle("Entries");
h_[ii][79] =(TH1F*) file[ii] ->Get("h_mz_1j_all_01234567");
h_[ii][79]->GetXaxis()->SetTitle("Z mass for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][79]->GetYaxis()->SetTitle("Entries");
h_[ii][80] =(TH1F*) file[ii] ->Get("h_mz_1j_c_01234567");
h_[ii][80]->GetXaxis()->SetTitle("Z mass for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][80]->GetYaxis()->SetTitle("Entries");
h_[ii][81] =(TH1F*) file[ii] ->Get("h_mz_1j_hf_01234567");
h_[ii][81]->GetXaxis()->SetTitle("Z mass for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][81]->GetYaxis()->SetTitle("Entries");
h_[ii][82] =(TH1F*) file[ii] ->Get("h_dielec_PT_01234567");
h_[ii][82]->GetXaxis()->SetTitle("Z P_{T}  of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][82]->GetYaxis()->SetTitle("Entries");
h_[ii][83] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_all_01234567");
h_[ii][83]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][83]->GetYaxis()->SetTitle("Entries");
h_[ii][84] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_c_01234567");
h_[ii][84]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][84]->GetYaxis()->SetTitle("Entries");
h_[ii][85] =(TH1F*) file[ii] ->Get("h_dielec_PT_1j_hf_01234567");
h_[ii][85]->GetXaxis()->SetTitle("Z P_{T} for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][85]->GetYaxis()->SetTitle("Entries");
h_[ii][86] =(TH1F*) file[ii] ->Get("h_dielecphi_01234567");
h_[ii][86]->GetXaxis()->SetTitle("Z #phi of MVA(e1)>0.5 MVA(e2)<0.5 ");
h_[ii][86]->GetYaxis()->SetTitle("Entries");
h_[ii][87] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_all_01234567");
h_[ii][87]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][87]->GetYaxis()->SetTitle("Entries");
h_[ii][88] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_c_01234567");
h_[ii][88]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][88]->GetYaxis()->SetTitle("Entries");
h_[ii][89] =(TH1F*) file[ii] ->Get("h_dielecphi_1j_hf_01234567");
h_[ii][89]->GetXaxis()->SetTitle("Z #phi for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][89]->GetYaxis()->SetTitle("Entries");
h_[ii][90] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_01234567");
h_[ii][90]->GetXaxis()->SetTitle("Z Y of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][90]->GetYaxis()->SetTitle("Entries");
h_[ii][91] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_all_01234567");
h_[ii][91]->GetXaxis()->SetTitle("Z Y for excl. Z+1 jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][91]->GetYaxis()->SetTitle("Entries");
h_[ii][92] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_c_01234567");
h_[ii][92]->GetXaxis()->SetTitle("Z Y for excl. Z+1 central jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][92]->GetYaxis()->SetTitle("Entries");
h_[ii][93] =(TH1F*) file[ii] ->Get("h_dielec_rapidity_1j_hf_01234567");
h_[ii][93]->GetXaxis()->SetTitle("Z Y for excl. Z+1 HF jet events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][93]->GetYaxis()->SetTitle("Entries");
h_[ii][94] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_01234567");
h_[ii][94]->GetXaxis()->SetTitle("#DeltaR(e,j) for central electrons #& central jets of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][94]->GetYaxis()->SetTitle("Entries");
h_[ii][95] =(TH1F*) file[ii] ->Get("h_deltar_elec_PFjet_hf_01234567");
h_[ii][95]->GetXaxis()->SetTitle("#DeltaR(e,j)for HF electrons #& HF jets of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][95]->GetYaxis()->SetTitle("Entries");
h_[ii][96] =(TH1F*) file[ii] ->Get("h_no_good_vtx_zevents_01234567");
h_[ii][96]->GetXaxis()->SetTitle("# of good vertices of Z events of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][96]->GetYaxis()->SetTitle("Entries");
h_[ii][97] =(TH1F*) file[ii] ->Get("h_mz_same_sign_01234567");
h_[ii][97]->GetXaxis()->SetTitle("same sign electron invariant mass of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][97]->GetYaxis()->SetTitle("Entries");
h_[ii][98] =(TH1F*) file[ii] ->Get("h_jet_pt_all_01234567");
h_[ii][98]->GetXaxis()->SetTitle("Jet pt for leading jet in all detector of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][98]->GetYaxis()->SetTitle("Entries");
h_[ii][99] =(TH1F*) file[ii] ->Get("h_jet_eta_all01234567");
h_[ii][99]->GetXaxis()->SetTitle("Jet #eta for leading jet in all detector of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][99]->GetYaxis()->SetTitle("Entries");
h_[ii][100] =(TH1F*) file[ii] ->Get("h_jet_pt_hbhe_01234567");
h_[ii][100]->GetXaxis()->SetTitle("Jet pt for leading jet in HBHE  of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][100]->GetYaxis()->SetTitle("Entries");
h_[ii][101] =(TH1F*) file[ii] ->Get("h_jet_eta_hbhe01234567");
h_[ii][101]->GetXaxis()->SetTitle("Jet #eta for leading jet in HBHE  of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][101]->GetYaxis()->SetTitle("Entries");
h_[ii][102] =(TH1F*) file[ii] ->Get("h_jet_pt_hf_01234567");
h_[ii][102]->GetXaxis()->SetTitle("Jet pt for leading jet in HF  of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][102]->GetYaxis()->SetTitle("Entries");
h_[ii][103] =(TH1F*) file[ii] ->Get("h_jet_eta_hf01234567");
h_[ii][103]->GetXaxis()->SetTitle("Jet #eta for leading jet in HF  of MVA(e1)>0.5 MVA(e2)<0.5");
h_[ii][103]->GetYaxis()->SetTitle("Entries");

h_[ii][104] =(TH1F*) file[ii] ->Get("h_elec_mva_l_cut");
h_[ii][104]->GetXaxis()->SetTitle("Leading Electron MVA");
h_[ii][104]->GetYaxis()->SetTitle("Entries");
h_[ii][105] =(TH1F*) file[ii] ->Get("h_elec_mva_nl_cut");
h_[ii][105]->GetXaxis()->SetTitle("Second Leading Electron MVA");
h_[ii][105]->GetYaxis()->SetTitle("Entries");

h_[ii][106] =(TH1F*) file[ii] ->Get("h_elec_mva");
h_[ii][106]->GetXaxis()->SetTitle("Electron MVA (normalized to area)");
h_[ii][106]->GetYaxis()->SetTitle("Entries");

h_[ii][107] =(TH1F*) file[ii] ->Get("h_elec_mva_pt25");
h_[ii][107]->GetXaxis()->SetTitle("Electron MVA for P_{T}(e)>25 GeV. (normalized to area)");
h_[ii][107]->GetYaxis()->SetTitle("Entries");




for(int jj=0;jj<number_of_histos;jj++)h_[ii][jj]->GetXaxis()->SetTitleSize(0.028);

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

bool integral_norm = false;
bool same_sign_0 = false;
bool same_sign_1 = false;
bool same_sign_4 = false;
bool same_sign_5 = false;
if((ii>1 && ii<11)|| ii==106 || ii==107 )integral_norm = true;
if(ii==35) same_sign_0=true; if(ii==66) same_sign_1=true; if(ii==97) same_sign_4=true; 

cout<<ii<<same_sign_0<<endl;
c[ii] = new TCanvas();

h_[0][ii]->Sumw2();
if(!integral_norm){
//h_[1][ii]->Sumw2();
h_[1][ii]->Scale((3503.75*5097.)/27759468.);
//h_[2][ii]->Sumw2();
h_[2][ii]->Scale((3503.75*5097.)/27759468.);
//h_[3][ii]->Sumw2();
h_[3][ii]->Scale((225.197*5097.)/6642654.);
//h_[4][ii]->Sumw2();
h_[4][ii]->Scale((36257.2*5097.)/16487626.);
//h_[5][ii]->Sumw2();
h_[5][ii]->Scale((54.838*5097.)/10000431.);
//h_[6][ii]->Sumw2();
h_[6][ii]->Scale((32.31*5097.)/9996622.);
//h_[7][ii]->Sumw2();
h_[7][ii]->Scale((17.654*5097.)/9799908.);

h_[8][ii]->Scale((5745.25*5097.)/6665597.);
}
if(integral_norm){
h_[1][ii]->Scale(h_[0][ii]->Integral()/h_[1][ii]->Integral());
}
if(!integral_norm){
TH1F *diboson = (TH1F*)h_[5][ii]->Clone();
		diboson ->Add(h_[6][ii]);
		diboson ->Add(h_[7][ii]);

TH1F*mc_bg = (TH1F*)diboson->Clone();
		mc_bg->Add(h_[1][ii]);
		mc_bg->Add(h_[2][ii]);
		mc_bg->Add(h_[3][ii]);
		mc_bg->Add(h_[4][ii]);
		mc_bg->Add(h_[8][ii]);
}

if(integral_norm)TH1F*mc_bg = (TH1F*)h_[1][ii]->Clone();



sprintf(name1,"histo_z_plot/hist_%i.root",ii);
sprintf(name2,"histo_z_plot/hist_%i.pdf",ii);
c[ii]->cd();
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
h_[0][ii]->Draw("E1sames");
h_[0][ii]->SetTitle("CMS Data 5.09 fb^{-1}");
h_[0][ii]->SetMarkerStyle(2);
h_[0][ii]->SetLineWidth(1);

if(same_sign_0){
qcd_bckg_0 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_0->Add(h_[1][ii],-1);
}
if(same_sign_1){
qcd_bckg_1 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_1->Add(h_[1][ii],-1);
}
if(same_sign_4){
qcd_bckg_4 = (TH1F*)h_[0][ii]->Clone();
qcd_bckg_4->Add(h_[1][ii],-1);
}







//TCanvas->Return(c);
c[ii]->Print(name1);
c[ii]->Print(name2);

}

TCanvas *c_qcd_0 = new TCanvas();
TCanvas *c_qcd_1 = new TCanvas();
TCanvas *c_qcd_4 = new TCanvas();
c_qcd_0->cd();
qcd_bckg_0->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG) of MVA(e1)<0.5 MVA(e2)<0.5");
qcd_bckg_0->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_0->GetYaxis()->SetTitle("Entries");
qcd_bckg_0->Draw("E1");
c_qcd_0->Print("histo_z_plot/qcd_bckground_0.root");
c_qcd_0->Print("histo_z_plot/qcd_bckground_0.pdf");
c_qcd_1->cd();
qcd_bckg_1->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG) of MVA(e1)>0.5 MVA(e2)>0.5");
qcd_bckg_1->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_1->GetYaxis()->SetTitle("Entries");
qcd_bckg_1->Draw("E1");
c_qcd_1->Print("histo_z_plot/qcd_bckground_1.root");
c_qcd_1->Print("histo_z_plot/qcd_bckground_1.pdf");
c_qcd_4->cd();
qcd_bckg_4->GetXaxis()->SetTitle("same-sign ee inv.M. Data-(MC+BG) of MVA(e1)>0.5 MVA(e2)<0.5");
qcd_bckg_4->GetXaxis()->SetTitleSize(0.03);
qcd_bckg_4->GetYaxis()->SetTitle("Entries");
qcd_bckg_4->Draw("E1");
c_qcd_4->Print("histo_z_plot/qcd_bckground_4.root");
c_qcd_4->Print("histo_z_plot/qcd_bckground_4.pdf");


}
