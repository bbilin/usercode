{

TFile * file = TFile::Open("mc.root");

TH2F * h_n_qjet = file ->Get("h_ptz_etajet_q");

TH2F * h_n_gjet = file ->Get("h_ptz_etajet_g");

TH2F * h_n_total = file ->Get("h_ptz_etajet");

TProfile2D * h_ptz_avg_jet = file ->Get("h_ptZ_etajet__average_PTZ");

TProfile2D * h_ptjet_avg_jet = file ->Get("h_ptZ_etajet__average_ETraw");

TProfile2D * h_ptz_avg_qjet = file ->Get("h_ptZ_etajet__average_PTZ_q");

TProfile2D * h_ptz_avg_gjet = file ->Get("h_ptZ_etajet__average_PTZ_g");

TProfile2D * h_ptjet_avg_qjet = file ->Get("h_ptZ_etajet__average_ETraw_q");

TProfile2D * h_ptjet_avg_gjet = file ->Get("h_ptZ_etajet__average_ETraw_g");

/////////////////////////
//&&&&&&&&&&&&&&&&&&&&&//
/////////////////////////

TFile * data = TFile::Open("../calib_anaslyzer/data.root");

TH2F * h_n_total_data = data ->Get("h_ptz_etajet");

TProfile2D * h_ptz_avg_data = data ->Get("h_ptZ_etajet__average_PTZ");

TProfile2D * h_ptjet_avg_data = data ->Get("h_ptZ_etajet__average_ETraw");

/////////////////////////
//&&&&&&&&&&&&&&&&&&&&&//
/////////////////////////

double omega_q[50][50]={};
double omega_g[50][50]={};
double k_q_g[50][50]={};
double C_data[50][50]={};
double c_q[50][50]={};
double data_avg_etraw[50][50]={};
double error_data_avg_etraw[50][50]={};
double c_q_err[50][50] ={};
double nevt[50][50]={};
  char output_file[200];
 ofstream outf[10];

for(int j = 1 ; j<11; j++){
//    sprintf(output_file,"file%i.txt",j-1);
//    outf[j-1].open(output_file);
	for (int i=1; i<11;i++){

C_data[i][j]=h_ptz_avg_data->GetBinContent(i,j)/h_ptjet_avg_data->GetBinContent(i,j);
data_avg_etraw[i][j] = h_ptjet_avg_data->GetBinContent(i,j);
error_data_avg_etraw[i][j] = h_ptjet_avg_data->GetBinError(i,j);

omega_q[i][j] =h_n_qjet->GetBinContent(i,j)/h_n_total->GetBinContent(i,j);
omega_g[i][j] =h_n_gjet->GetBinContent(i,j)/h_n_total->GetBinContent(i,j);

k_q_g[i][j]=(h_ptz_avg_qjet->GetBinContent(i,j)/h_ptjet_avg_qjet->GetBinContent(i,j))/(h_ptz_avg_gjet->GetBinContent(i,j)/h_ptjet_avg_gjet->GetBinContent(i,j));
nevt[i][j]=h_n_total_data->GetBinContent(i,j);
//if(h_ptjet_avg_data->GetBinContent(i,j)==0||h_n_total->GetBinContent(i,j)==0||h_ptjet_avg_gjet->GetBinContent(i,j)==0||h_ptjet_avg_qjet->GetBinContent(i,j)==0||h_ptz_avg_gjet->GetBinContent(i,j)==0||nevt[i][j]<10||h_n_total->GetBinContent(i,j)<10) continue;
if(/*nevt[i][j]<100||*/h_n_total->GetBinContent(i,j)<100) continue;
//c_q[i][j]= C_data[i][j] *(omega_q[i][j]  + k_q_g[i][j]*omega_g[i][j]);
c_q[i][j]=(h_ptz_avg_qjet->GetBinContent(i,j)/h_ptjet_avg_qjet->GetBinContent(i,j)) *(omega_q[i][j]  + k_q_g[i][j]*omega_g[i][j]);

c_q_err[i][j]=c_q[i][j]/sqrt(nevt[i][j]);



//index j corr to eta

//cout<<
}
}

/////////////////////////
//&&&&&&&&&&&&&&&&&&&&&//
/////////////////////////

float c_draw_yq[10][10]={};
float c_draw_x[10][10]={};
float err_x[10][10]={};
float err_y[10][10]={};
cout<<"double ";
for(int i=0;i<10;i++){
	for(int j=0;j<10;j++){
c_draw_yq[i][j]=c_q[j+1][i+1];
cout<<"c["<<i<<"]"<<"["<<j<<"]"<<"="<<c_draw_yq[i][j]<<",";
c_draw_x[i][j]=data_avg_etraw[j+1][i+1];
err_x[i][j]=error_data_avg_etraw[j+1][i+1];
err_y[i][j]=c_q_err[j+1][i+1];
//cout<<c_draw_x[i][j]<<endl;
//cout<<i<<"  "<<j<<""<<c_draw_yq[i][j]<<endl;
}
}
cout<<endl;
float title[]={0,0.519,1.038,1.557,2.076,2.595,3.114,3.633,4.152,4.671,5.191};
TCanvas*c[50];
TGraphErrors*gr[50];
for(int i=0;i<10;i++){


    TString name1 = "|#eta|=";
    name1 += (i)*0.519f;
    name1 += "-";
    name1 += (i+1)*0.519f;

     TString name2 = "quark_ratio_plots/eta";
    name2 += (i)*0.5;
    name2 += "-";
    name2 += (i+1.)*0.5;
    name2 += ".root";

     TString name3 = "quark_ratio_plots/eta";
    name3 += (i)*0.5;
    name3 += "-";
    name3 += (i+1)*0.5;
    name3 += ".pdf";
    


c[i]=new TCanvas;
c[i]->cd();
   gr[i] = new TGraphErrors(10,c_draw_x[i],c_draw_yq[i],err_x[i],err_y[i]);
   
    gr[i] ->GetXaxis()->SetTitle("<E_{T}^{raw}> (GeV)");
    gr[i] ->GetYaxis()->SetTitle("<C_{q}>");
    gr[i]->SetTitle(name1);
   gr[i]->Draw("A*");
 
c[i]->SaveAs(name2);
c[i]->SaveAs(name3);
}

/*  TProfile *h_ratio[50];
for(j = 1;j<11;j++){
    TString name1 = "quark raw. |#eta|=";
    name1 += (j-1)*0.591;
    name1 += "-";
    name1 += (j)*0.591;
 h_ratio[j] = new TProfile(name1,name1,8,20,180);

//	for(i = 1;i<11;j++){
//	 h_ratio[j]->Fill(data_avg_etraw[i][j],c_q[i][j]);
//}
//cout<<"burda mı"<<endl;
}
for(i = 1;i<11;j++){
h_ratio[1]->Fill(data_avg_etraw[i][1],c_q[i][1],1);
}

/*TCanvas *c[100];
for(j = 1;j<11;j++){

     TString name2 = "|#eta|=";
    name2 += (j-1)*0.591;
    name2 += "-";
    name2 += (j)*0.591;
    name2 += ".root";

c[j]=new TCanvas;
c[j]->cd();
    h_ratio[j]->GetXaxis()->SetTitle("<E_{T}^{raw}> (GeV)");
    h_ratio[j]->GetYaxis()->SetTitle("<C_{q}>");
    h_ratio[j]->Draw();
c[j]->SaveAs(name2);
}

TCanvas *c= new TCanvas;
c->cd();
    h_ratio[1]->GetXaxis()->SetTitle("<E_{T}^{raw}> (GeV)");
    h_ratio[1]->GetYaxis()->SetTitle("<C_{q}>");
    h_ratio[1]->Draw();
c[1]->SaveAs("test.root");

*/
}
