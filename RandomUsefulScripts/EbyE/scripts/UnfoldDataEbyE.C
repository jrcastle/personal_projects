#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "TLatex.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>

bool iter256 = false;
bool resp_data = false;

//int centbin = 0;  // 0 - 5 %
//int centbin = 1;  // 5 - 10 %
//int centbin = 2;  // 10 - 15 %
//int centbin = 3;  // 15 - 20 %
int centbin = 4;  // 20 - 25 %
//int centbin = 5;  // 25 - 30 %
//int centbin = 6;  // 30 - 35 %
//int centbin = 7;  // 35 - 40 %
//int centbin = 8;  // 40 - 45 %
//int centbin = 9;  // 45 - 50 %
//int centbin = 10; // 50 - 55 %
//int centbin = 11; // 55 - 60 %
//int centbin = 12; // 60 - 65 %
//int centbin = 13; // 65 - 70 %
//int centbin = 14; // 70 - 75 %
//int centbin = 15; // 75 - 80 %
//int centbin = 16; // 80 - 85 %
//int centbin = 17; // 85 - 90 %
//int centbin = 18; // 90 - 95 %
//int centbin = 19; // 95 - 100 %

int Nbins = 150;
double v2Min = 0.0;
double v2Max = 0.6;
int VN = 2;
int sw = 1;
static const int NCENT40 = 40;
static const int NCENT = 20;
//double v2_2SEdiffMax[NCENT40] = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6 }; //DATA
double v2_2SEdiffMax[NCENT40] = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6 }; //MC

//-- Unfolding Histos
TH2D * h2D[NCENT];
TH2D * h2Dsub0[NCENT];
TH2D * h2Dsub1[NCENT];
TH2D * h2Dx[NCENT];

TH1D * h1D[NCENT];
TH1D * h1Dsub0[NCENT];
TH1D * h1Dsub1[NCENT];
TH1D * h1Dx[NCENT];
TH1D * h1Dy[NCENT];

TH1D * hw[NCENT];

TF1 * f1Dx[NCENT];
TF1 * f1Dy[NCENT];

double sigma_x[NCENT];
double sigma_xe[NCENT];
double mean_x[NCENT];
double chi2_x[NCENT];
double ndf_x[NCENT];
double sigma_y[NCENT];
double sigma_ye[NCENT];
double mean_y[NCENT];
double chi2_y[NCENT];
double ndf_y[NCENT];
double sigma_2SE[NCENT];
double Vn_mean[NCENT];
double Vn_rms[NCENT];

TH2D * hresp[NCENT];
TH1D * hreco1[NCENT];
TH1D * hreco2[NCENT];
TH1D * hreco4[NCENT];
TH1D * hreco8[NCENT];
TH1D * hreco16[NCENT];
TH1D * hreco32[NCENT];
TH1D * hreco64[NCENT];
TH1D * hreco128[NCENT];
TH1D * hreco256[NCENT];

TH1D * h1Dr[NCENT];
TH1D * hreco1r[NCENT];
TH1D * hreco2r[NCENT];
TH1D * hreco4r[NCENT];
TH1D * hreco8r[NCENT];
TH1D * hreco16r[NCENT];
TH1D * hreco32r[NCENT];
TH1D * hreco64r[NCENT];

//-- Data Histos
static const int NVn = 7;
TH2D * hVn2Dfull[NVn][NCENT40];
TH2D * hVn2Dsub0[NVn][NCENT40];
TH2D * hVn2Dsub1[NVn][NCENT40];
TH2D * hVn2D0v1[NVn][NCENT40];
TH1D * hVnFull[NVn][NCENT40];
TH1D * hVnSub0[NVn][NCENT40];
TH1D * hVnSub1[NVn][NCENT40];

void UnfoldDataEbyE(){
  
  TFile * fsave = new TFile(Form("~/root/macros/EbyE/macros/txt/PbPb_2015/data%i.root", VN), "recreate");

  //--Get histos from datafile
  TFile * fData = new TFile("/Users/jcastle/root/macros/EbyE/macros/data/PbPb_2015/CastleEbyE.root");
  for(int ivn = 1; ivn <=4; ivn++){
    if(ivn != VN) continue;

    for(int icent = 0; icent < NCENT40; icent++){

      hVn2Dfull[ivn][icent]  = (TH2D*) fData->Get(Form("qwebye/hVn2Dfull_%i_%i",ivn,icent));
      hVn2Dsub0[ivn][icent]  = (TH2D*) fData->Get(Form("qwebye/hVn2Dsub0_%i_%i",ivn,icent));
      hVn2Dsub1[ivn][icent]  = (TH2D*) fData->Get(Form("qwebye/hVn2Dsub1_%i_%i",ivn,icent));
      hVn2D0v1[ivn][icent]   = (TH2D*) fData->Get(Form("qwebye/hVn2D0v1_%i_%i",ivn,icent));
      hVnFull[ivn][icent]    = (TH1D*) fData->Get(Form("qwebye/hVnFull_%i_%i",ivn,icent));
      hVnSub0[ivn][icent]    = (TH1D*) fData->Get(Form("qwebye/hVnSub0_%i_%i",ivn,icent));
      hVnSub1[ivn][icent]    = (TH1D*) fData->Get(Form("qwebye/hVnSub1_%i_%i",ivn,icent));

    }


  }

  //-- Load unfolding library and declare unfolding objects
  //gSystem->Load("~/root/RooUnfold-1.1.1/libRooUnfold");
  RooUnfoldResponse * response[NCENT];
  RooUnfoldBayes * unfold1[NCENT];
  RooUnfoldBayes * unfold2[NCENT];
  RooUnfoldBayes * unfold4[NCENT];
  RooUnfoldBayes * unfold8[NCENT];
  RooUnfoldBayes * unfold16[NCENT];
  RooUnfoldBayes * unfold32[NCENT];
  RooUnfoldBayes * unfold64[NCENT];
  RooUnfoldBayes * unfold128[NCENT];
  RooUnfoldBayes * unfold256[NCENT];

  TLatex latex;
  latex.SetNDC();

  //-- Begin Centrality loop
  for(int c = 0; c < NCENT; c++){
    
    if(c != centbin) continue;
    std::cout<<"!! Processing Cent = "<<c<<std::endl;


    //-- Initiate unfolding histos and condense from 40 to 20 cent bins
    h2D[c] = (TH2D*)hVn2Dfull[VN][2*c]->Clone(Form("h2D_%i_%i", VN, c));
    h2D[c]->Add(hVn2Dfull[VN][2*c+1]);
    if ( h2D[c]->GetEntries() < 1000 ) continue;
    h2Dsub0[c] = (TH2D*)hVn2Dsub0[VN][2*c]->Clone(Form("h2Dsub0_%i_%i", VN, c));
    h2Dsub0[c]->Add(hVn2Dsub0[VN][2*c+1]);
    h2Dsub1[c] = (TH2D*)hVn2Dsub1[VN][2*c]->Clone(Form("h2Dsub1_%i_%i", VN, c));
    h2Dsub1[c]->Add(hVn2Dsub1[VN][2*c+1]);
    h2Dx[c] = (TH2D*)hVn2D0v1[VN][2*c]->Clone(Form("h2Dx_%i_%i", VN, c));
    h2Dx[c]->Add(hVn2D0v1[VN][2*c+1]);

    h1D[c] = (TH1D*)hVnFull[VN][2*c]->Clone(Form("h1D_%i_%i", VN, c));
    h1D[c]->Add(hVnFull[VN][2*c+1]);
    h1Dsub0[c] = (TH1D*)hVnSub0[VN][2*c]->Clone(Form("h1Dsub0_%i_%i", VN, c));
    h1Dsub0[c]->Add(hVnSub0[VN][2*c+1]);
    h1Dsub1[c] = (TH1D*)hVnSub1[VN][2*c]->Clone(Form("h1Dsub1_%i_%i", VN, c));
    h1Dsub1[c]->Add(hVnSub1[VN][2*c+1]);

    h1Dx[c] = h2Dx[c]->ProjectionX(Form("h1Dx_%i_%i", VN, c));
    h1Dy[c] = h2Dx[c]->ProjectionY(Form("h1Dy_%i_%i", VN, c));

    f1Dx[c] = new TF1(Form("f1Dx_%i_%i", VN, c), "gaus", -v2_2SEdiffMax[2*c], v2_2SEdiffMax[2*c]);
    f1Dy[c] = new TF1(Form("f1Dy_%i_%i", VN, c), "gaus", -v2_2SEdiffMax[2*c], v2_2SEdiffMax[2*c]);

    h1Dx[c]->Fit(f1Dx[c], "NLM");
    h1Dy[c]->Fit(f1Dy[c], "NLM");

    sigma_x[c] = f1Dx[c]->GetParameter("Sigma");
    sigma_xe[c] = f1Dx[c]->GetParError(f1Dx[c]->GetParNumber("Sigma"));
    mean_x[c]  = f1Dx[c]->GetParameter("Mean");
    chi2_x[c]  = f1Dx[c]->GetChisquare();

    sigma_y[c] = f1Dy[c]->GetParameter("Sigma");
    sigma_ye[c] = f1Dy[c]->GetParError(f1Dy[c]->GetParNumber("Sigma"));
    mean_y[c]  = f1Dy[c]->GetParameter("Mean");
    chi2_y[c]  = f1Dy[c]->GetChisquare();

    sigma_2SE[c] = 0.5*(sigma_x[c]+sigma_y[c]);
    ndf_x[c] = f1Dx[c]->GetNDF();
    ndf_y[c] = f1Dy[c]->GetNDF();

    double sigma = sigma_2SE[c]/2.;

    //Set up response function
    hresp[c] = new TH2D(Form("hresp_%i_%i", VN, c), "hresp", Nbins, 0., v2Max, Nbins, 0., v2Max);

    if ( sw ) {
      hw[c] = (TH1D*) h1D[c]->Clone(Form("hw_%i", c));
    }
    for ( int i = 1; i <= Nbins; i++ ) {
      for ( int j = 1; j <= Nbins; j++ ) {
	double w = 1.;
	if ( sw ) {
	  w = hw[c]->GetBinContent(j);
	}
	double v_mess = hresp[c]->GetXaxis()->GetBinCenter(i);
	double v_true = hresp[c]->GetYaxis()->GetBinCenter(j);
	double resp = v_mess * TMath::Gaus(sqrt(v_mess*v_mess + v_true*v_true), 0, sigma) * TMath::BesselI0( v_mess*v_true/sigma/sigma );
	//if ( i == 1 ) cout << "!!! i = " << i << "\t j = " << j << "\t resp = " << resp << endl;                                                                                                                     
	if ( TMath::IsNaN(resp) || TMath::Infinity()==resp ) resp = 0;
	hresp[c]->SetBinContent(i, j, resp*w);

      }
    }

    response[c] = new RooUnfoldResponse(0, 0, hresp[c], Form("response_%i_%i", VN, c));

    //-- Unfold!
    unfold1[c] = new RooUnfoldBayes( response[c], h1D[c], 1 );
    unfold2[c] = new RooUnfoldBayes( response[c], h1D[c], 2 );
    unfold4[c] = new RooUnfoldBayes( response[c], h1D[c], 4 );
    unfold8[c] = new RooUnfoldBayes( response[c], h1D[c], 8 );
    unfold16[c] = new RooUnfoldBayes( response[c], h1D[c], 16 );
    unfold32[c] = new RooUnfoldBayes( response[c], h1D[c], 32 );
    unfold64[c] = new RooUnfoldBayes( response[c], h1D[c], 64 );
    unfold128[c] = new RooUnfoldBayes( response[c], h1D[c], 128 );
    //unfold128[c] = new RooUnfoldBayes( response[c], h1D[c], 2048 );
    if(iter256) unfold256[c] = new RooUnfoldBayes( response[c], h1D[c], 256 );
    else unfold256[c] = 0;

    hreco1[c] = (TH1D*) unfold1[c]->Hreco();
    hreco1[c]->SetName(Form("hreco1_%i_%i", VN, c));

    hreco2[c] = (TH1D*) unfold2[c]->Hreco();
    hreco2[c]->SetName(Form("hreco2_%i_%i", VN, c));

    hreco4[c] = (TH1D*) unfold4[c]->Hreco();
    hreco4[c]->SetName(Form("hreco4_%i_%i", VN, c));

    hreco8[c] = (TH1D*) unfold8[c]->Hreco();
    hreco8[c]->SetName(Form("hreco8_%i_%i", VN, c));

    hreco16[c] = (TH1D*) unfold16[c]->Hreco();
    hreco16[c]->SetName(Form("hreco16_%i_%i", VN, c));

    hreco32[c] = (TH1D*) unfold32[c]->Hreco();
    hreco32[c]->SetName(Form("hreco32_%i_%i", VN, c));

    hreco64[c] = (TH1D*) unfold64[c]->Hreco();
    hreco64[c]->SetName(Form("hreco64_%i_%i", VN, c));

    hreco128[c] = (TH1D*) unfold128[c]->Hreco();
    hreco128[c]->SetName(Form("hreco128_%i_%i", VN, c));

    if(iter256){

      hreco256[c] = (TH1D*) unfold256[c]->Hreco();
      hreco256[c]->SetName(Form("hreco256_%i_%i", VN, c));
      hreco256[c]->GetXaxis()->SetTitle("v_{2}^{unfold}");

    }  

    //Save Files
    fsave->cd();
    h2D[c]->Write();
    h2Dsub0[c]->Write();
    h2Dsub1[c]->Write();
    h2Dx[c]->Write();

    h1D[c]->Write();
    h1Dsub0[c]->Write();
    h1Dsub1[c]->Write();
    h1Dx[c]->Write();
    h1Dy[c]->Write();

    f1Dx[c]->Write();
    f1Dy[c]->Write();

    hresp[c]->Write();
    hreco1[c]->Write();
    hreco2[c]->Write();
    hreco4[c]->Write();
    hreco8[c]->Write();
    hreco16[c]->Write();
    hreco32[c]->Write();
    hreco64[c]->Write();
    hreco128[c]->Write();
    if(iter256) hreco256[c]->Write();

  }




}
