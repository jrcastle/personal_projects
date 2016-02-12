#include "/home/j550c590/tdrstyle.C"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>

//unsigned int selectedRun = 182972;
double vtxCut = 15.;
TString fileName = "";

const int NCENT_PBPB = 40;
double cent_min[NCENT_PBPB] = {0.0 ,2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5};
double cent_max[NCENT_PBPB] = {2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0};
double v2_2SEdiffMax[NCENT_PBPB] = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6 };

static const int Nbins = 150;
static const int Nbins_2SEDiff = 300;
static const double v2Min = 0.0;
static const double v2Max = 0.6;
static const int VN = 2;

static const int nptbinsDefault = 11;
static const double ptbinsDefault[]={1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 7.00, 8.00};
static const int netabinsDefault = 14;
static const double etabinsDefault[]= {-2.4, -2.0, -1.6, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4};

static const int ptBinMin = 4;  //-- for(int ipt = ptBinMin; ipt <= ptBinMax; ipt++) --> GetBinContent(ipt+1,ieta+1);
static const int ptBinMax = 4;
static const int etaBinMin = 4; //-- for(int ieta = etaBinMin; ieta <= etaBinMax; ieta++) --> GetBinContent(ipt+1,ieta+1); 
static const int etaBinMax = 9;


TTree * tree;
unsigned int runno_;
double centval;
double vtx;
TH2D * sumw;
TH2D * sumwqx;
TH2D * sumwqy;

TFile * tf;
TFile * fV2Det;

TH1D * hV2DetX_TkEff_0;
TH1D * hV2DetY_TkEff_0;
TH1D * hV2DetX_TkEff_1;
TH1D * hV2DetY_TkEff_1;
TH1D * hV2DetX_TkEff_full;
TH1D * hV2DetY_TkEff_full;

//
// MAIN
//
void ReadTree_normDet(){
    
    setTDRStyle();
    TH1D::SetDefaultSumw2();

    tf = new TFile(fileName);
    fV2Det = new TFile("V2Det.root");

    sumwqx = new TH2D("sumwqx","sumwqx",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    sumwqy = new TH2D("sumwqy","sumwqy",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    sumw =  new TH2D("sumw","sumw",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);

    tree = (TTree *) tf->Get("ebyeana/tree");
    tree->SetBranchAddress("Run",         &runno_);
    tree->SetBranchAddress("Cent",        &centval);
    tree->SetBranchAddress("Vtx",         &vtx);
    tree->SetBranchAddress("sumw",        &sumw);
    tree->SetBranchAddress("sumwqx",      &sumwqx);
    tree->SetBranchAddress("sumwqy",      &sumwqy);
    
    hV2DetX_TkEff_0 = (TH1D *) fV2Det->Get("SubEvt_0/hV2DetX_TkEff_0");
    hV2DetY_TkEff_0 = (TH1D *) fV2Det->Get("SubEvt_0/hV2DetY_TkEff_0");
    hV2DetX_TkEff_1 = (TH1D *) fV2Det->Get("SubEvt_1/hV2DetX_TkEff_1");
    hV2DetY_TkEff_1 = (TH1D *) fV2Det->Get("SubEvt_1/hV2DetY_TkEff_1");
    hV2DetX_TkEff_full = (TH1D *) fV2Det->Get("FullEvt/hV2DetX_TkEff_full");
    hV2DetY_TkEff_full = (TH1D *) fV2Det->Get("FullEvt/hV2DetY_TkEff_full");
    
    TFile * fHists = new TFile("CastleEbyE.root","recreate");
    
    TH2D * hVn2Dfull[7][NCENT_PBPB];
    TH2D * hVn2Dsub0[7][NCENT_PBPB];
    TH2D * hVn2Dsub1[7][NCENT_PBPB];
    TH2D * hVn2D0v1[7][NCENT_PBPB];
    
    TH1D * hVnFull[7][NCENT_PBPB];
    TH1D * hVnSub0[7][NCENT_PBPB];
    TH1D * hVnSub1[7][NCENT_PBPB];
    
    TDirectory * qwebye = fHists->mkdir("qwebye");
    
    for ( int c = 0; c < NCENT_PBPB; c++ ) {

        for ( int n = 1; n < 7; n++ ) {

	  if(n != VN) continue;
	  qwebye->cd();
	  hVn2Dfull[n][c] = new TH2D(Form("hVn2Dfull_%i_%i", n, c), Form("hVn2Dfull_%i_%i", n, c), Nbins, -v2Max, v2Max, Nbins, -v2Max, v2Max);
	  hVn2Dfull[n][c]->SetOption("colz");
	  hVn2Dfull[n][c]->GetXaxis()->SetTitle(Form("v_{%i,x}^{obs}",n));
	  hVn2Dfull[n][c]->GetYaxis()->SetTitle(Form("v_{%i,y}^{obs}",n));

	  hVn2Dsub0[n][c] = new TH2D(Form("hVn2Dsub0_%i_%i", n, c), Form("hVn2Dsub0_%i_%i", n, c), Nbins, -v2Max, v2Max, Nbins, -v2Max, v2Max);
	  hVn2Dsub0[n][c]->SetOption("colz");
	  hVn2Dsub0[n][c]->GetXaxis()->SetTitle(Form("v_{%i,x}^{obs,a}",n));
	  hVn2Dsub0[n][c]->GetYaxis()->SetTitle(Form("v_{%i,y}^{obs,a}",n));

	  hVn2Dsub1[n][c] = new TH2D(Form("hVn2Dsub1_%i_%i", n, c), Form("hVn2Dsub1_%i_%i", n, c), Nbins, -v2Max, v2Max, Nbins, -v2Max, v2Max);
	  hVn2Dsub1[n][c]->SetOption("colz");
	  hVn2Dsub1[n][c]->GetXaxis()->SetTitle(Form("v_{%i,x}^{obs,b}",n));
	  hVn2Dsub1[n][c]->GetYaxis()->SetTitle(Form("v_{%i,y}^{obs,b}",n));

	  hVn2D0v1[n][c] = new TH2D(Form("hVn2D0v1_%i_%i", n, c), Form("hVn2D0v1_%i_%i", n, c), Nbins_2SEDiff, -v2_2SEdiffMax[c], v2_2SEdiffMax[c], Nbins_2SEDiff, -v2_2SEdiffMax[c], v2_2SEdiffMax[c]);
	  hVn2D0v1[n][c]->SetOption("colz");
	  hVn2D0v1[n][c]->GetXaxis()->SetTitle(Form("v_{%i,x}^{obs,a} - v_{%i,x}^{obs,b}",n,n));
	  hVn2D0v1[n][c]->GetYaxis()->SetTitle(Form("v_{%i,y}^{obs,a} - v_{%i,y}^{obs,b}",n,n));

	  hVnFull[n][c] = new TH1D(Form("hVnFull_%i_%i", n, c), Form("hVnFull_%i_%i", n, c), Nbins, v2Min, v2Max);
	  hVnFull[n][c]->GetXaxis()->SetTitle(Form("v_{%i}^{obs}",n));

	  hVnSub0[n][c] = new TH1D(Form("hVnSub0_%i_%i", n, c), Form("hVnSub0_%i_%i", n, c), Nbins, v2Min, v2Max);
	  hVnSub0[n][c]->GetXaxis()->SetTitle(Form("v_{%i}^{obs,a}",n));

	  hVnSub1[n][c] = new TH1D(Form("hVnSub1_%i_%i", n, c), Form("hVnSub1_%i_%i", n, c), Nbins, v2Min, v2Max);
	  hVnSub1[n][c]->GetXaxis()->SetTitle(Form("v_{%i}^{obs,b}",n));
        }
        
    }
    
    double VnRaw_x_0;
    double VnRaw_y_0;
    double VnRaw_x_1;
    double VnRaw_y_1;
    double VnRaw_x_full;
    double VnRaw_y_full;

    double sumw_0;
    double sumw_1;
    double sumw_full;
    
    double VnCorrected_x_0;
    double VnCorrected_y_0;
    double VnCorrected_x_1;
    double VnCorrected_y_1;
    double VnCorrected_x_full;
    double VnCorrected_y_full;
    
    //
    // Tree Loop
    //
    
    cout<<"Begin EVENT loop, contains "<<tree->GetEntries()<<" Events"<<endl;
    
    //for(int ievent = 0; ievent<1000000; ievent++) {
    for(int ievent = 0; ievent<tree->GetEntries(); ievent++) {
        
      if((ievent+1)% 500000 == 0) cout<<"Processing Event "<<ievent+1<<"\t"<<100.*(ievent+1)/(double)tree->GetEntries()<<"% Completed"<<endl;
      
        tree->GetEntry(ievent);

        //if(runno_ != selectedRun) continue;
	if(TMath::Abs(vtx) > vtxCut) continue;
      
        //-- Calculate cent bin
	int icent = (centval - cent_min[0]) / 2.5;

	//-- Reset raw and sumw values
	VnRaw_x_0    = 0;
	VnRaw_y_0    = 0;
	VnRaw_x_1    = 0;
	VnRaw_y_1    = 0;
	VnRaw_x_full = 0;
	VnRaw_y_full = 0;                

	sumw_0       = 0;
	sumw_1       = 0;
	sumw_full    = 0;

	for(int ipt = ptBinMin; ipt <= ptBinMax; ipt++){
	  for(int ieta = etaBinMin; ieta <= etaBinMax; ieta++){

	    if(sumw->GetBinContent(ipt+1,ieta+1) !=0){

	      //-- Subevent 0 (eta >= 0)
	      if(etabinsDefault[ieta] >= 0){
		VnRaw_x_0     += sumwqx->GetBinContent(ipt+1,ieta+1);
		VnRaw_y_0     += sumwqy->GetBinContent(ipt+1,ieta+1);
		sumw_0        += sumw->GetBinContent(ipt+1,ieta+1);
	      }
	      //-- Subevent 1 (eta < 0)
	      else{
		VnRaw_x_1     += sumwqx->GetBinContent(ipt+1,ieta+1);
		VnRaw_y_1     += sumwqy->GetBinContent(ipt+1,ieta+1);
		sumw_1        += sumw->GetBinContent(ipt+1,ieta+1);
	      }
	      //-- Full Event
	      VnRaw_x_full    += sumwqx->GetBinContent(ipt+1,ieta+1);
	      VnRaw_y_full    += sumwqy->GetBinContent(ipt+1,ieta+1);
	      sumw_full       += sumw->GetBinContent(ipt+1,ieta+1);

	    }

	  }

	}

	//-- Fill Histograms
	if(sumw_full == 0 || sumw_0 == 0 || sumw_1 == 0) continue;

	VnRaw_x_full /= sumw_full;
	VnRaw_y_full /= sumw_full;

	VnRaw_x_0 /= sumw_0;
	VnRaw_y_0 /= sumw_0;

	VnRaw_x_1 /= sumw_1;
	VnRaw_y_1 /= sumw_1;
                
	//-- Full Tracker
	VnCorrected_x_full = VnRaw_x_full - hV2DetX_TkEff_full->GetBinContent(icent+1);
	VnCorrected_y_full = VnRaw_y_full - hV2DetY_TkEff_full->GetBinContent(icent+1);

	hVnFull[2][icent]->Fill(TMath::Sqrt( VnCorrected_x_full * VnCorrected_x_full +  VnCorrected_y_full * VnCorrected_y_full));
	hVn2Dfull[2][icent]->Fill(VnCorrected_x_full,VnCorrected_y_full);
                
	//-- SubEvt 0 (Eta > 0)
	VnCorrected_x_0 = VnRaw_x_0 - hV2DetX_TkEff_0->GetBinContent(icent+1);
	VnCorrected_y_0 = VnRaw_y_0 - hV2DetY_TkEff_0->GetBinContent(icent+1);
        
	hVnSub0[2][icent]->Fill( TMath::Sqrt( VnCorrected_x_0 * VnCorrected_x_0 + VnCorrected_y_0 * VnCorrected_y_0 ) );
	hVn2Dsub0[2][icent]->Fill(VnCorrected_x_0, VnCorrected_y_0);
                
	//-- SubEvt 1 (Eta < 0)
	VnCorrected_x_1 = VnRaw_x_1 - hV2DetX_TkEff_1->GetBinContent(icent+1);
	VnCorrected_y_1 = VnRaw_y_1 - hV2DetY_TkEff_1->GetBinContent(icent+1);
                
	hVnSub1[2][icent]->Fill( TMath::Sqrt( VnCorrected_x_1 * VnCorrected_x_1 + VnCorrected_y_1 *VnCorrected_y_1 ) );
	hVn2Dsub1[2][icent]->Fill(VnCorrected_x_1, VnCorrected_y_1);
                
	//-- SubEvt Difference
	hVn2D0v1[2][icent]->Fill(VnCorrected_x_0 - VnCorrected_x_1, VnCorrected_y_0 - VnCorrected_y_1);
        
    }
        
    fHists->Write();
    
    cout<<"File written, process completed"<<endl;
    
}
    
    
    
    
    
    
    
    
    
    
    
