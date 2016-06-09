#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/interface/differentialBinning.h"
#include "/home/j550c590/tdrstyle.C"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include <iostream>


void makeMCTree() {

  int ptBin   = 0;
  int centBin = 0;
  bool realisticNEvt = 0;

  int N = 1600000;
  //int N = 10000;
  int n0;
  int n1;

  double v2Max = 0.5;
  int NBins    = 125;

  TString fileName = "MCTree.root";

  double nOrder_ = 2.;
  double pi = TMath::Pi();

  TTree * tree;
  TDirectory * qwebyeDir;
  double centval;
  double vtx;
  TH2D * sumw;
  TH2D * sumwqx;
  TH2D * sumwqy;
  TH2I * mult;

  TH1D * hV2True;
  double v2True;
  Double_t etaTrack[500];

  setTDRStyle();
    
  TFile * tf = new TFile(fileName,"recreate");
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  TH2I::SetDefaultSumw2();
    
  tf->cd();
  qwebyeDir = tf->mkdir("ebyeana","ebyeana");
  qwebyeDir->cd();
  tree    = new TTree("tree","EbyE Tree");
  hV2True = new TH1D("hV2True","hV2True", NBins, 0, v2Max);
  sumwqx  = new TH2D("sumwqx","sumwqx",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
  sumwqy  = new TH2D("sumwqy","sumwqy",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
  sumw    = new TH2D("sumw","sumw",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
  mult   = new TH2I("mult","mult",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);

  sumwqx->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqx->GetYaxis()->SetTitle("#eta");
  sumwqx->SetOption("colz");

  sumwqy->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumwqy->GetYaxis()->SetTitle("#eta");
  sumwqy->SetOption("colz");

  sumw->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  sumw->GetYaxis()->SetTitle("#eta");
  sumw->SetOption("colz");

  mult->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  mult->GetYaxis()->SetTitle("#eta");
  mult->SetOption("colz");

  tree->Branch("Cent",   &centval, "cent/D");
  tree->Branch("v2True", &v2True,  "v2True/D");
  tree->Branch("sumwqx", "TH2D",   &sumwqx, 128000, 0);
  tree->Branch("sumwqy", "TH2D",   &sumwqy, 128000, 0);
  tree->Branch("sumw",   "TH2D",   &sumw,   128000, 0);
  tree->Branch("mult",   "TH2I",   &mult,   128000, 0);
    
  TRandom3 * ran = new TRandom3(0);

  TF1 * bessGauss = new TF1("bessGauss","[0]*(x/([2]*[2]))*TMath::Exp(-(x*x+[1]*[1])/(2*[2]*[2]))*TMath::BesselI0(([1]*x)/([2]*[2]))", 0.00, v2Max);
  bessGauss->SetParameters(1, v2BGmcMU[ptBin][centBin], v2BGmcDELTA[ptBin][centBin]);

  TF1 * dNdphi = new TF1("dNdphi","1 + 2 * ( [0] * TMath::Cos( 1 * (x - [6]) ) + [1] * TMath::Cos( 2 * (x - [6]) ) + [2] * TMath::Cos( 3 * (x - [6]) ) + [3] * TMath::Cos( 4 * (x - [6]) ) + [4] * TMath::Cos( 5 * (x - [6]) ) + [5] * TMath::Cos( 6 * (x - [6]) ) )",-pi,pi);

  Double_t v1  = 0.;
  //Double_t v2  = 0.;    
  Double_t v3  = 0.;
  Double_t v4  = 0.;
  Double_t v5  = 0.;
  Double_t v6  = 0.;
  Double_t Psi = 0.;
  Int_t multGauss;

  int nevt;
  if(realisticNEvt) nevt = NEVENTS[ptBin][centBin];
  else              nevt = N;

  for(Int_t i=0; i<nevt; i++){
        
    if((i+1)% 20000 == 0) std::cout<<"Processing Event "<<i+1<<"\t"<<(double)100*(i+1)/nevt<<"% Completed"<<std::endl;
        
    mult->Reset();
    sumw->Reset();
    sumwqx->Reset();
    sumwqy->Reset();

    //-- Set Event Parameters
    multGauss = 0;
    while(multGauss < 4) multGauss = ran->Gaus(meanMult[ptBin][centBin], stdDevMult[ptBin][centBin]);

    v2True = 0.6;
    while(v2True >=0.5) v2True = bessGauss->GetRandom();

    Psi = ran->Uniform(-pi,pi);
    hV2True->Fill(v2True);
    centval = ran->Uniform(cent_min[centBin],cent_max[centBin]);
        
    dNdphi->SetParameters(v1,v2True,v3,v4,v5,v6,Psi);

    n0 = 0;
    n1 = 0;

    while(n0 < 2 || n1 < 2){

      n0 = 0;
      n1 = 0;
      for(int iMult = 0; iMult<multGauss; iMult++){
	etaTrack[iMult] = ran->Uniform(-1.,1.);
	if(etaTrack[iMult] > 0) n0++;
	else                    n1++;
      }
      //if(n0 < 2 || n1 < 2) std::cout<<"Rethrow etaTracks..."<<std::endl;
    }
    
    for(int iMult = 0; iMult<multGauss; iMult++){
            
      Double_t phiTrack = dNdphi->GetRandom();
      Double_t ptTrack = ran->Uniform(pt_min[ptBin], pt_max[ptBin]);
      double w = 1.;
      
      mult->Fill(ptTrack, etaTrack[iMult], 1);    
      sumw->Fill(ptTrack, etaTrack[iMult], 1./w);
      sumwqx->Fill(ptTrack, etaTrack[iMult], TMath::Cos(nOrder_*phiTrack)/w);
      sumwqy->Fill(ptTrack, etaTrack[iMult], TMath::Sin(nOrder_*phiTrack)/w);
      
    } //-- End mult loop

    tree->Fill();

  } //-- End event loop

  tf->Write();
  cout<<"Process Complete"<<endl;

} //-- End macro
