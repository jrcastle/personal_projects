#include "/home/j550c590/tdrstyle.C"
#include "DAgostiniUnfold.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/interface/differentailBinning.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfoldBayes.h"
#include "/home/j550c590/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/macros/unfolding/RooUnfold/src/RooUnfold.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h" 
#include "TFile.h"
#include "TCanvas.h"


void testDAgostiniUnfold(){

  setTDRStyle();

  int ptBin   = 0;
  int centBin = 1;
  int VN      = 2;

  TFile * fMC   = new TFile( Form("/rfs/jcastle/PbPb2015/mimicDataMC/cent%i-%i/pt%.2f-%.2f/MCTree.root", cent_min[centBin], cent_max[centBin], pt_min[ptBin], pt_max[ptBin]) );
  TFile * fData = new TFile( Form("/rfs/jcastle/PbPb2015/mimicDataMC/cent%i-%i/pt%.2f-%.2f/CastleEbyE.root", cent_min[centBin], cent_max[centBin], pt_min[ptBin], pt_max[ptBin]) );

  //-- CastleEbyE.root data members
  static const int NVn = 7;
  TH1D * hVnFull[NVn][NCENT_PBPB];
  TH2D * hResp[NVn][NCENT_PBPB];

  //-- Unfolding data
  TH2D * hresp[NCENT];
  TH1D * h1D[NCENT];
  TH1D * hreco1[NCENT];
  TH1D * hreco1_RU[NCENT];
  RooUnfoldResponse * RUresp[NCENT];
  RooUnfoldBayes * RUBayes[NCENT];


  //-- Get DATA
  for(int ivn = 1; ivn <=4; ivn++){
    if(ivn != VN) continue;
    for(int icent = 0; icent < NCENT_PBPB; icent++){
      if(icent == NCENT40_cutoff) break;
      hVnFull[ivn][icent]   = (TH1D*) fData->Get(Form("qwebye/hVnFull_%i_%i",ivn,icent));
      hResp[ivn][icent]     = (TH2D*) fData->Get(Form("qwebye/hResp_%i_%i",ivn,icent));
    }
  }

  //-- Condense from 2.5 to 10% cent binning
  for(int c = 0; c < NCENT; c++){

    hreco1[c]    = 0;
    hreco1_RU[c] = 0;

    if(c != centBin) continue;
    std::cout<<"!! Processing Cent = "<<c<<std::endl;


    //-- Observed distribution
    h1D[c] = (TH1D*)hVnFull[VN][4*c]->Clone(Form("h1D_%i_%i", VN, c));
    h1D[c]->Add(hVnFull[VN][4*c+1]);
    h1D[c]->Add(hVnFull[VN][4*c+2]);
    h1D[c]->Add(hVnFull[VN][4*c+3]);

    //-- Response matrix
    hresp[c] = (TH2D*) hResp[VN][4*c]->Clone(Form("hresp_%i_%i", VN, c));
    hresp[c]->SetOption("colz");
    hresp[c]->Add(hResp[VN][4*c+1]);
    hresp[c]->Add(hResp[VN][4*c+2]);
    hresp[c]->Add(hResp[VN][4*c+3]);

    DAgostiniUnfold unfolder(hresp[c]);
    unfolder.Unfold(h1D[c], 1);
    hreco1[c] = (TH1D*) unfolder.GetHReco();
    hreco1[c]->GetXaxis()->SetTitle("v_{2}");
    hreco1[c]->GetYaxis()->SetTitle("Events");
    hreco1[c]->SetMarkerStyle(20);


    RUresp[c] = new RooUnfoldResponse(0, 0, hresp[c], Form("response_%i_%i", VN, c));
    RUBayes[c] = new RooUnfoldBayes(RUresp[c], h1D[c], 1);
    hreco1_RU[c] = (TH1D*) RUBayes[c]->Hreco();
    hreco1_RU[c]->SetLineColor(2);
    hreco1_RU[c]->SetMarkerColor(2);
    hreco1_RU[c]->SetMarkerStyle(24);

  }

  TH1D * hratio = (TH1D*) hreco1[centBin]->Clone("ratio");
  hratio->Divide(hreco1_RU[centBin]);
  hratio->GetXaxis()->SetTitle("v_{2}");
  hratio->GetYaxis()->SetTitle("Ratio");
  TLine * line = new TLine(0., 1., v2Max[ptBin][centBin], 1.);
  line->SetLineColor(2);
  line->SetLineStyle(2);
  line->SetLineWidth(2);

  TCanvas * c = new TCanvas("c","c", 500, 1000);
  c->Divide(1,2);
  c->cd(1);
  c->cd(1)->SetLogy();
  hreco1[centBin]->Draw();
  hreco1_RU[centBin]->Draw("same");

  c->cd(2);
  hratio->Draw();
  line->Draw("same");

}
