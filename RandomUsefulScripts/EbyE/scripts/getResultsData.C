#include "/Users/jcastle/root/macros/tdrstyle.C"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TGraph.h"
#include <iostream>

bool iter256 = false;

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

static const int NCENT = 20;
static const double cent_min[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95};
static const double cent_max[] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};

void getResultsData(){
    
    setTDRStyle();
    TH1D::SetDefaultSumw2();
    
    TFile * f_MC = new TFile("CastleEbyE.root");
    TFile * f_Unfold = new TFile("data2.root");
    
    TLatex * latex[NCENT];

    //-- 2SE Difference Tilt Classifier
    TF1 * SEDiffTiltTrend[NCENT];
    
    //-- Gaussian Fits for 2SE
    TF1 * fitx[NCENT];
    TF1 * fity[NCENT];
    
    //-- Analyzer Hists
    TH2D * h2D_Full[NCENT];
    TH1D * h1D_Full[NCENT];
    TH2D * h2D_Sub0[NCENT];
    TH2D * h2D_Sub1[NCENT];
    TH2D * h2D_Sub[NCENT];
    TH1D * h1D_SubX[NCENT];
    TH1D * h1D_SubY[NCENT];
    
    //-- Unfolding Hists
    TH2D * hresp_Obs[NCENT];
    TH1D * hUnfold1[NCENT];
    TH1D * hUnfold2[NCENT];
    TH1D * hUnfold4[NCENT];
    TH1D * hUnfold8[NCENT];
    TH1D * hUnfold16[NCENT];
    TH1D * hUnfold32[NCENT];
    TH1D * hUnfold64[NCENT];
    TH1D * hUnfold128[NCENT];
    if(iter256) TH1D * hUnfold256[NCENT];
    
    //-- Ratio Hists
    TH1D * hratio2_1[NCENT];
    TH1D * hratio4_2[NCENT];
    TH1D * hratio8_4[NCENT];
    TH1D * hratio16_8[NCENT];
    TH1D * hratio32_16[NCENT];
    TH1D * hratio64_32[NCENT];
    TH1D * hratio128_64[NCENT];
    
    TH1D * hratio1_128[NCENT];
    TH1D * hratio2_128[NCENT];
    TH1D * hratio4_128[NCENT];
    TH1D * hratio8_128[NCENT];
    TH1D * hratio16_128[NCENT];
    TH1D * hratio32_128[NCENT];
    TH1D * hratio64_128[NCENT];
    
    TCanvas * c0[NCENT];
    TCanvas * c1[NCENT];
    TCanvas * c2[NCENT];
    TCanvas * c3[NCENT];
    TCanvas * c4[NCENT];
    TCanvas * c5[NCENT];
    TCanvas * c6[NCENT];
    TCanvas * c7[NCENT];
    TCanvas * c8[NCENT];
    TCanvas * c9[NCENT];
    TCanvas * c10[NCENT];
    TCanvas * c11[NCENT];
    TCanvas * c12[NCENT];
    TCanvas * c13[NCENT];
    TCanvas * c15[NCENT];
    TCanvas * c16[NCENT];
    
    TLegend * leg2[NCENT];
    TLegend * leg10[NCENT];
    TLine * line[NCENT];
    
    TPaveText * tp_cmspreliminary[NCENT];
    
    
    for(int icent = 0; icent < NCENT; icent++){
    
        if(icent != centbin) continue;
    
        latex[icent] = new TLatex();
    
	//-- 2SE Difference Tilt Classifier
        SEDiffTiltTrend[icent] = new TF1(Form("SEDiffTiltTrend_%i",icent),"[0]*x",-0.6,0.6);
        SEDiffTiltTrend[icent]->SetLineColor(kRed);
        SEDiffTiltTrend[icent]->SetParameter(0,1.);

        //-- Gaussian Fits for 2SE
        fitx[icent] = new TF1(Form("fitx_%i",icent),"gaus",-0.6,0.6);
        fity[icent] = new TF1(Form("fity_%i",icent),"gaus",-0.6,0.6);
    
        //-- Analyzer Hists
        h2D_Full[icent] = (TH2D*) f_Unfold->Get(Form("h2D_2_%i",icent));
        h1D_Full[icent] = (TH1D*) f_Unfold->Get(Form("h1D_2_%i",icent));
        h2D_Sub[icent]  = (TH2D*) f_Unfold->Get(Form("h2Dx_2_%i",icent));
        h1D_SubX[icent] = (TH1D*) f_Unfold->Get(Form("h1Dx_2_%i",icent));
        h1D_SubY[icent] = (TH1D*) f_Unfold->Get(Form("h1Dy_2_%i",icent));
        h2D_Sub0[icent] = (TH2D*) f_Unfold->Get(Form("h2Dsub0_2_%i",icent));
        h2D_Sub1[icent] = (TH2D*) f_Unfold->Get(Form("h2Dsub1_2_%i",icent));
        
        h2D_Sub0[icent]->GetXaxis()->SetTitle("v_{2,x}^{obs,a}");
        h2D_Sub1[icent]->GetXaxis()->SetTitle("v_{2,x}^{obs,b}");
        
        h2D_Sub0[icent]->GetYaxis()->SetTitle("v_{2,y}^{obs,a}");
        h2D_Sub1[icent]->GetYaxis()->SetTitle("v_{2,y}^{obs,b}");
    
        //-- Unfolding Hists
        hresp_Obs[icent]  = (TH2D*) f_Unfold->Get(Form("hresp_2_%i",icent));
        hUnfold1[icent]   = (TH1D*) f_Unfold->Get(Form("hreco1_2_%i",icent));
        hUnfold2[icent]   = (TH1D*) f_Unfold->Get(Form("hreco2_2_%i",icent));
        hUnfold4[icent]   = (TH1D*) f_Unfold->Get(Form("hreco4_2_%i",icent));
        hUnfold8[icent]   = (TH1D*) f_Unfold->Get(Form("hreco8_2_%i",icent));
        hUnfold16[icent]  = (TH1D*) f_Unfold->Get(Form("hreco16_2_%i",icent));
        hUnfold32[icent]  = (TH1D*) f_Unfold->Get(Form("hreco32_2_%i",icent));
        hUnfold64[icent]  = (TH1D*) f_Unfold->Get(Form("hreco64_2_%i",icent));
        hUnfold128[icent] = (TH1D*) f_Unfold->Get(Form("hreco128_2_%i",icent));
        if(iter256) hUnfold256[icent] = (TH1D*) f_Unfold->Get(Form("hreco256_2_%i",icent));

	//-- Normalize unfolded histograms    
        h1D_Full[icent]->Scale(1./h1D_Full[icent]->Integral());
        hUnfold1[icent]->Scale(1./hUnfold1[icent]->Integral());
        hUnfold2[icent]->Scale(1./hUnfold2[icent]->Integral());
        hUnfold4[icent]->Scale(1./hUnfold4[icent]->Integral());
        hUnfold8[icent]->Scale(1./hUnfold8[icent]->Integral());
        hUnfold16[icent]->Scale(1./hUnfold16[icent]->Integral());
        hUnfold32[icent]->Scale(1./hUnfold32[icent]->Integral());
        hUnfold64[icent]->Scale(1./hUnfold64[icent]->Integral());
        hUnfold128[icent]->Scale(1./hUnfold128[icent]->Integral());
        if(iter256) hUnfold256[icent]->Scale(1./hUnfold256[icent]->Integral());
    
        //-- Ratio Hists
        hratio2_1[icent] = (TH1D*) hUnfold2[icent]->Clone();
        hratio4_2[icent] = (TH1D*) hUnfold4[icent]->Clone();
        hratio8_4[icent] = (TH1D*) hUnfold8[icent]->Clone();
        hratio16_8[icent] = (TH1D*) hUnfold16[icent]->Clone();
        hratio32_16[icent] = (TH1D*) hUnfold32[icent]->Clone();
        hratio64_32[icent] = (TH1D*) hUnfold64[icent]->Clone();
        hratio128_64[icent] = (TH1D*) hUnfold128[icent]->Clone();
    
        hratio1_128[icent] = (TH1D*) hUnfold1[icent]->Clone();
        hratio2_128[icent] = (TH1D*) hUnfold2[icent]->Clone();
        hratio4_128[icent] = (TH1D*) hUnfold4[icent]->Clone();
        hratio8_128[icent] = (TH1D*) hUnfold8[icent]->Clone();
        hratio16_128[icent] = (TH1D*) hUnfold16[icent]->Clone();
        hratio32_128[icent] = (TH1D*) hUnfold32[icent]->Clone();
        hratio64_128[icent] = (TH1D*) hUnfold64[icent]->Clone();
    
        hratio2_1[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 2}/v_{2}^{iter = 1}");
        hratio4_2[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 4}/v_{2}^{iter = 2}");
        hratio8_4[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 8}/v_{2}^{iter = 4}");
        hratio16_8[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 16}/v_{2}^{iter = 8}");
        hratio32_16[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 32}/v_{2}^{iter = 16}");
        hratio64_32[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 64}/v_{2}^{iter = 32}");
        hratio128_64[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 128}/v_{2}^{iter = 64}");
    
        hratio1_128[icent]->GetXaxis()->SetTitle("v_{2}");
        //hratio2_128[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 2}/v_{2}^{iter = 128}");
        //hratio4_128[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 4}/v_{2}^{iter = 128}");
        //hratio8_128[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 8}/v_{2}^{iter = 128}");
        //hratio16_128[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 16}/v_{2}^{iter = 128}");
        //hratio32_128[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 32}/v_{2}^{iter = 128}");
        //hratio64_128[icent]->GetXaxis()->SetTitle("v_{2}^{iter = 64}/v_{2}^{iter = 128}");
    
        hratio2_1[icent]->Divide(hUnfold1[icent]);
        hratio4_2[icent]->Divide(hUnfold2[icent]);
        hratio8_4[icent]->Divide(hUnfold4[icent]);
        hratio16_8[icent]->Divide(hUnfold8[icent]);
        hratio32_16[icent]->Divide(hUnfold16[icent]);
        hratio64_32[icent]->Divide(hUnfold32[icent]);
        hratio128_64[icent]->Divide(hUnfold64[icent]);
    
        hratio1_128[icent]->Divide(hUnfold128[icent]);
        hratio2_128[icent]->Divide(hUnfold128[icent]);
        hratio4_128[icent]->Divide(hUnfold128[icent]);
        hratio8_128[icent]->Divide(hUnfold128[icent]);
        hratio16_128[icent]->Divide(hUnfold128[icent]);
        hratio32_128[icent]->Divide(hUnfold128[icent]);
        hratio64_128[icent]->Divide(hUnfold128[icent]);
    
        hratio1_128[icent]->SetMarkerColor(kOrange-2);
        hratio2_128[icent]->SetMarkerColor(kGreen+3);
        hratio4_128[icent]->SetMarkerColor(kCyan);
        hratio8_128[icent]->SetMarkerColor(kMagenta);
        hratio16_128[icent]->SetMarkerColor(kViolet-1);
        hratio32_128[icent]->SetMarkerColor(kBlue);
        hratio64_128[icent]->SetMarkerColor(kRed);
    
        hratio1_128[icent]->SetMarkerStyle(21);
        hratio2_128[icent]->SetMarkerStyle(26);
        hratio4_128[icent]->SetMarkerStyle(32);
        hratio8_128[icent]->SetMarkerStyle(4);
        hratio16_128[icent]->SetMarkerStyle(25);
        hratio32_128[icent]->SetMarkerStyle(27);
        hratio64_128[icent]->SetMarkerStyle(28);
    
        hratio1_128[icent]->SetLineColor(kOrange-2);
        hratio2_128[icent]->SetLineColor(kGreen+3);
        hratio4_128[icent]->SetLineColor(kCyan);
        hratio8_128[icent]->SetLineColor(kMagenta);
        hratio16_128[icent]->SetLineColor(kViolet-1);
        hratio32_128[icent]->SetLineColor(kBlue);
        hratio64_128[icent]->SetLineColor(kRed);
    
        hratio1_128[icent]->GetXaxis()->SetRange(1,63);
        hratio2_128[icent]->GetXaxis()->SetRange(1,63);
        hratio4_128[icent]->GetXaxis()->SetRange(1,63);
        hratio8_128[icent]->GetXaxis()->SetRange(1,63);
        hratio16_128[icent]->GetXaxis()->SetRange(1,63);
        hratio32_128[icent]->GetXaxis()->SetRange(1,63);
        hratio64_128[icent]->GetXaxis()->SetRange(1,63);
    
        hratio1_128[icent]->SetMinimum(0.8);
        hratio1_128[icent]->SetMaximum(1.2);
    
        hUnfold1[icent]->SetMarkerColor(kOrange-2);
        hUnfold2[icent]->SetMarkerColor(kGreen+3);
        hUnfold4[icent]->SetMarkerColor(kCyan);
        hUnfold8[icent]->SetMarkerColor(kMagenta);
        hUnfold16[icent]->SetMarkerColor(kViolet-1);
        hUnfold32[icent]->SetMarkerColor(kBlue);
        hUnfold64[icent]->SetMarkerColor(kRed);
        hUnfold128[icent]->SetMarkerColor(kBlack);
        if(iter256) hUnfold256[icent]->SetMarkerColor(kCyan+2);
     
        hUnfold1[icent]->SetLineColor(kOrange-2);
        hUnfold2[icent]->SetLineColor(kGreen+3);
        hUnfold4[icent]->SetLineColor(kCyan);
        hUnfold8[icent]->SetLineColor(kMagenta);
        hUnfold16[icent]->SetLineColor(kViolet-1);
        hUnfold32[icent]->SetLineColor(kBlue);
        hUnfold64[icent]->SetLineColor(kRed);
        hUnfold128[icent]->SetLineColor(kBlack);
        if(iter256) hUnfold256[icent]->SetLineColor(kCyan+2);
    
        hUnfold1[icent]->SetMarkerStyle(21);
        hUnfold2[icent]->SetMarkerStyle(26);
        hUnfold4[icent]->SetMarkerStyle(32);
        hUnfold8[icent]->SetMarkerStyle(4);
        hUnfold16[icent]->SetMarkerStyle(25);
        hUnfold32[icent]->SetMarkerStyle(27);
        hUnfold64[icent]->SetMarkerStyle(28);
        hUnfold128[icent]->SetMarkerStyle(30);
        hUnfold128[icent]->SetMarkerStyle(8);
    
        hUnfold1[icent]->GetXaxis()->SetRange(1,2*63);
        hUnfold2[icent]->GetXaxis()->SetRange(1,2*63);
        hUnfold4[icent]->GetXaxis()->SetRange(1,2*63);
        hUnfold8[icent]->GetXaxis()->SetRange(1,2*63);
        hUnfold16[icent]->GetXaxis()->SetRange(1,2*63);
        hUnfold32[icent]->GetXaxis()->SetRange(1,2*63);
        hUnfold64[icent]->GetXaxis()->SetRange(1,2*63);
        hUnfold128[icent]->GetXaxis()->SetRange(1,2*63);
        if(iter256) hUnfold256[icent]->GetXaxis()->SetRange(1,2*63);
        h1D_Full[icent]->GetXaxis()->SetRange(1,2*63);
    
    //====================================================================
        
        tp_cmspreliminary[icent] = new TPaveText(0.65,0.75,0.94,0.94,"NDC");
        tp_cmspreliminary[icent]->SetFillColor(0);
        tp_cmspreliminary[icent]->AddText("CMS Preliminary");
        tp_cmspreliminary[icent]->AddText("2011 PbPb, #sqrt{s} = 2.76 TeV");
        tp_cmspreliminary[icent]->AddText("p_{T} > 0.5 GeV/c, |#eta| < 2.4");
        tp_cmspreliminary[icent]->AddText(Form("%.0f - %.0f %% Centrality",cent_min[icent],cent_max[icent]));
        
        
        //-- 2SE Fits
        c0[icent] = new TCanvas(Form("c0_%i",icent),Form("c0_%i",icent),600,600);
        c0[icent]->cd();
        c0[icent]->SetLogy();
        h1D_SubX[icent]->Fit(Form("fitx_%i",icent));
        double sigma_x = fitx[icent]->GetParameter("Sigma");
        double sigma_xe = fitx[icent]->GetParError(fitx[icent]->GetParNumber("Sigma"));
        double mean_x  = fitx[icent]->GetParameter("Mean");
        double chi2_x  = fitx[icent]->GetChisquare();
        double ndf_x = fitx[icent]->GetNDF();
        h1D_SubX[icent]->Draw();
        latex[icent]->DrawLatex(-0.2, 500, Form("#splitline{#chi^{2}/NDF = %f}{#splitline{#delta_{2SE} = %f}{#pm%f}}", chi2_x/ndf_x, sigma_x, sigma_xe));
        latex[icent]->DrawLatex(-0.2, 50, Form("Mean = %f", mean_x));
        //tp_cmspreliminary[icent]->Draw();
        c0[icent]->SaveAs(Form("plots/SubXFit_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c0[icent]->Close();
    
        c1[icent] = new TCanvas(Form("c1_%i",icent),Form("c1_%i",icent),600,600);
        c1[icent]->cd();
        c1[icent]->SetLogy();
        h1D_SubY[icent]->Fit(Form("fity_%i",icent));
        double sigma_y = fity[icent]->GetParameter("Sigma");
        double sigma_ye = fity[icent]->GetParError(fity[icent]->GetParNumber("Sigma"));
        double mean_y  = fity[icent]->GetParameter("Mean");
        double chi2_y  = fity[icent]->GetChisquare();
        double ndf_y = fity[icent]->GetNDF();
        h1D_SubY[icent]->Draw();
        latex[icent]->DrawLatex(-0.2, 500, Form("#splitline{#chi^{2}/NDF = %f}{#splitline{#delta_{2SE} = %f}{#pm%f}}", chi2_y/ndf_y, sigma_y, sigma_ye));
        latex[icent]->DrawLatex(-0.2, 50, Form("Mean = %f", mean_y));
        //tp_cmspreliminary[icent]->Draw();
        c1[icent]->SaveAs(Form("plots/SubYFit_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c1[icent]->Close();

    
    
        //-- Unfolding Results
        c2[icent] = new TCanvas(Form("c2_%i",icent),Form("c2_%i",icent),600,600);
        c2[icent]->cd();
        c2[icent]->SetLogy();
        h1D_Full[icent]->Draw();
        hUnfold1[icent]->Draw("same");
        hUnfold2[icent]->Draw("same");
        hUnfold4[icent]->Draw("same");
        hUnfold8[icent]->Draw("same");
        hUnfold16[icent]->Draw("same");
        hUnfold32[icent]->Draw("same");
        hUnfold64[icent]->Draw("same");
        hUnfold128[icent]->Draw("same");
        leg2[icent] = new TLegend(0.3, 0.35, 0.6, 0.65);
        leg2[icent]->SetFillColor(kWhite);
        leg2[icent]->SetBorderSize(0);
        leg2[icent]->AddEntry(h1D_Full[icent],"Measured","pl");
        leg2[icent]->AddEntry(hUnfold1[icent],"Bayes 1","pl");
        leg2[icent]->AddEntry(hUnfold2[icent],"Bayes 2","pl");
        leg2[icent]->AddEntry(hUnfold4[icent],"Bayes 4","pl");
        leg2[icent]->AddEntry(hUnfold8[icent],"Bayes 8","pl");
        leg2[icent]->AddEntry(hUnfold16[icent],"Bayes 16","pl");
        leg2[icent]->AddEntry(hUnfold32[icent],"Bayes 32","pl");
        leg2[icent]->AddEntry(hUnfold64[icent],"Bayes 64","pl");
        leg2[icent]->AddEntry(hUnfold128[icent],"Bayes 128","pl");
        leg2[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c2[icent]->SaveAs(Form("plots/unfolding_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c2[icent]->Close();
    

    
        //-- Ratio Plots
        c3[icent] = new TCanvas(Form("c3_%i",icent),Form("c3_%i",icent),600,600);
        c3[icent]->cd();
        hratio2_1[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c3[icent]->SaveAs(Form("plots/ratio2_1_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c3[icent]->Close();
    
        c4[icent] = new TCanvas(Form("c4_%i",icent),Form("c4_%i",icent),600,600);
        c4[icent]->cd();
        hratio4_2[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c4[icent]->SaveAs(Form("plots/ratio4_2_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c4[icent]->Close();
    
        c5[icent] = new TCanvas(Form("c5_%i",icent),Form("c5_%i",icent),600,600);
        c5[icent]->cd();
        hratio8_4[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c5[icent]->SaveAs(Form("plots/ratio8_4_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c5[icent]->Close();
    
        c6[icent] = new TCanvas(Form("c6_%i",icent),Form("c6_%i",icent),600,600);
        c6[icent]->cd();
        hratio16_8[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c6[icent]->SaveAs(Form("plots/ratio16_8_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c6[icent]->Close();
    
        c7[icent] = new TCanvas(Form("c7_%i",icent),Form("c7_%i",icent),600,600);
        c7[icent]->cd();
        hratio32_16[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c7[icent]->SaveAs(Form("plots/ratio32_16_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c7[icent]->Close();
    
        c8[icent] = new TCanvas(Form("c8_%i",icent),Form("c8_%i",icent),600,600);
        c8[icent]->cd();
        hratio64_32[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c8[icent]->SaveAs(Form("plots/ratio64_32_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c8[icent]->Close();
    
        c9[icent] = new TCanvas(Form("c9_%i",icent),Form("c9_%i",icent),600,600);
        c9[icent]->cd();
        hratio128_64[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c9[icent]->SaveAs(Form("plots/ratio128_64_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c9[icent]->Close();
    
	//-- Different ratios
        c10[icent] = new TCanvas(Form("c10_%i",icent),Form("c10_%i",icent),600,600);
        c10[icent]->cd();
        hratio1_128[icent]->Draw();
        hratio2_128[icent]->Draw("same");
        hratio4_128[icent]->Draw("same");
        hratio8_128[icent]->Draw("same");
        hratio16_128[icent]->Draw("same");
        hratio32_128[icent]->Draw("same");
        hratio64_128[icent]->Draw("same");
        leg10[icent] = new TLegend(0.42, 0.65, 0.57, 0.9);
        leg10[icent]->SetFillColor(kWhite);
        leg10[icent]->SetBorderSize(0);
        leg10[icent]->AddEntry(hratio1_128[icent],"v_{2}^{iter = 1}/v_{2}^{iter = 128}","pl");
        leg10[icent]->AddEntry(hratio2_128[icent],"v_{2}^{iter = 2}/v_{2}^{iter = 128}","pl");
        leg10[icent]->AddEntry(hratio4_128[icent],"v_{2}^{iter = 4}/v_{2}^{iter = 128}","pl");
        leg10[icent]->AddEntry(hratio8_128[icent],"v_{2}^{iter = 8}/v_{2}^{iter = 128}","pl");
        leg10[icent]->AddEntry(hratio16_128[icent],"v_{2}^{iter = 16}/v_{2}^{iter = 128}","pl");
        leg10[icent]->AddEntry(hratio32_128[icent],"v_{2}^{iter = 32}/v_{2}^{iter = 128}","pl");
        leg10[icent]->AddEntry(hratio64_128[icent],"v_{2}^{iter = 64}/v_{2}^{iter = 128}","pl");
        leg10[icent]->Draw();
        line[icent] = new TLine(0,1,0.25,1);
        line[icent]->SetLineColor(kBlue+4);
        line[icent]->SetLineStyle(9);
        line[icent]->SetLineWidth(3.5);
        line[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c10[icent]->SaveAs(Form("plots/ratioTo128_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c10[icent]->Close();
    
    
    
	//-- Observed flow vectors
        c11[icent] = new TCanvas(Form("c11_%i",icent),Form("c11_%i",icent),600,600);
        c11[icent]->cd();
        h2D_Full[icent]->Draw("colz");
        //tp_cmspreliminary[icent]->Draw();
        c11[icent]->SaveAs(Form("plots/h2D_Full_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c11[icent]->Close();
    
        c12[icent] = new TCanvas(Form("c12_%i",icent),Form("c12_%i",icent),600,600);
        c12[icent]->cd();
        h1D_Full[icent]->Draw();
        //tp_cmspreliminary[icent]->Draw();
        c12[icent]->SaveAs(Form("plots/h1D_Full_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c12[icent]->Close();
    
        c13[icent] = new TCanvas(Form("c13_%i",icent),Form("c13_%i",icent),600,600);
        c13[icent]->cd();
	//h2D_Sub[icent]->Fit(Form("SEDiffTiltTrend_%i",icent));
        h2D_Sub[icent]->Draw("colz");
        //tp_cmspreliminary[icent]->Draw();
	//double tiltSlope = SEDiffTiltTrend[icent]->GetParameter(0);
        //double tiltAngle = (180./TMath::Pi()) * TMath::ATan2(tiltSlope,1);
	//latex[icent]->DrawLatex(-0.55,-0.55,Form("Tilt angle = %.2f#circ",tiltAngle));
        c13[icent]->SaveAs(Form("plots/h2D_Sub_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c13[icent]->Close();
        
        c15[icent] = new TCanvas(Form("c15_%i",icent),Form("c15_%i",icent),600,600);
        c15[icent]->cd();
        h2D_Sub0[icent]->Draw("colz");
        //tp_cmspreliminary[icent]->Draw();
        c15[icent]->SaveAs(Form("plots/h2D_Sub0_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c15[icent]->Close();
        
        c16[icent] = new TCanvas(Form("c16_%i",icent),Form("c16_%i",icent),600,600);
        c16[icent]->cd();
        h2D_Sub1[icent]->Draw("colz");
        //tp_cmspreliminary[icent]->Draw();
        c16[icent]->SaveAs(Form("plots/h2D_Sub1_cent_%.0f_%.0f.png",cent_min[icent],cent_max[icent]));
        c16[icent]->Close();

    
    }
    
    //-- ATLAS Compare Plot 20 - 25 %
    double ATLASv2[27] = {0.004, 
			  0.0119,
			  0.0199,
			  0.0279,
			  0.0358,
			  0.0438,
			  0.0518,
			  0.0597,
			  0.0677,
			  0.0757,
			  0.0836,
			  0.0916,
			  0.0996,
			  0.1075,
			  0.1155,
			  0.1235,
			  0.1314,
			  0.1394,
			  0.1473,
			  0.1553,
			  0.1633,
			  0.1712,
			  0.1792,
			  0.1872,
			  0.1951,
			  0.2031,
			  0.2111
    };
    double ATLASPv2[27] = {0.1396,
			   0.4543,
			   0.8843,
			   1.4561,
			   2.2996,
			   3.4501,
			   4.7409,
			   6.2473,
			   7.8973,
			   9.464,
			   10.7785,
			   11.6046,
			   11.8612,
			   11.4698,
			   10.3253,
			   8.968,
			   7.3163,
			   5.6882,
			   4.11,
			   2.724,
			   1.689,
			   0.967,
			   0.5206,
			   0.2637,
			   0.1324,
			   0.0575,
			   0.0168
    };
    
    double ATLASPv2ERR[27] = {0.012465151,
			      0.026816413,
			      0.037190321,
			      0.048488349,
			      0.069370383,
			      0.100239563,
			      0.1325581,
			      0.168334102,
			      0.205278835,
			      0.23712347,
			      0.260353951,
			      0.269828464,
			      0.265256706,
			      0.246495801,
			      0.213249267,
			      0.199405366,
			      0.192192091,
			      0.17253985,
			      0.19056697,
			      0.153966522,
			      0.109862505,
			      0.071437875,
			      0.04370492,
			      0.055932817,
			      0.04552373,
			      0.03212787,
			      0.014603082
    };
    
    TH1D * hATLAS = new TH1D("hATLAS","hATLAS",150,0,0.25);
    hATLAS->GetXaxis()->SetTitle("v_{2}");
    for(int i=0;i<27;i++){
        hATLAS->Fill(ATLASv2[i],ATLASPv2[i]);
        hATLAS->SetBinError(hATLAS->FindBin(ATLASv2[i]),ATLASPv2ERR[i]);
    }

    hATLAS->Scale(1./hATLAS->Integral());
    hATLAS->GetYaxis()->SetTitle("P(v_{2})");
    
    hATLAS->SetLineColor(kRed);
    hATLAS->SetMarkerColor(kRed);
    hATLAS->SetMarkerStyle(4);

    //-- Observed v2 distribution for reference
    h1D_Full[centbin]->SetLineColor(kBlue);
    h1D_Full[centbin]->SetMarkerColor(kBlue);
    h1D_Full[centbin]->SetMarkerStyle(25);    
    h1D_Full[centbin]->Scale(hATLAS->GetMaximum()/h1D_Full[centbin]->GetMaximum());

    //-- Compare 128th iteration of unfolding to ATLAS
    hUnfold128[centbin]->Scale(hATLAS->GetMaximum()/hUnfold128[centbin]->GetMaximum());
    if(iter256) hUnfold256[centbin]->Scale(hATLAS->GetMaximum()/hUnfold256[centbin]->GetMaximum());
    
    TCanvas * cATLAS = new TCanvas("cATLAS","cATLAS",600,600);
    cATLAS->cd();
    cATLAS->SetLogy();
    hATLAS->Draw();
    h1D_Full[centbin]->Draw("same");
    if(iter256) hUnfold256[centbin]]->Draw("same");
    else hUnfold128[centbin]->Draw("same");
    TLegend * legATLAS = new TLegend(0.3, 0.35, 0.6, 0.65);
    legATLAS->SetFillColor(kWhite);
    legATLAS->SetBorderSize(0);
    legATLAS->AddEntry(h1D_Full[centbin],"CMS v_{2}^{obs}","pl");
    if(iter256) legATLAS->AddEntry(hUnfold256[centbin],"CMS iter = 256","pl");
    else legATLAS->AddEntry(hUnfold128[centbin],"CMS iter = 128","pl");
    legATLAS->AddEntry(hATLAS,"ATLAS P(v_{2})","pl");
    legATLAS->Draw();
    cATLAS->SaveAs("plots/ATLASCompare.png");
    cATLAS->Close();

//-- Bessel-Gaussian fits to ATLAS and CMS
TF1 * ATLASbessGauss = new TF1("ATLASbessGauss","[0]*(x/([2]*[2]))*TMath::Exp(-(x*x+[1]*[1])/(2*[2]*[2]))*TMath::BesselI0(([1]*x)/([2]*[2]))",0.00,0.22);
TF1 * CMSbessGauss = new TF1("CMSbessGauss","[0]*(x/([2]*[2]))*TMath::Exp(-(x*x+[1]*[1])/(2*[2]*[2]))*TMath::BesselI0(([1]*x)/([2]*[2]))",0.00,0.22);
ATLASbessGauss->SetParameters(1,0.09,0.03);
CMSbessGauss->SetParameters(1,0.09,0.03);
ATLASbessGauss->SetLineColor(kRed);
CMSbessGauss->SetLineColor(kBlack);
TCanvas * cFit = new TCanvas("cFit","cFit",600,600);
cFit->SetLogy();
cFit->cd();
hATLAS->Fit("ATLASbessGauss","N","N",0,0.22);
hUnfold128[centbin]->Fit("CMSbessGauss","N","N",0,0.22);
ATLASbessGauss->Draw();
CMSbessGauss->Draw("same");
legATLAS->Draw();    
cFit->SaveAs("plots/BGFits.png");    
cFit->Close();

//-- Bessel-Gaussian fit comparison to ATLAS and CMS unfolded distns 
TCanvas * cComp = new TCanvas("cComp","cComp",600,600);
cComp->cd();
TH1D * hComp = new TH1D("hComp","hComp",2,1,2);
hComp->SetMinimum(0.02);
hComp->SetMaximum(0.1);
hComp->GetXaxis()->SetBinLabel(1,"<v_{2}>");
hComp->GetXaxis()->SetBinLabel(2,"#sigma_{v_{2}}");
hComp->SetBinContent(1,CMSbessGauss->GetParameter(1));
hComp->SetBinContent(2,CMSbessGauss->GetParameter(2));
hComp->SetBinError(1,CMSbessGauss->GetParError(1));
hComp->SetBinError(2,CMSbessGauss->GetParError(2));
hComp->Draw();
TBox * ATLASMu = new TBox(1.0, ATLASbessGauss->GetParameter(1) - 0.019 * ATLASbessGauss->GetParameter(1), 1.5, ATLASbessGauss->GetParameter(1) + 0.019 * ATLASbessGauss->GetParameter(1));
ATLASMu->SetFillStyle(3001);
ATLASMu->SetFillColor(2);
ATLASMu->Draw("same");
TBox * ATLASSig = new TBox(1.5, ATLASbessGauss->GetParameter(2) - 0.019 * ATLASbessGauss->GetParameter(2), 2., ATLASbessGauss->GetParameter(2) + 0.019 * ATLASbessGauss->GetParameter(2));
ATLASSig->SetFillStyle(3001);
ATLASSig->SetFillColor(2);
ATLASSig->Draw("same");
TLegend * legcomp = new TLegend(0.2, 0.52, 0.47, 0.65);
legcomp->SetFillColor(kWhite);
//legcomp->SetBorderSize(0);
legcomp->AddEntry(hUnfold128[centbin],"CMS iter = 128","p");
legcomp->AddEntry(hATLAS,"ATLAS P(v_{2})","l");
legcomp->Draw();
cComp->SaveAs("plots/FitComp.png");
cComp->Close();
    
}
