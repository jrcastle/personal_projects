// -*- C++ -*-
//
// Package:    TrackEfficiencyMC
// Class:      TrackEfficiencyMC
//
/*
 
 Description: Interpolates track efficiencies for pt, centrality, phi, and eta
 and spits out an efficiency for each.
 
 Implementation:
 TrackEfficiencyMC * tkEff = new TrackEfficiencyMC("MC");
 
 //Enter Event Loop
 double cent_ =   (however you get centrality);
 double centbin = (convert centrality to a centrality bin MUST INCORPORATE THE NEW 200 CENTRALITY BINS STANDARD, i.e. cent = 50 ==> centbin = 100.)
 
 //Enter Track Loop
 double trackPt_  (however you get pt);
 double trackPhi_ (however you get phi);
 double trackEta_ (however you get Eta);
 
 double eff_cent   = 1.;
 double eff_pt     = 1.;
 double eff_accept = 1.;
 
 double eff = tkEff->getEfficiencies(trackPt_, cent_, trackPhi_, trackEta_, eff_pt, eff_cent, eff_accept);

*/
//
// Original Author:  James Castle
//         Created:  Mon Jan 26 2015
#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TDirectory.h"

/*static const int npt = 29;
double ptmin[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.55, 0.55, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.65, 0.65, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 8.0};
double ptmax[] = {0.55, 0.55, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.65, 0.65, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 3.0, 8.0, 8.0, 8.0, 300};
double cent_min[] = {0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 0};
double cent_max[] = {20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 200, 200};
*/
class TrackEfficiencyMC{
    
    public:
        TrackEfficiencyMC(TString effOption);
        double getEfficiencies(double pt, double cent, double phi, double eta);
    
    private:
        TFile * f_eff;
        TDirectory * effDir[npt];
        TProfile * p_eff_pt[npt];
        TProfile * p_eff_rmin[npt];
        TProfile2D * p_eff_accept[npt];
        TH2D * p_eff_acceptDATA[npt];
        TH1D * minimumCount;
        
    
        double eff;
        double eff_pt;
        double eff_accept;
    
        double eff_pt_low;
        double eff_pt_mid;
        double eff_pt_high;
        double pt_low;
        double pt_mid;
        double pt_high;
    
        double a0;
        double a1;
        double a2;
    
        double phi_low;
        double phi_high;
        double eta_low;
        double eta_high;
    
        int phi_bin_low;
        int phi_bin_high;
        int eta_bin_low;
        int eta_bin_high;
    
        double Q11;
        double Q21;
        double Q12;
        double Q22;
    
};

//-----------------------------------------------------------------------------
//                                  METHODS
//-----------------------------------------------------------------------------

TrackEfficiencyMC::TrackEfficiencyMC(TString effOption){
    
    
    f_eff = new TFile(effOption);
    
    for(int ipt=0; ipt<npt;ipt++){
        
        effDir[ipt]= (TDirectory*) f_eff->Get(Form("eff_pt%d_%d_cent%d_%d",(int)(100*ptmin[ipt]),(int)(100*ptmax[ipt]),(int)(cent_min[ipt]/2),(int)(cent_max[ipt]/2)));
        effDir[ipt]->cd();
        p_eff_pt[ipt]=(TProfile*) effDir[ipt]->Get("p_eff_pt");
        p_eff_accept[ipt]=(TProfile2D*) effDir[ipt]->Get("p_eff_acceptance");
        p_eff_rmin[ipt]=(TProfile*) effDir[ipt]->Get("p_eff_rmin");
        
    }
    
    eff_pt_low = 1;
    eff_pt_mid = 1;
    eff_pt_high = 1;
    pt_low = 1;
    pt_mid = 1;
    pt_high = 1;
    
    a0 = 1;
    a1 = 1;
    a2 = 1;
    
    phi_low = 1;
    phi_high = 1;
    eta_low = 1;
    eta_high = 1;
    
    phi_bin_low = 1;
    phi_bin_high = 1;
    eta_bin_low = 1;
    eta_bin_high = 1;
    
    Q11 = 1;
    Q21 = 1;
    Q12 = 1;
    Q22 = 1;
    
}
//-----------------------------------------------------------------------------

double TrackEfficiencyMC::getEfficiencies(double pt, double cent, double phi, double eta){
    
    for(int ipt=0;ipt<npt;ipt++){
        
        if(pt>=ptmin[ipt] && pt<ptmax[ipt] && cent>=cent_min[ipt] && cent<cent_max[ipt]){
            
            // 3-point interpolate eff_pt
            eff_pt_low = p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt)-1);
            eff_pt_mid = p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
            eff_pt_high = p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt)+1);
            
            pt_low = p_eff_pt[ipt]->GetBinCenter(p_eff_pt[ipt]->FindBin(pt) - 1 );
            pt_mid = p_eff_pt[ipt]->GetBinCenter(p_eff_pt[ipt]->FindBin(pt) );
            pt_high = p_eff_pt[ipt]->GetBinCenter(p_eff_pt[ipt]->FindBin(pt) + 1 );
            
            if(pt_low < 0.5){
                
                eff_pt_low = p_eff_pt[ipt]->GetBinContent(1);
                eff_pt_mid = p_eff_pt[ipt]->GetBinContent(2);
                eff_pt_high = p_eff_pt[ipt]->GetBinContent(3);
                
                pt_low = p_eff_pt[ipt]->GetBinCenter(1);
                pt_mid = p_eff_pt[ipt]->GetBinCenter(2);
                pt_high = p_eff_pt[ipt]->GetBinCenter(3);
                
            }
            
            if(pt_high > 300){
                
                eff_pt_low = p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->GetNbinsX() - 2);
                eff_pt_mid = p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->GetNbinsX() - 1);
                eff_pt_high = p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->GetNbinsX());
                
                pt_low = p_eff_pt[ipt]->GetBinCenter(p_eff_pt[ipt]->GetNbinsX() - 2);
                pt_mid = p_eff_pt[ipt]->GetBinCenter(p_eff_pt[ipt]->GetNbinsX() - 1);
                pt_high = p_eff_pt[ipt]->GetBinCenter(p_eff_pt[ipt]->GetNbinsX());
                
                
            }
            
            if(pt_low < ptmin[ipt]){
                
                for(int jpt=0;jpt<npt;jpt++){
                    
                    if(pt_low>=ptmin[jpt] && pt_low<ptmax[jpt] && cent>=cent_min[jpt] && cent<cent_max[jpt]){
                        
                        eff_pt_low = p_eff_pt[jpt]->GetBinContent(p_eff_pt[jpt]->FindBin(pt_low));
                        
                    }
                    
                }
                
            }
            
            if(pt_high > ptmax[ipt]){
                
                for(int jpt=0;jpt<npt;jpt++){
                    
                    if(pt_high>=ptmin[jpt] && pt_high<ptmax[jpt] && cent>=cent_min[jpt] && cent<cent_max[jpt]){
                        
                        eff_pt_high = p_eff_pt[jpt]->GetBinContent(p_eff_pt[jpt]->FindBin(pt_high));
                        
                    }
                    
                }
                
            }
            
            a0 = eff_pt_low / ( (pt_low - pt_mid) * (pt_low - pt_high) );
            a1 = eff_pt_mid / ( (pt_mid - pt_low ) * (pt_mid - pt_high) );
            a2 = eff_pt_high / ( (pt_high - pt_low) * (pt_high - pt_mid) );
            
            eff_pt = a0 * (pt - pt_mid) * (pt - pt_high) + a1 * (pt - pt_low) * (pt - pt_high) + a2 * (pt - pt_low) * (pt - pt_mid);
            
            if(eff_pt_low == 0 || eff_pt_mid == 0 || eff_pt_high == 0) eff_pt = 0;
            
            // 5-point interpolate eff_accept
            if(phi > p_eff_accept[ipt]->GetXaxis()->GetBinCenter(p_eff_accept[ipt]->GetXaxis()->FindBin(phi))){
                
                phi_low = p_eff_accept[ipt]->GetXaxis()->GetBinCenter(p_eff_accept[ipt]->GetXaxis()->FindBin(phi));
                phi_high = p_eff_accept[ipt]->GetXaxis()->GetBinCenter(p_eff_accept[ipt]->GetXaxis()->FindBin(phi) + 1);
                phi_bin_low = p_eff_accept[ipt]->GetXaxis()->FindBin(phi);
                phi_bin_high = p_eff_accept[ipt]->GetXaxis()->FindBin(phi) + 1;
            }
            else{
                
                phi_low = p_eff_accept[ipt]->GetXaxis()->GetBinCenter(p_eff_accept[ipt]->GetXaxis()->FindBin(phi) - 1);
                phi_high = p_eff_accept[ipt]->GetXaxis()->GetBinCenter(p_eff_accept[ipt]->GetXaxis()->FindBin(phi));
                phi_bin_low = p_eff_accept[ipt]->GetXaxis()->FindBin(phi) - 1;
                phi_bin_high = p_eff_accept[ipt]->GetXaxis()->FindBin(phi);
                
            }
            
            if( phi_high > TMath::Pi() ){
                
                phi_bin_high = 1;
                
            }
            
            if( phi_low < -TMath::Pi() ){
                
                phi_bin_low = p_eff_accept[ipt]->GetNbinsX();
                
            }
            
            if(eta > p_eff_accept[ipt]->GetYaxis()->GetBinCenter(p_eff_accept[ipt]->GetYaxis()->FindBin(eta))){
                
                eta_low = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
                eta_high = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(p_eff_accept[ipt]->GetYaxis()->FindBin(eta) + 1);
                eta_bin_low = p_eff_accept[ipt]->GetYaxis()->FindBin(eta);
                eta_bin_high = p_eff_accept[ipt]->GetYaxis()->FindBin(eta) + 1;
            }
            else{
                
                eta_low = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(p_eff_accept[ipt]->GetYaxis()->FindBin(eta) - 1);
                eta_high = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
                eta_bin_low = p_eff_accept[ipt]->GetYaxis()->FindBin(eta) - 1;
                eta_bin_high = p_eff_accept[ipt]->GetYaxis()->FindBin(eta);
                
            }
            
            if(eta_low < -2.4){
                
                eta_low = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(1);
                eta_high = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(2);
                eta_bin_low = 1;
                eta_bin_high = 2;
                
            }
            
            if(eta_high > 2.4){
                
                eta_low = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(p_eff_accept[ipt]->GetNbinsY() - 1 );
                eta_high = p_eff_accept[ipt]->GetYaxis()->GetBinCenter(p_eff_accept[ipt]->GetNbinsY());
                eta_bin_low = p_eff_accept[ipt]->GetNbinsY() - 1;
                eta_bin_high = p_eff_accept[ipt]->GetNbinsY();
                
            }
            
            Q11 = p_eff_accept[ipt]->GetBinContent(phi_bin_low,eta_bin_low);
            Q12 = p_eff_accept[ipt]->GetBinContent(phi_bin_low,eta_bin_high);
            Q21 = p_eff_accept[ipt]->GetBinContent(phi_bin_high,eta_bin_low);
            Q22 = p_eff_accept[ipt]->GetBinContent(phi_bin_high,eta_bin_high);
            
            eff_accept = ( 1/ ( (phi_high - phi_low) * (eta_high - eta_low) ) ) * ( Q11 * (phi_high - phi) * (eta_high - eta) + Q21 * (phi - phi_low) * (eta_high - eta) + Q12 * (phi_high - phi) * (eta - eta_low) + Q22 * (phi - phi_low) * (eta - eta_low) );
            
            if(Q11 == 0 || Q12 == 0 || Q21 == 0 || Q22 == 0 ) eff_accept = 0;
            
            eff = eff_pt * eff_accept;
            
            if( eff==0 ){
                if(pt > 100) eff = 0.8;
                else eff = 1;
            }
            
        }
        
    }
    
    return eff;
    
}

//-----------------------------------------------------------------------------





    