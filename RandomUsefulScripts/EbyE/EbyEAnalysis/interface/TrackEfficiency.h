#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TDirectory.h"

static const int NIOV = 10;
static const int npt = 29;
double ptmin[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.55, 0.55, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.65, 0.65, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 8.0};
double ptmax[] = {0.55, 0.55, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.65, 0.65, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0, 3.0, 8.0, 8.0, 8.0, 300};
double cent_min[] = {0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 60, 100, 0, 20, 40, 0};
double cent_max[] = {20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 60, 100, 200, 20, 40, 200, 200};
static const int NCENT = 10;
double centBinWidth = 10.;
int cmin[NCENT] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
int cmax[NCENT] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
static const int Nvtx = 10;
double vtxBinWidth = 5.;
int vmin[Nvtx] = {-25, -20, -15, -10, -5, 0, 5, 10, 15, 20};
int vmax[Nvtx] = {-20, -15, -10, -5, 0, 5, 10, 15, 20, 25};
unsigned int runs[] = {181611, 181671, 181701, 181801, 182001, 182201, 182401, 182601, 182801, 182951, 183101};


class TrackEfficiency{
  
public:
  TrackEfficiency(TString effOption);
  double getEfficiencies(double pt, double cent, double phi, double eta);
  double getPtEfficiency(double pt, double cent, double phi, double eta);
  double getEffAcceptVtxBins(double cent, double vtx, double phi, double eta);
  int setRunNumber(unsigned int runno_);
  
private:
  TFile * f_eff;
  TFile * f_effPt;
  TDirectory * effDir[npt];
  TH2D * eff[10];
  TH2D * effVtx[NCENT][Nvtx];
  int IOVDir;
  bool DATA;
  TString eOption;
    
  TProfile * p_eff_pt[npt];
  double eff_pt_low;
  double eff_pt_mid;
  double eff_pt_high;
  double pt_low;
  double pt_mid;
  double pt_high;
    
  double a0;
  double a1;
  double a2;
    
  double eff_pt;
  
};

//-----------------------------------------------------------------------------
//                                  METHODS
//-----------------------------------------------------------------------------

TrackEfficiency::TrackEfficiency(TString effOption){
        
    eOption = effOption;
    f_eff = new TFile(eOption);
	setRunNumber(runs[0]);    
 
    f_effPt = new TFile("akVs3Calo_20140920.root");
    
    for(int ipt=0; ipt<npt;ipt++){
        
        effDir[ipt]= (TDirectory*) f_effPt->Get(Form("eff_pt%d_%d_cent%d_%d",(int)(100*ptmin[ipt]),(int)(100*ptmax[ipt]),(int)(cent_min[ipt]/2),(int)(cent_max[ipt]/2)));
        effDir[ipt]->cd();
        p_eff_pt[ipt]=(TProfile*) effDir[ipt]->Get("p_eff_pt");
        
    }
    
    
}
//-----------------------------------------------------------------------------

int TrackEfficiency::setRunNumber(unsigned int runno_){
    
  IOVDir = -1;
  for(int iIOV = 0; iIOV < NIOV; iIOV++){
    if(runno_ >= runs[iIOV] && runno_ < runs[iIOV+1]) {
      IOVDir = iIOV;
      break;
    }
  }
  if(IOVDir<0) return 0;
  
  if(eOption == "Eff.root"){
    for(int i = 0; i<10; i++) {
      eff[i] = (TH2D *) f_eff->Get(Form("%d/Eff_%d_%d",runs[IOVDir],10*i,10*i+10));
    }
  }
    
  if(eOption == "EffVtx.root"){
    for(int icent = 0; icent < 10; icent++){
      for(int ivtx = 0; ivtx < Nvtx; ivtx++){
        effVtx[icent][ivtx] = (TH2D*) f_eff->Get(Form("%d/vtx_%i_%i/eff_c%i_%i",runs[IOVDir],vmin[ivtx],vmax[ivtx],cmin[icent],cmax[icent]));
      }
    }
  }
  return 1;
}

//-----------------------------------------------------------------------------

double TrackEfficiency::getEfficiencies(double pt, double cent, double phi, double eta){
    
  int icent = cent/10;
  if(icent>=0 && icent<10) {
    return eff[icent]->GetBinContent( eff[icent]->GetXaxis()->FindBin(phi),eff[icent]->GetYaxis()->FindBin(eta));
    
  }
  
  return 1.;
}

//-----------------------------------------------------------------------------


double TrackEfficiency::getEffAcceptVtxBins(double cent, double vtx, double phi, double eta){
    
    int icent = cent/10;
    int ivtx = (vtx - vmin[0]) / vtxBinWidth;
    if(icent>=0 && icent<10 && ivtx>=0 && ivtx<10) {
        return effVtx[icent][ivtx]->GetBinContent( effVtx[icent][ivtx]->GetXaxis()->FindBin(phi),effVtx[icent][ivtx]->GetYaxis()->FindBin(eta));
        
    }
    
    return 1.;
}

//-----------------------------------------------------------------------------


double TrackEfficiency::getPtEfficiency(double pt, double cent, double phi, double eta){
    
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
            
            if(eff_pt_low == 0 || eff_pt_mid == 0 || eff_pt_high == 0) return 0.;
            else return eff_pt;
        
        }
        
    }
    
    return 1.;
    
}



    
