#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>

//unsigned int selectedRun = 182972;
double vtxCut= 15.;
TString fileName = "/rfs/jcastle/PbPb2015/EbyETree_InitialTest_TkEffBlackBox_Runs262620-263035.root";

const int NCENT_PBPB = 40;
double cent_min[NCENT_PBPB] = {0.0 ,2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5};
double cent_max[NCENT_PBPB] = {2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0};
static const double centBinsDefault[] = {0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0,42.5, 45.0, 47.5, 50.0, 52.5, 55.0, 57.5, 60.0, 62.5, 65.0, 67.5, 70.0, 72.5, 75.0, 77.5, 80.0, 82.5, 85.0, 87.5, 90.0, 92.5, 95.0, 97.5, 100.0};

static const int nptbinsDefault = 16;
static const double ptbinsDefault[]={0.2,  0.3,  0.4,  0.5,  0.6,  0.8,  1.0,  1.2,  1.6,  2.0, 2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  8.0};
static const int netabinsDefault = 12;
static const double etabinsDefault[]={-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0,  0.4,  0.8, 1.2,  1.6,  2.0,  2.4};

static const int ptBinMin = 7;  //-- for(int ipt = ptBinMin; ipt <= ptBinMax; ipt++) --> GetBinContent(ipt+1,ieta+1);
static const int ptBinMax = 16;
static const int etaBinMin = 0; //-- for(int ieta = etaBinMin; ieta <= etaBinMax; ieta++) --> GetBinContent(ipt+1,ieta+1);
static const int etaBinMax = 12;

TFile * tf;

TFile * fRes;
TDirectory * SubEvt_0;
TDirectory * SubEvt_1;
TDirectory * FullEvt;

TTree * tree;

unsigned int runno_;
double centval;
double vtx;
TH2D * sumw;
TH2D * sumwqx;
TH2D * sumwqy;

double V2RawX_TkEff_0[NCENT_PBPB+1];
double V2RawY_TkEff_0[NCENT_PBPB+1];
double V2RawX_TkEff_1[NCENT_PBPB+1];
double V2RawY_TkEff_1[NCENT_PBPB+1];
double V2RawX_TkEff_full[NCENT_PBPB+1];
double V2RawY_TkEff_full[NCENT_PBPB+1];

double sumw_0[NCENT_PBPB+1];
double sumw_1[NCENT_PBPB+1];
double sumw_full[NCENT_PBPB+1];

double V2DetX_TkEff_0[NCENT_PBPB+1];
double V2DetY_TkEff_0[NCENT_PBPB+1];
double V2DetX_TkEff_1[NCENT_PBPB+1];
double V2DetY_TkEff_1[NCENT_PBPB+1];
double V2DetX_TkEff_full[NCENT_PBPB+1];
double V2DetY_TkEff_full[NCENT_PBPB+1];

int Nevents[NCENT_PBPB+1];
int NFails[NCENT_PBPB+1];

//
// MAIN
//
void makeV2Det(){

    tf = new TFile(fileName);
    tree = (TTree *) tf->Get("ebyeana/tree");

    sumwqx = new TH2D("sumwqx","sumwqx",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    sumwqy = new TH2D("sumwqy","sumwqy",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    sumw =  new TH2D("sumw","sumw",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    
    tree->SetBranchAddress("Run",         &runno_);
    tree->SetBranchAddress("Cent",        &centval);
    tree->SetBranchAddress("Vtx",         &vtx);
    tree->SetBranchAddress("sumw",        &sumw);
    tree->SetBranchAddress("sumwqx",      &sumwqx);
    tree->SetBranchAddress("sumwqy",      &sumwqy);
    
    fRes = new TFile("V2Det.root","recreate");
    SubEvt_0 = fRes->mkdir("SubEvt_0");
    SubEvt_1 = fRes->mkdir("SubEvt_1");
    FullEvt = fRes->mkdir("FullEvt");
    
    SubEvt_0->cd();
    TH1D * hV2DetX_TkEff_0 = new TH1D("hV2DetX_TkEff_0","hV2DetX_TkEff_0",NCENT_PBPB,centBinsDefault);
    TH1D * hV2DetY_TkEff_0 = new TH1D("hV2DetY_TkEff_0","hV2DetY_TkEff_0",NCENT_PBPB,centBinsDefault);
    
    SubEvt_1->cd();
    TH1D * hV2DetX_TkEff_1 = new TH1D("hV2DetX_TkEff_1","hV2DetX_TkEff_1",NCENT_PBPB,centBinsDefault);
    TH1D * hV2DetY_TkEff_1 = new TH1D("hV2DetY_TkEff_1","hV2DetY_TkEff_1",NCENT_PBPB,centBinsDefault);
    
    FullEvt->cd();
    TH1D * hV2DetX_TkEff_full = new TH1D("hV2DetX_TkEff_full","hV2DetX_TkEff_full",NCENT_PBPB,centBinsDefault);
    TH1D * hV2DetY_TkEff_full = new TH1D("hV2DetY_TkEff_full","hV2DetY_TkEff_full",NCENT_PBPB,centBinsDefault);
    
    
    //-- initialize all variables
    for(int icent = 0; icent<(NCENT_PBPB+1); icent++){
        
	V2DetX_TkEff_0[icent] = 0.;
        V2DetY_TkEff_0[icent] = 0.;
        V2DetX_TkEff_1[icent] = 0.;
        V2DetY_TkEff_1[icent] = 0.;
        V2DetX_TkEff_full[icent] = 0.;
        V2DetY_TkEff_full[icent] = 0.;

        Nevents[icent] = 0;
        NFails[icent] = 0;
        
    }
    
    
    //
    // Calculate Vn_det
    //
    
    cout<<"Begin DETECTOR loop, contains "<<tree->GetEntries()<<" Events"<<endl;
    
    //for(int ievent = 0; ievent<10000; ievent++) {
    for(int ievent = 0; ievent<tree->GetEntries(); ievent++) {
        
      if((ievent+1)% 500000 == 0) cout<<"Processing Event "<<ievent+1<<"\t"<<(100.*(ievent+1)/tree->GetEntries())<<"% Completed"<<endl;
        
      tree->GetEntry(ievent);
      
      //-- Vertex Cut
      if(TMath::Abs(vtx) > vtxCut) continue;
      
      //-- Calculate centbin
      int icent = (centval - cent_min[0]) / 2.5;

      //-- Reset Raw and sumw values
      V2RawX_TkEff_0[icent]    = 0.;
      V2RawY_TkEff_0[icent]    = 0.;
      V2RawX_TkEff_1[icent]    = 0.;
      V2RawY_TkEff_1[icent]    = 0.;
      V2RawX_TkEff_full[icent] = 0.;
      V2RawY_TkEff_full[icent] = 0.;            
      
      sumw_0[icent]    = 0;
      sumw_1[icent]    = 0;
      sumw_full[icent] = 0;
      
      Nevents[icent]++;
      Nevents[40]++;
      
      for(int ipt = ptBinMin; ipt <= ptBinMax; ipt++){
	for(int ieta = etaBinMin; ieta <= etaBinMax; ieta++){
	  
	  if(sumw->GetBinContent(ipt+1,ieta+1) !=0){
	    
	    //-- Subevent 0 (eta >= 0)
	    if(etabinsDefault[ieta] >= 0){
	      V2RawX_TkEff_0[icent]     += sumwqx->GetBinContent(ipt+1,ieta+1);
	      V2RawY_TkEff_0[icent]     += sumwqy->GetBinContent(ipt+1,ieta+1);
	      V2RawX_TkEff_0[40]        += sumwqx->GetBinContent(ipt+1,ieta+1);
	      V2RawY_TkEff_0[40]        += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	      sumw_0[icent]             += sumw->GetBinContent(ipt+1,ieta+1);
	      sumw_0[40]                += sumw->GetBinContent(ipt+1,ieta+1);
	    }
	    //-- Subevent 1 (eta < 0)
	    else{
	      V2RawX_TkEff_1[icent]     += sumwqx->GetBinContent(ipt+1,ieta+1);
	      V2RawY_TkEff_1[icent]     += sumwqy->GetBinContent(ipt+1,ieta+1);
	      V2RawX_TkEff_1[40]        += sumwqx->GetBinContent(ipt+1,ieta+1);
	      V2RawY_TkEff_1[40]        += sumwqy->GetBinContent(ipt+1,ieta+1);
	      
	      sumw_1[icent]             += sumw->GetBinContent(ipt+1,ieta+1);
	      sumw_1[40]                += sumw->GetBinContent(ipt+1,ieta+1);
	    }
	    //-- Full Event
	    V2RawX_TkEff_full[icent]    += sumwqx->GetBinContent(ipt+1,ieta+1);
	    V2RawY_TkEff_full[icent]    += sumwqy->GetBinContent(ipt+1,ieta+1);
	    V2RawX_TkEff_full[40]       += sumwqx->GetBinContent(ipt+1,ieta+1);
	    V2RawY_TkEff_full[40]       += sumwqy->GetBinContent(ipt+1,ieta+1);
	    
	    sumw_full[icent]            += sumw->GetBinContent(ipt+1,ieta+1);
	    sumw_full[40]               += sumw->GetBinContent(ipt+1,ieta+1);
	  }
	  
	}
	
      }
      
      //-- Only use events that have tracks in all subevents
      if(sumw_0[icent] == 0 || sumw_1[icent] == 0 || sumw_full[icent] == 0){
	NFails[icent]++;
	NFails[40]++;
      }
      else{      
	V2DetX_TkEff_0[icent]    += V2RawX_TkEff_0[icent] / sumw_0[icent];
	V2DetY_TkEff_0[icent]    += V2RawY_TkEff_0[icent] / sumw_0[icent];
	V2DetX_TkEff_1[icent]    += V2RawX_TkEff_1[icent] / sumw_1[icent];
	V2DetY_TkEff_1[icent]    += V2RawY_TkEff_1[icent] / sumw_1[icent];
	V2DetX_TkEff_full[icent] += V2RawX_TkEff_full[icent] / sumw_full[icent];
	V2DetY_TkEff_full[icent] += V2RawY_TkEff_full[icent] / sumw_full[icent];
	
	V2DetX_TkEff_0[40]       += V2RawX_TkEff_0[icent] / sumw_0[icent];
	V2DetY_TkEff_0[40]       += V2RawY_TkEff_0[icent] / sumw_0[icent];
	V2DetX_TkEff_1[40]       += V2RawX_TkEff_1[icent] / sumw_1[icent];
	V2DetY_TkEff_1[40]       += V2RawY_TkEff_1[icent] / sumw_1[icent];
	V2DetX_TkEff_full[40]    += V2RawX_TkEff_full[icent] / sumw_1[icent];
	V2DetY_TkEff_full[40]    += V2RawY_TkEff_full[icent] / sumw_1[icent];
      }

    }
    
    //-- Average V2Det over all events for each centrality bin
    for(int icent = 0; icent < (NCENT_PBPB+1); icent++){
    
        V2DetX_TkEff_0[icent]    /= ( (double) Nevents[icent] - (double) NFails[icent] );
        V2DetY_TkEff_0[icent]    /= ( (double) Nevents[icent] - (double) NFails[icent] );
        V2DetX_TkEff_1[icent]    /= ( (double) Nevents[icent] - (double) NFails[icent] );
        V2DetY_TkEff_1[icent]    /= ( (double) Nevents[icent] - (double) NFails[icent] );
        V2DetX_TkEff_full[icent] /= ( (double) Nevents[icent] - (double) NFails[icent] );
        V2DetY_TkEff_full[icent] /= ( (double) Nevents[icent] - (double) NFails[icent] );
        
        if(icent==40) continue;
        
	//-- Populate histograms that will be used by ReadTree_normDet.C
        hV2DetX_TkEff_0->SetBinContent(icent+1,V2DetX_TkEff_0[icent]);
        hV2DetX_TkEff_1->SetBinContent(icent+1,V2DetX_TkEff_1[icent]);
        hV2DetX_TkEff_full->SetBinContent(icent+1,V2DetX_TkEff_full[icent]);
        
        hV2DetY_TkEff_0->SetBinContent(icent+1,V2DetY_TkEff_0[icent]);
        hV2DetY_TkEff_1->SetBinContent(icent+1,V2DetY_TkEff_1[icent]);
        hV2DetY_TkEff_full->SetBinContent(icent+1,V2DetY_TkEff_full[icent]);
        
    }
    //-- Average over all events and centralities (to be spat out by the subsequent lines of code as a quick check 
    cout<<"double V2DetX_TkEff_0 = "<<V2DetX_TkEff_0[40]<<";"<<"\t//(over all events and centralities)"<<endl;
    cout<<"double V2DetY_TkEff_0 = "<<V2DetY_TkEff_0[40]<<";"<<"\t//(over all events and centralities)"<<endl;
    
    cout<<"double V2DetX_TkEff_1 = "<<V2DetX_TkEff_1[40]<<";"<<"\t//(over all events and centralities)"<<endl;
    cout<<"double V2DetY_TkEff_1 = "<<V2DetY_TkEff_1[40]<<";"<<"\t//(over all events and centralities)"<<endl;
    
    cout<<"double V2DetX_TkEff_full = "<<V2DetX_TkEff_full[40]<<";"<<"\t//(over all events and centralities)"<<endl;
    cout<<"double V2DetY_TkEff_full = "<<V2DetY_TkEff_full[40]<<";"<<"\t//(over all events and centralities)"<<endl;
    
    
    fRes->Write();
    cout<<"File written, process completed"<<endl;

}









    
    
