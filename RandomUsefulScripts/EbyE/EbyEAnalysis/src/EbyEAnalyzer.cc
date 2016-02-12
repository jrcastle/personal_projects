// -*- C++ -*-
//
// Package:    EbyEAnalyzer
// Class:      EbyEAnalyzer
//
/**\class EbyEAnalyzer EbyEAnalyzer.cc RecoHI/EbyEAnalyzer/src/EbyEAnalyzer.cc
 
 Description: <one line class summary>
 
 Implementation:
 <Notes on implementation>
 */
//
// Original Author:  Sergey Petrushanko
//         Created:  Fri Jul 11 10:05:00 2008
// $Id: EbyEAnalyzer.cc,v 1.18 2011/10/07 09:41:29 yilmaz Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
//#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"
//#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
//#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
//#include "HeavyIonsAnalysis/EbyEAnalysis/interface/TrackEfficiency.h"
//#include "HeavyIonsAnalysis/EbyEAnalysis/interface/TrackEfficiencyMC.h"
#include "/afs/cern.ch/work/j/jcastle/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/TrackingCode/HIRun2015Ana/macros/TrackCorrector3D.h"

#include "TROOT.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>

#include <vector>
#include <iostream>
using std::vector;
using std::rand;
using namespace std;
//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
//#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"

/*
static const int nptbinsDefault = 16;
static const double ptbinsDefault[]={
    0.2,  0.3,  0.4,  0.5,  0.6,  0.8,  1.0,  1.2,  1.6,  2.0,
    2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  8.0};
static const int netabinsDefault = 12;
static const double etabinsDefault[]={-2.4, -2.0, -1.6, -1.2,
				      -0.8, -0.4, 0.0,  0.4,  0.8,
				      1.2,  1.6,  2.0,  2.4};
*/
static const int nptbinsDefault = 11;                                                                                                                                                                                            
static const double ptbinsDefault[]={1.00, 1.25, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00, 6.00, 7.00, 8.00};
static const int netabinsDefault = 14;
static const double etabinsDefault[]= {-2.4, -2.0, -1.6, -1.2, -1.0, -0.8, -0.4, 0.0, 0.4, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4};

static const int NCENT = 10;
double centBinWidth = 10.;
int cmin[NCENT] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
int cmax[NCENT] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
static const int Nvtx = 10;
double vtxBinWidth = 5.;
int vmin[Nvtx] = {-25, -20, -15, -10, -5, 0, 5, 10, 15, 20};
int vmax[Nvtx] = {-20, -15, -10, -5, 0, 5, 10, 15, 20, 25};


//
// class decleration
//

class EbyEAnalyzer : public edm::EDAnalyzer {
public:
    explicit EbyEAnalyzer(const edm::ParameterSet&);
    ~EbyEAnalyzer();
    
private:
    //edm::Service<TFileService> fs;
    
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    
    // ----------member data ---------------------------
    int nOrder_;

    unsigned int runno_;
    
    edm::Service<TFileService> fs;
    
    edm::InputTag vertexTag_;
    edm::Handle<std::vector<reco::Vertex>> vertex_;
    
    edm::InputTag trackTag_;
    edm::Handle<reco::TrackCollection> trackCollection_;
    
    edm::EDGetTokenT<reco::Centrality> CentralityTag_;
    edm::EDGetTokenT<int> CentralityBinTag_;

    double caloCentRef_;
    double caloCentRefWidth_;
    int caloCentRefMinBin_;
    int caloCentRefMaxBin_;
    int CentBinCompression_;
    double centval;
    
    int vs_sell;   // vertex collection size
    float vzr_sell;
    double vtx;
    
    bool loadDB_;
    double minpt_;
    double maxpt_;
    double etaMax_;

    bool useTeff_;
    TrackCorrector3D * teff;
    
    TTree * tree;
    
    TRandom3 * ran;
    
    bool FirstEvent_;
    string effTable_;
    
    bool Branch_Cent;
    bool Branch_Vtx;
    bool Branch_sumw;
    bool Branch_sumwqx;
    bool Branch_sumwqy;
    bool Branch_Run;
    
    bool Subevent_Standard;
    bool countHisto;
    
    TH2D * sumw;
    TH2D * sumwqx;
    TH2D * sumwqy;

    
    TH2D * count[NCENT][Nvtx];

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EbyEAnalyzer::EbyEAnalyzer(const edm::ParameterSet& iConfig):
  runno_(0),
  CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralityTag_"))),
  CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinTag_")))
{
    runno_ = 0;
    loadDB_ = kTRUE;
    FirstEvent_ = kTRUE;
    
    ran = new TRandom3(0);
    
    vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
    trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
    nOrder_ = iConfig.getUntrackedParameter<int>("EPOrder_",2);
    caloCentRef_ = iConfig.getUntrackedParameter<double>("caloCentRef_",80.);
    caloCentRefWidth_ = iConfig.getUntrackedParameter<double>("caloCentRefWidth_",5.);
    CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
    useTeff_ = iConfig.getUntrackedParameter<bool>("useTeff_");
    effTable_ = iConfig.getParameter<std::string>("effTable_");
    teff = 0;

    if(useTeff_){
        teff = new TrackCorrector3D(effTable_);
	teff->load("HITrackCorrections");
    }
    
    loadDB_ = iConfig.getUntrackedParameter<bool>("loadDB_",true);
    
    Branch_Cent = iConfig.getUntrackedParameter<bool>("Branch_Cent",true);
    Branch_Vtx = iConfig.getUntrackedParameter<bool>("Branch_Vtx",true);
    Branch_sumw = iConfig.getUntrackedParameter<bool>("Branch_sumw",true);
    Branch_sumwqx = iConfig.getUntrackedParameter<bool>("Branch_sumwqx",true);
    Branch_sumwqy = iConfig.getUntrackedParameter<bool>("Branch_sumwqy",true);
    Branch_Run = iConfig.getUntrackedParameter<bool>("Branch_Run",true);
    Subevent_Standard = iConfig.getUntrackedParameter<bool>("Subevent_Standard",true);
    countHisto = iConfig.getUntrackedParameter<bool>("countHisto",false);
    
    if(Subevent_Standard) std::cout<<"Standard subevent selection will be used"<<std::endl;
    if(countHisto) std::cout<<"Setting up counting histograms..."<<std::endl;

    tree = fs->make<TTree>("tree","EP tree");
    
    if(countHisto){
        for(int ivtx = 0; ivtx < Nvtx; ivtx++){
            for( int icent = 0; icent < NCENT; icent++){
                count[icent][ivtx] = fs->make<TH2D>(Form("count_c%i_%i_v%i_%i",cmin[icent],cmax[icent],vmin[ivtx],vmax[ivtx]),Form("count_c%i_%i_v%i_%i",cmin[icent],cmax[icent],vmin[ivtx],vmax[ivtx]),50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
                count[icent][ivtx]->GetXaxis()->SetTitle("#phi [rad]");
                count[icent][ivtx]->GetYaxis()->SetTitle("#eta");
                count[icent][ivtx]->SetOption("colz");
            }
        }
    }else{
        for(int ivtx = 0; ivtx < Nvtx; ivtx++){
            for( int icent = 0; icent < NCENT; icent++){
                count[icent][ivtx]   = 0;
            }
        }
    }
    
    sumwqx = fs->make<TH2D>("sumwqx","sumwqx",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    sumwqy = fs->make<TH2D>("sumwqy","sumwqy",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    sumw =  fs->make<TH2D>("sumw","sumw",nptbinsDefault,ptbinsDefault, netabinsDefault, etabinsDefault);
    
    sumwqx->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    sumwqx->GetYaxis()->SetTitle("#eta");
    sumwqx->SetOption("colz");
    
    sumwqy->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    sumwqy->GetYaxis()->SetTitle("#eta");
    sumwqy->SetOption("colz");
    
    sumw->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    sumw->GetYaxis()->SetTitle("#eta");
    sumw->SetOption("colz");
    
    if(Branch_Cent)              tree->Branch("Cent",        &centval,     "cent/D");
    if(Branch_Vtx)               tree->Branch("Vtx",         &vtx,         "vtx/D");
    if(Branch_sumwqx)            tree->Branch("sumwqx", "TH2D",  &sumwqx, 128000, 0);
    if(Branch_sumwqy)            tree->Branch("sumwqy", "TH2D",  &sumwqy, 128000, 0);
    if(Branch_sumw)              tree->Branch("sumw",   "TH2D",  &sumw,   128000, 0);
    if(Branch_Run)               tree->Branch("Run",         &runno_,      "run/i");
    
    minpt_ = iConfig.getUntrackedParameter<double>("minpt_",0.0);
    maxpt_ = iConfig.getUntrackedParameter<double>("maxpt_",8.0);
    etaMax_ = iConfig.getUntrackedParameter<double>("etaMax_",2.4);
    FirstEvent_ = kTRUE;

    
}


EbyEAnalyzer::~EbyEAnalyzer()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called to produce and analyze the data  ------------
void
EbyEAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using namespace reco;
    
    Bool_t newrun = kFALSE;
    if(runno_ != iEvent.id().run()) newrun = kTRUE;
    runno_ = iEvent.id().run();
    if(FirstEvent_ || newrun) {
        FirstEvent_ = kFALSE;
        newrun = kFALSE;
        //
        //Get flattening parameter file.
        //
        //edm::ESHandle<RPFlatParams> flatparmsDB_;
        //iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
        //LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
        //if(!db->IsSuccess()) {
        //    cout<<"Failed to load DB!"<<endl;
        //    loadDB_ = kFALSE;
        //}
    }//First event
    
    //
    //Get Centrality
    //
    edm::Handle<int> cbin_;
    iEvent.getByToken(CentralityBinTag_,cbin_);
    int hiBin = *cbin_; //HF tower centrality
    int hiBinHF = hiBin;
    centval = 0.5 * (hiBinHF + 0.5);
    //
    //Get Vertex
    //
    iEvent.getByLabel(vertexTag_,vertex_);
    const reco::VertexCollection * vertices3 = vertex_.product();
    vs_sell = vertices3->size();
    if(vs_sell>0) {
        vzr_sell = vertices3->begin()->z();
    } else vzr_sell = -999.9;
        
    vtx = vzr_sell;
    
    const VertexCollection * recoVertices = vertex_.product();
    int primaryvtx = 0;
    math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
    double vxError = (*recoVertices)[primaryvtx].xError();
    double vyError = (*recoVertices)[primaryvtx].yError();
    double vzError = (*recoVertices)[primaryvtx].zError();
    
    if( TMath::Abs(vtx) < 25.){
    
        //
        //Tracking part
        //
        sumw->Reset();
        sumwqx->Reset();
        sumwqy->Reset();

        iEvent.getByLabel(trackTag_,trackCollection_);

        for(TrackCollection::const_iterator itTrack = trackCollection_->begin();
            itTrack != trackCollection_->end();
            ++itTrack) {
            
            if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
            if ( itTrack->charge() == 0 ) continue;
            double d0 = -1.* itTrack->dxy(v1);
            double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
            double dz=itTrack->dz(v1);
            double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
            if ( fabs(itTrack->eta()) > etaMax_ ) continue;
            if ( fabs( dz/dzerror ) > 3. ) continue;
            if ( fabs( d0/derror ) > 3. ) continue;
            if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;
	    //if ( itTrack->chi2()/itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) continue;
	    //if ( itTrack->algo() < 4 || itTrack->algo() > 7 ) continue;
	    //if ( itTrack->numberOfValidHits() < 10 ) continue;
            if(itTrack->pt() < minpt_) continue;
            if(itTrack->pt() > maxpt_) continue;
            
            if(countHisto){
                int vtxbin = (vtx - vmin[0]) / vtxBinWidth;
                int centbin = (centval - cmin[0]) / centBinWidth;
                count[centbin][vtxbin]->Fill(itTrack->phi(),itTrack->eta());
            }
            
	    double w = 1.;
	    if(useTeff_) w = teff->getWeight(itTrack->pt(), itTrack->eta(), hiBinHF);
	    if(w > 0.0) {
                    
	        //Do standard subevent selection
	        if(Subevent_Standard){
                        
		    sumw->Fill(itTrack->pt(), itTrack->eta(), 1./w);
		    sumwqx->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(nOrder_*itTrack->phi())/w);
		    sumwqy->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(nOrder_*itTrack->phi())/w);
                        
		}
                    
	    }
                
        }
        tree->Fill();
    }
}

// ------------ method called once each job just before starting event loop  ------------
void
EbyEAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
EbyEAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(EbyEAnalyzer);
