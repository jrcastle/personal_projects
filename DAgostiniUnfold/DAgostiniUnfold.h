//-- Unfolding procedure follows that proposed by D'Agostini in: 
//--       "A multidimensional unfolding method based on Bayes' theorem", D'Agostini, NIM-A 362, 487
//-- Error propagation based on improved approach to D'Agostini as outlined in:
//--       http://hepunx.rl.ac.uk/~adye/software/unfold/bayes_errors.pdf

#include <iostream>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"

class DAgostiniUnfold{

 public:
  //-- Constructor
  DAgostiniUnfold(TH2D * hresp);
  ~DAgostiniUnfold();

  //-- Methods
  int computeObsCovMatrix(int iter);
  int computeObsErrMatrix(int iter);
  int computeRespCovMatrix(int iter);
  int computeRespErrMatrix(int iter);
  double kronDelta(int i, int j);
  void makeEfficiencyVector();
  int Unfold(TH1D * hObserved, int niter, TH1D * hInitialPrior = 0);

  //-- Getters
  TH1D * GetHReco();
  TH2D * GetHRecoCov();
  TH1D * GetObserved();
  TH2D * GetResponse();

  //-- Setters
  void SetPrior(TH1D * hPrior, int iter);

 private:
  TH2D * respErrMat_;          //-- The error matrix for each iteration that arises from the uncertainties on the respose matrix elements
  TH2D * hresp_;               //-- Response matrix p( ei | cj ) with dimensions truth VS observed
  TH2D * Mij_;                 //-- The unfolding matrix with dimensions truth vs observed
  TH2D * obsErrMat_;           //-- The error matrix for each iteration that arises from the statistical uncertainties in the observed distributions
  TH2D * obsCovMat_;           //-- The covariance matrix for each iteration that arises from the statistical uncertainties in the observed distributions 
  TH2D * respCovMat_;          //-- The covariance matrix for each iteration that arises from the statistical uncertainties in the response matrix elements
  TH2D * recoCovMat_;          //-- The covariance matrix for the final, unfolded distribution.  Equal to the sum of obsCovMat_ + respCovMat_
  TH1D * hObs_;                //-- Observed, smeared histogram
  TH1D * hPrior_;              //-- "Prior" histogram, before iterations it will be the guess at the true, underlying histogram
  TH1D * hReco_;               //-- Unfolded histogram

  int Nobs;                    //-- Number of bins in the observed histogram
  int Ntru;                    //-- Number of bins in the true histogram
  std::vector<double> effVec_; //-- Detection efficiency vector = Sum[ P( ei | cj ) , { i, 1, Nobs} ]

};


//-- =======================================
//--               CONSTRUCTOR             
//-- =======================================

DAgostiniUnfold::DAgostiniUnfold(TH2D * hresp){

  std::cout<<"Constructing a new DAgostiniUnfold object..."<<std::endl;

  //-- Store the response matrix and observed distribution into class members
  hresp_  = (TH2D*) hresp->Clone("hresp_");

  Nobs = hresp_->GetNbinsX();
  Ntru = hresp_->GetNbinsY();

  //-- Normalize the response matrix by the projection onto the truth axis
  TH1D * hy = (TH1D*) hresp_->ProjectionY();
  SetPrior(hy, 0);
  for(int y = 1; y <= Ntru; y++){
    double ntrue = hy->GetBinContent(y);
    if(ntrue==0){
      for(int x = 1; x <= Nobs; x++){
	hresp_->SetBinContent(x, y, 0.);
      }
    }
    else{
      for(int x = 1; x <= Nobs; x++){
	double binCont   = hresp_->GetBinContent(x, y);
	double binErr    = hresp_->GetBinError(x, y);
	hresp_->SetBinContent(x, y, binCont / ntrue);
	hresp_->SetBinError(x, y, binErr / ntrue);
      }
    }
  }
  delete hy;  

  //-- Create a std::vector that contains the detection efficiencies for each truth bin
  makeEfficiencyVector();

  //-- Initialize the TH2D that corresponds to the error propagation matrix for the response matrix
  //-- In this case we'll collaps the i and k elements onto the x axis of the TH2D.  Details in the
  //-- computeObsErrMatrix method
  //--
  //-- dn(cj)/dP( ei | ck )
  //-- The respective axes correspond to:
  //-- x <==> ei ck
  //-- y <==> cj

  int nx = Nobs;
  int ny = Ntru;

  double xL = hresp_->ProjectionX()->GetBinLowEdge(1);
  double xH = hresp_->ProjectionX()->GetBinLowEdge(Nobs) + hresp_->ProjectionX()->GetBinWidth(Nobs);
  double yL = hresp_->ProjectionY()->GetBinLowEdge(1);
  double yH = hresp_->ProjectionY()->GetBinLowEdge(Ntru) + hresp_->ProjectionY()->GetBinWidth(Ntru);

  respErrMat_  = new TH2D("respErrMat_0", "respErrMat_0", nx*ny, xL, xH, ny, yL, yH);

  //-- Initialize the covariance matrices for the respective error propagation matrices, which will build
  //-- the final covariance martrix for the unfolded distribution Cov_unfold = Cov_obs + Cov_resp
  obsCovMat_  = new TH2D("obsCovMat_",  "obsCovMat_",  Ntru, yL, yH, Ntru, yL, yH);
  respCovMat_ = new TH2D("respCovMat_", "respCovMat_", Ntru, yL, yH, Ntru, yL, yH);
  recoCovMat_ = new TH2D("recoCovMat_", "recoCovMat_", Ntru, yL, yH, Ntru, yL, yH);

  //-- Initialize what will be the final unfolded distribution
  hReco_ = new TH1D("hReco_", "hReco_", Ntru, yL, yH);

  //-- Initialize what will be the unfolding matrix
  Mij_ = new TH2D("Mij_", "Mij_", Nobs, xL, xH, Ntru, yL, yH);

}

//-- =======================================
//--                DESTRUCTOR
//-- =======================================

DAgostiniUnfold::~DAgostiniUnfold(){

  std::cout<<"Destroying DAgostiniUnfold object..."<<std::endl;
  /*  
  delete respErrMat_;
  delete hresp_;
  delete Mij_;
  delete obsErrMat_;
  delete obsCovMat_;
  delete respCovMat_;
  delete recoCovMat_;
  delete hObs_;
  delete hPrior_;
  delete hReco_;
  */
}

//-- =======================================
//--            computeObsCovMatrix
//-- =======================================

int DAgostiniUnfold::computeObsCovMatrix(int iter){

  std::cout<<"Converting the observed error propagation matrix to a covariance matrix for iteration "<<iter<<std::endl;

  //-- Check name of obsErrMat_.  If it doesn't match the iteration, return 0
  string name = obsErrMat_->GetName();
  if(name != Form("obsErrMat_%i", iter)){
    std::cout<<"Using the wrong iteration of the observed error propagation matrix to construct its covariance matrix"<<std::endl;
    std::cout<<"We are on iteration "<<iter<<std::endl;
    std::cout<<"While the name of obsErrMat_ is "<< obsErrMat_->GetName() <<std::endl;

    return 0;

  }

  obsCovMat_->Reset();

  for(int k = 1; k <= Ntru; k++){
    for(int L = 1; L <= Ntru; L++){

      double doubSum = 0;
      for(int i = 1; i <= Nobs; i++){
	for(int j = 1; j <= Nobs; j++){

	  double errMat_ik = obsErrMat_->GetBinContent(i, k);
	  double covEiEj   = hObs_->GetBinContent(i) * kronDelta(i, j);
	  double errMat_jL = obsErrMat_->GetBinContent(j, L);

	  doubSum += errMat_ik * covEiEj * errMat_jL;

	}
      } //-- End double sum loop

      obsCovMat_->SetBinContent(k, L, doubSum);

    }
  }

  return 1;

}

//-- =======================================
//--            computeObsErrMatrix
//-- =======================================

int DAgostiniUnfold::computeObsErrMatrix(int iter){

  std::cout<<"Computing the observed error propagation matrix for iteration "<<iter<<std::endl;

  if(iter == 0){
    obsErrMat_ = (TH2D*) Mij_->Clone( Form("obsErrMat_%i", iter) );
    if(!obsErrMat_) return 0;
    return 1;
  }
  else{

    string name = obsErrMat_->GetName();
    if(name != Form("obsErrMat_%i", iter-1)){
      std::cout<<"obsErrMatPrev_ will not be set to the previous iteration"<<std::endl;
      std::cout<<"We are on iteration "<<iter<<std::endl;
      std::cout<<"The name of obsErrMat_ is "<< obsErrMat_->GetName() <<std::endl;
      std::cout<<"While it should be "<<Form("obsErrMat_%i", iter-1)<<std::endl;
      return 0;
    }

    TH2D obsErrMatPrev_ = *obsErrMat_;
    obsErrMat_->SetName( Form("obsErrMat_%i", iter) );


    for(int x = 1; x <= Nobs; x++){
      for(int y = 1; y <= Ntru; y++){

	double Mij         = Mij_->GetBinContent(x, y);
	double newEst      = hReco_->GetBinContent(y);
	double oldEst      = hPrior_->GetBinContent(y);
	double oldErrMatji = obsErrMatPrev_.GetBinContent(x, y);
	double sum         = 0.;

	for(int k = 1; k <= Nobs; k++){
	  for(int l = 1; l <= Ntru; l++){

	    double nEk         = hObs_->GetBinContent(k);
	    double n0Cl        = hPrior_->GetBinContent(l);
	    double epsl        = effVec_[l-1];
	    double Mki         = Mij_->GetBinContent(k, y);
	    double Mlk         = Mij_->GetBinContent(k, l);
	    double oldErrMatjl = obsErrMatPrev_.GetBinContent(x, l);

	    sum += ( nEk * epsl / n0Cl ) * Mki * Mlk * oldErrMatjl;

	  }
	} //-- End double sum loops

	double obsErrMatji = Mij + (newEst / oldEst) * oldErrMatji - sum;
	obsErrMat_->SetBinContent(x, y, obsErrMatji);

      }
    } //-- End matrix bins loops

    return 1;

  } // End else 

  return 0; //-- This return should never be reached, but is placed here to keep the compiler happy

}

//-- =======================================
//--            computeRespCovMatrix
//-- =======================================

int DAgostiniUnfold::computeRespCovMatrix(int iter){

  std::cout<<"Converting the response error propagation matrix to a covariance matrix for iteration "<<iter<<std::endl;

  //-- Check name of respErrMat_.  If it doesn't match the iteration, return 0
  string name = respErrMat_->GetName();
  if(name != Form("respErrMat_%i", iter)){
    std::cout<<"Using the wrong iteration of the observed error propagation matrix to construct its covariance matrix"<<std::endl;
    std::cout<<"We are on iteration "<<iter<<std::endl;
    std::cout<<"While the name of respErrMat_ is "<< respErrMat_->GetName() <<std::endl;

    return 0;

  }

  respCovMat_->Reset();

  //-- Cov_KL = A_Kj * B_jr * (A^T)_rL
  //-- A_Kj = response error propagation matrix
  //-- B_jr = diagonal matrix of the variances of each bin of the response matrix.
  //-- K = 1...Ntru
  //-- L = 1...Ntru
  //-- j = 1...Nobs*Ntru
  //-- r = 1...Nobs*Nrtu

  for(int k = 1; k <= Ntru; k++){
    for(int L = 1; L <= Ntru; L++){

      int    binjC = 1;
      int    binjE = 0;
      double sumj = 0;

      for(int j = 1; j <= Nobs*Ntru; j++){

	int    binrC = 1;
	int    binrE = 0;
	double sumr = 0;

	//-----
	for(int r = 1; r <= Nobs*Ntru; r++){
	
	  binrE++;
	  if( binrE > Nobs){
	    binrE = 1;
	    binrC++;
	  }

	  double respVar_xy = pow( hresp_->GetBinError(binrE, binrC), 2 );
	  double respErrMatT_rL = respErrMat_->GetBinContent(L, r);
	  sumr += respVar_xy * respErrMatT_rL; 

	}
	//----

	binjE++;
	if( binjE > Nobs){
	  binjE = 1;
	  binjC++;
	}

	double respVar_xy = pow( hresp_->GetBinError(binjE, binjC), 2 );
	double respErrMat_kj = respErrMat_->GetBinContent(k, j);
	sumj += respErrMat_kj * respVar_xy * sumr;

      }

      respCovMat_->SetBinContent(k, L , sumj);
      
    }
  }

  return 1;

}

//-- =======================================
//--           computeRespErrMatrix
//-- =======================================

int DAgostiniUnfold::computeRespErrMatrix(int iter){

  std::cout<<"Computing the response error propagation matrix for iteration "<<iter<<std::endl;

  //-- dn(cj)/dP( ei | ck )
  //-- Using a TH2D to represent this, so the respective axes correspond to:
  //-- x <==> ei ck 
  //-- y <==> cj

  if(iter == 0){

    respErrMat_->SetName( Form("respErrMat_%i", iter) );

    for(int i = 1; i <= Nobs; i++){
      for(int j = 1; j <= Ntru; j++){
	for(int k = 1; k <= Ntru; k++){

	  double epsj   = effVec_[j-1];
	  double n0Cj   = hPrior_->GetBinContent(j);
	  double nEi    = hObs_->GetBinContent(i);
	  double nHatCj = hReco_->GetBinContent(j);
	  double n0Ck   = hPrior_->GetBinContent(k);
	  double Mij    = Mij_->GetBinContent(i, j);

	  double sumfi = 0.;
          for(int L = 1; L <= Ntru; L++) sumfi += hresp_->GetBinContent(i, L) * hPrior_->GetBinContent(L);
          double fi     = sumfi;
          if(fi == 0) return 0;

	  double dncj_dPeick = (1/epsj) * ( n0Cj * nEi / fi - nHatCj ) * kronDelta(j,k) - ( n0Ck * nEi / fi ) * Mij;

	  //-- Collapse i and k to the "x" dimension of the matrix
	  //-- 
	  //--  y=1  __                                                         __
	  //--   .   |                                                           |
	  //--   .   |                                                           |
	  //--   .   |  k=1          k=2           k=3        ...         k=nC   |
	  //--   .   |                                                           |
	  //--   .   |                                                           |
	  //--  y=nC --                                                         --
	  //--        i=1...nE   i=nE+1...2nE   i=2nE+1...3nE ... i=(nC-1)nE+1...nCnE


	  int bin = i + ( k-1 ) * Nobs;
	  respErrMat_->SetBinContent(bin, j, dncj_dPeick);

	}
      }
    }

    return 1;

  } //-- End if(iter == 0)
  else{

    string name = respErrMat_->GetName();
    if(name != Form("respErrMat_%i", iter-1)){
      std::cout<<"respErrMatPrev_ will not be set to the previous iteration"<<std::endl;
      std::cout<<"We are on iteration "<<iter<<std::endl;
      std::cout<<"The name of respErrMat_ is "<< respErrMat_->GetName() <<std::endl;
      std::cout<<"While it should be "<<Form("respErrMat_%i", iter-1)<<std::endl;
      return 0;
    }

    TH2D respErrMatPrev_ = *respErrMat_;
    respErrMat_->SetName( Form("respErrMat_%i", iter) );

    for(int i = 1; i <= Nobs; i++){
      for(int j = 1; j <= Ntru; j++){
        for(int k = 1; k <= Ntru; k++){

          double epsj         = effVec_[j-1];
          double n0Cj         = hPrior_->GetBinContent(j);
          double nEi          = hObs_->GetBinContent(i);
          double nHatCj       = hReco_->GetBinContent(j);
          double n0Ck         = hPrior_->GetBinContent(k);
          double Mij          = Mij_->GetBinContent(i, j);

	  int bin = i + ( k-1 ) * Nobs;
	  double prevIter_ijk = respErrMatPrev_.GetBinContent(bin, j);

	  double sumfi = 0.;
	  for(int L = 1; L <= Ntru; L++) sumfi += hresp_->GetBinContent(i, L) * hPrior_->GetBinContent(L);
          double fi     = sumfi;
          if(fi == 0) return 0;

	  double doubSum = 0.;
	  for(int L = 1; L <= Nobs; L++){
	    for(int r = 1; r <= Ntru; r++){

	      double nEL          = hObs_->GetBinContent(L);
	      double MLj          = Mij_->GetBinContent(L, j);
	      double MLr          = Mij_->GetBinContent(L, r);

	      int bin = i + ( k-1 ) * Nobs;
	      double prevIter_irk = respErrMatPrev_.GetBinContent(bin, r);

	      doubSum += nEL * MLj * MLr * prevIter_irk;

	    }
	  } //-- End doubSum loop

          double dncj_dPeick = (1/epsj) * ( n0Cj * nEi / fi - nHatCj ) * kronDelta(j,k) - ( n0Ck * nEi / fi ) * Mij + (nHatCj / n0Cj) * prevIter_ijk - (epsj / n0Cj) * doubSum;
	  respErrMat_->SetBinContent(bin, j, dncj_dPeick);

        }
      }
    } //-- End TH2D loop

    return 1;

  } //-- End else

  return 0; //-- This return should never be reached, but is placed here to keep the compiler happy

}

//-- =======================================
//--                kronDelta
//-- =======================================

double DAgostiniUnfold::kronDelta(int i, int j){

  //-- Kronecker delta
  if(i == j) return 1.;
  return 0.;

}

//-- =======================================
//--          makeEfficiencyVector
//-- =======================================

void DAgostiniUnfold::makeEfficiencyVector(){

  //-- Make the detection efficiency vector = Sum[ P( ei | cj ) , { i, 1, Nobs} ]

  for(int y = 1; y <= Ntru; y++){

    double sum = 0.;
    for(int m = 1; m <= Nobs; m++) sum += hresp_->GetBinContent(m, y);
    effVec_.push_back(sum);

  }

}

//-- =======================================
//--                  Unfold
//-- =======================================

int DAgostiniUnfold::Unfold(TH1D * hObserved, int niter, TH1D * hInitialPrior){

  hObs_  = (TH1D*) hObserved->Clone("hresp_");

  //-- If the initial prior is not specified, use the projection onto the truth axis as the starting point
  if(hInitialPrior){
    std::cout<<"Manually setting the starting distribution"<<std::endl;
    SetPrior(hInitialPrior, 0);
  }
  else std::cout<<"Using the projection of the response onto the truth axis as the starting distribution"<<std::endl;

  for(int k = 0; k < niter; k++){

    //-- Check if the prior is set (important for first iteration)
    if(!hPrior_){
      std::cout<<"Prior not set, abort unfolding procedure"<<std::endl;
      return 0;
    }

    //-- Update the prior for subsequent iterations
    if(k > 0) SetPrior(hReco_, k);

    //-- Construct the unfolding matrix
    std::cout<<"Constructing the unfolding matrix for iteration "<<k<<"..."<<std::endl;
    for(int x = 1; x <= Nobs; x++){
      for(int y = 1; y <= Ntru; y++){

	double epsj  = effVec_[y-1];
	if(epsj == 0) continue;

	double PeicLc0L = 0.;
	for(int L = 1; L <= Ntru; L++) PeicLc0L += hresp_->GetBinContent(x, L) * hPrior_->GetBinContent(L);
	if(PeicLc0L == 0) continue;

	double Peicj = hresp_->GetBinContent(x, y);
	double c0j   = hPrior_->GetBinContent(y);

	double Mij = Peicj * c0j / PeicLc0L / epsj;
	Mij_->SetBinContent(x, y, Mij);

      }
    } //-- End unfolding matrix construction loop

    //-- Unfold!
    for(int y = 1; y <= Ntru; y++){

      double cHat_j = 0;
      for(int x = 1; x <= Nobs; x++) cHat_j += Mij_->GetBinContent(x, y) * hObs_->GetBinContent(x);
      hReco_->SetBinContent(y, cHat_j);

    } //-- End unfolding loop

    //-- Set up the error propagation matrices
    
    int obsErrCheck  = computeObsErrMatrix(k);
    int respErrCheck = computeRespErrMatrix(k);

    if(!obsErrCheck){
      std::cout<<"Procedure broke when constructing the observed error propagation matrix"<<std::endl;
      return 0;
    }
    if(!respErrCheck){
      std::cout<<"Procedure broke when constructing the response error propagation matrix"<<std::endl;
      return 0;
    }
    
    //-- Set up the covaciance matrices
    int obsCovCheck  = computeObsCovMatrix(k);
    //int respCovCheck = computeRespCovMatrix(k);

    if(!obsCovCheck){
      std::cout<<"Procedure broke when constructing the observed covariance matrix"<<std::endl;
      return 0;
    }
    //if(!respCovCheck){
    //  std::cout<<"Procedure broke when constructing the response covariance matrix"<<std::endl;
    //  return 0;
    // }

    //-- Make the covariance matrix for the unfolded distribution
    recoCovMat_->Reset();
    recoCovMat_->Add(obsCovMat_);
    recoCovMat_->Add(respCovMat_);

    //Error bars on hReco_ set as the square root of the diagonal elements of the cov matrix for hReco_
    for(int i = 1; i<=Ntru; i++){
      double err = TMath::Sqrt( recoCovMat_->GetBinContent(i, i) );
      hReco_->SetBinError(i, err);
    }
    

  } //-- End iteration loop

  return 1;
}

//-- =======================================
//--                GetHReco
//-- =======================================
TH1D * DAgostiniUnfold::GetHReco(){

  //-- Returns the unfolded distribution as a TH1D object
  return hReco_;

}

//-- =======================================
//--              GetHRecoCov
//-- =======================================

TH2D * DAgostiniUnfold::GetHRecoCov(){

  //-- Returns the unfolded distribution's covariance matrix as a TH2D object
  return recoCovMat_;

}

//-- =======================================
//--               GetObserved
//-- =======================================

TH1D * DAgostiniUnfold::GetObserved(){

  //-- Returns the observed, smeared distribution as a TH1D object
  return hObs_;

}

//-- =======================================
//--               GetResponse
//-- =======================================

TH2D * DAgostiniUnfold::GetResponse(){

  //-- Returns the response matrix as a TH2D object
  return hresp_;

}

//-- =======================================
//--                SetPrior
//-- =======================================

void DAgostiniUnfold::SetPrior(TH1D * hPrior, int iter){

  //-- Sets the Prior distribution that will be used to construct the unfolding matrix for each iteration
  hPrior_ = (TH1D*) hPrior->Clone( Form("hPrior_%i", iter) );

}
