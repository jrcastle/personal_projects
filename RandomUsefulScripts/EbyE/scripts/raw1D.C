/*
 * =====================================================================================
 *
 *       Filename:  raw1D.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/18/2013 10:36:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

{
    
//	int s1 = 18;
//	int VN = 4;
	int SAVE = 1;
#include "label.h"
#include "style.h"
	SetStyle();
	gStyle->SetOptTitle(0);
	TFile* f = new TFile(Form("~/root/macros/EbyE/macros/data/%s", fname[s1]));
	gROOT->Macro("ebye.C");


	TH1D * h1D[20];

	if ( SAVE ) {
		TFile * fsave = new TFile(Form("~/root/macros/EbyE/macros/txt/%s/weighted/raw1D%i.root", txtfname[s1], VN), "recreate");
	}
	for ( int c = 0; c < 20; c++ ) {
		cout << "!! processing Cent = " << c << endl;
		h1D[c] = (TH1D*)hVnFull[VN][2*c]->Clone(Form("h1D_%i_%i", VN, c));
		h1D[c]->Add(hVnFull[VN][2*c+1]);

		if ( SAVE ) {
			fsave->cd();
			h1D[c]->Write();
		}
	}
}
