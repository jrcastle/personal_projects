/*
 * =====================================================================================
 *
 *       Filename:  ebye.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/02/2013 16:06:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

{
	const int NCENT_PBPB = 40;
	TH2D * hVn2Dfull[7][NCENT_PBPB];
	TH2D * hVn2Dsub0[7][NCENT_PBPB];
	TH2D * hVn2Dsub1[7][NCENT_PBPB];
	TH2D * hVn2D0v1[7][NCENT_PBPB];

	TH1D * hVnFull[7][NCENT_PBPB];
	TH1D * hVnSub0[7][NCENT_PBPB];
	TH1D * hVnSub1[7][NCENT_PBPB];
    
    //TH1D * hVn1DXsub01[7][NCENT_PBPB];  //JAMES
    //TH1D * hVn1DYsub01[7][NCENT_PBPB];  //JAMES
    //TH1D * hinput;
    //TH2D * hinVsm;
	//TNtupleD * hNtVn[7][NCENT_PBPB];

	for ( int n = VN; n < VN+1; n++ ) {
		for ( int c = 0; c < NCENT_PBPB; c++ ) {
//			cout << " n = " << n << "\tc = " << c << endl;
			hVn2Dfull[n][c] = (TH2D*) f->Get(Form("qwebye/hVn2Dfull_%i_%i", n, c));
			hVn2Dsub0[n][c] = (TH2D*) f->Get(Form("qwebye/hVn2Dsub0_%i_%i", n, c));
			hVn2Dsub1[n][c] = (TH2D*) f->Get(Form("qwebye/hVn2Dsub1_%i_%i", n, c));
			hVn2D0v1[n][c] = (TH2D*) f->Get(Form("qwebye/hVn2D0v1_%i_%i", n, c));

			hVnFull[n][c] = (TH1D*) f->Get(Form("qwebye/hVnFull_%i_%i", n, c));
			hVnSub0[n][c] = (TH1D*) f->Get(Form("qwebye/hVnSub0_%i_%i", n, c));
			hVnSub1[n][c] = (TH1D*) f->Get(Form("qwebye/hVnSub1_%i_%i", n, c));
            
            //JAMES ADDTION!!!!!!!!!---------------------------------
            
            //hVn1DXsub01[n][c] = (TH1D*) f->Get(Form("qwebye/hVn1DXsub01_%i_%i",n,c));
            //hVn1DYsub01[n][c] = (TH1D*) f->Get(Form("qwebye/hVn1DYsub01_%i_%i",n,c));
            
            //END JAMES ADDITION!!!!---------------------------------

			hVn2Dfull[n][c]->GetXaxis()->SetTitle("V_{2,x}^{obs}");
			hVn2Dsub0[n][c]->GetXaxis()->SetTitle("V_{2,x}^{obs}");
			hVn2Dsub1[n][c]->GetXaxis()->SetTitle("V_{2,x}^{obs}");
			hVn2D0v1[n][c]->GetXaxis()->SetTitle("V_{2,x}^{obs,a} - V_{2,x}^{obs,b}");

			hVnFull[n][c]->GetXaxis()->SetTitle("V_{2}^{obs}");
			hVnSub0[n][c]->GetXaxis()->SetTitle("V_{2}^{obs}");
			hVnSub1[n][c]->GetXaxis()->SetTitle("V_{2}^{obs}");

			hVn2Dfull[n][c]->GetYaxis()->SetTitle("V_{2,y}^{obs}");
			hVn2Dsub0[n][c]->GetYaxis()->SetTitle("V_{2,y}^{obs}");
			hVn2Dsub1[n][c]->GetYaxis()->SetTitle("V_{2,y}^{obs}");
			hVn2D0v1[n][c]->GetYaxis()->SetTitle("V_{2,y}^{obs,a} - V_{2,y}^{obs,b}");

		}
	}
    
    //JAMES ADDTION!!!!!!!!!---------------------------------
    
    //hinput = (TH1D*) f->Get("inputV2");
    //hinVsm = (TH2D*) f->Get("inputVSsmear");
    
    //END JAMES ADDITION!!!!---------------------------------
}


