char* fname[] = {
	"./test.root", 		// 0, reserved
	"/",			// 1, default
	"/",		// 2, default
	"/", 	// 3, pt and eta
	"/", 	// 4, simu v2 = 0.1
	"/", 		// 5, including pix tracks
	"/", 	// 6, with acc and eff
	"/", 	// 7, with eff
	"/", 	// 8, standard
	"EbyE_pix_eff_subevt1_nt/ebyevn.root", 			// 9, rnd sub event * 
	"EbyE_pix_acc_eff_nt/ebyevn_newbin.root", 		// 10, new bin from here STANDARD run
	"EbyE_pix_eff_pt05_pt15_nt/ebyevn_nt.root", 		// 11, 0.5<pT<1.5
	"EbyE_pix_eff_pt05_pt15_subevt1_nt/ebyevn_nt.root", 	// 12, 0.5<pT<1.5 rnd subevt
	"EbyE_pix_eff_pt15_nt/ebyevn_nt.root", 			// 13, pT>1.5
	"EbyE_pix_eff_pt15_subevt1_nt/ebyevn_nt.root", 		// 14, pT>1.5 rnd subevt
	"EbyE_pPb_eff_nt/ebyevn_nt.root", 			// 15, pPb eta>0 eta<0
	"EbyE_pPb_eff_subevt1_nt/ebyevn_nt.root", 		// 16, pPb rnd subevt
	"EbyE_pPb_eff_nt/ebyevn_nt_Rcumu4_f0.root", 		// 17, pPb v2 shift by cumu v24
	"EbyE_pix_acc_eff_nt/ebyevn_newbin_Rcumu4_f0.root", 	// 18, PbPb v2 shift cumu v24
	"EbyE_pix_eff_subevt1_nt/ebyevn_Rcumu4_f0.root", 	// 19, PbPb v2 shift cumu v24 rnd sub event
	"EbyE_pPb_eff_subevt1_nt/ebyevn_Rcumu4_f0.root", 		// 20, pPb v2 shift cumu v24 rnd subevt
	"EbyE_pix_acc_eff_MatchMult_nt/ebyevn_bin150.root", 		// 21, PbPb match mult
	"EbyE_pPb_eff_MatchMult_nt/ebyevn_nt.root", 			// 22, pPb match mult
	"EbyE_pix_ppreco_acc_eff_nt/ebyevn_nt2.root", 		// 23, ppreco PbPb
	"EbyE_pPb_eff_nt_align/ebyevn_nt2.root", 		// 24, pPb HM remove bad alignment
	"EbyE_Pbp_eff_nt_align/ebyevn_nt2.root", 		// 25, Pbp HM remove bad alignment
	"EbyE_pPb_MB_eff_nt_align/ebyevn_nt2.root", 		// 26, pPb MB remove bad alignment
	"EbyE_Pbp_MB_eff_nt_align/ebyevn_nt2.root", 		// 27, Pbp MB remove bad alignment
	"EbyE_pPb_merge_eff_nt_align/ebyevn_nt2.root", 		// 28, pPb+Pbp 24+25 HM remove bad alignment
	"EbyE_pPb_MB_merge_eff_nt_align/ebyevn_nt2.root", 		// 29, pPb+Pbp 26+27 MB remove bad alignment
	"EbyE_pPb_MB_HM_eff_nt_align/ebyevn_nt2.root", 		// 30, pPb MB+HM 28+29 MB remove bad alignment
	"EbyE_pPb_eff_MatchMult_nt_align/ebyevn_nt2.root", 		// 31, pPb HM force multiplicity remove bad alignment
	"EbyE_Pbp_eff_MatchMult_nt_align/ebyevn_nt2.root", 		// 32, pPb HM force multiplicity remove bad alignment
	"EbyE_pPb_MB_eff_MatchMult_nt_align/ebyevn_nt2.root", 		// 33, pPb MB force multiplicity remove bad alignment
	"EbyE_Pbp_MB_eff_MatchMult_nt_align/ebyevn_nt2.root", 		// 34, pPb MB force multiplicity remove bad alignment
	"EbyE_pix_ppreco_MatchMult_acc_eff_nt/ebyevn_nt2.root", 	// 35, PbPb ppreco force multiplicity remove bad alignment
	"EbyE_PbPb_Hijing/ebyevn_nt2.root", 	// 36, PbPb hijing 2760 GeV |eta|<2.4 stable rnd v2=0.05, fluct_v2=0.02
	"Hijing_GEN_v2_15_fv2_2/ebyevn_nt2.root", 	// 37, PbPb hijing 2760 GeV |eta|<2.4 stable rnd v2=0.15, fluct_v2=0.02
	"hijing_nonflow_v2_15_fv2_2/ebyevn_nt2.root", 	// 38, PbPb hijing 2760 GeV |eta|<2.4 stable with non-flow v2=0.15, fluct_v2=0.02
	"hijing_nonflow_v2_5_fv2_2/ebyevn_nt2.root", 	// 39, PbPb hijing 2760 GeV |eta|<2.4 stable with non-flow v2=0.05, fluct_v2=0.02
	"Toy_m500_v2_15_fv2_2/ebyevn_nt2.root", 	// 40, ToyMC mult=500 with flow v2=0.15, fluct_v2=0.02
	"Toy_m500_v2_05_fv2_2/ebyevn_nt2.root", 	// 41, ToyMC mult=500 with flow v2=0.05, fluct_v2=0.02
	"Toy_m500_f50_v2_15_fv2_2/ebyevn_nt2.root", 	// 42, ToyMC mult=500 mult_sigma=50 with flow v2=0.15, fluct_v2=0.02
	"Toy_m500_f50_v2_05_fv2_2/ebyevn_nt2.root", 	// 43, ToyMC mult=500 mult_sigma=50 with flow v2=0.05, fluct_v2=0.02
	"Toy_m2000_f500_v2_15_fv2_2/ebyevn_nt2.root", 	// 44, ToyMC mult=2000 mult_sigma=500 with flow v2=0.15, fluct_v2=0.02
	"Toy_m2000_f500_v2_05_fv2_2/ebyevn_nt2.root", 	// 45, ToyMC mult=2000 mult_sigma=500 with flow v2=0.05, fluct_v2=0.02
	"Hijing_mod1_v2_15_fv2_2/ebyevn_nt2.root", 	// 46, Hijing b=0-10 without nonflow v2=0.15, fluct_v2=0.02
	"Hijing_mod1_v2_05_fv2_2/ebyevn_nt2.root", 	// 47, Hijing b=0-10 without nonflow v2=0.05, fluct_v2=0.02
	"CastleEbyE_m500_v2_15_fv2_2/CastleEbyE.root",    // 48, Castle Attempts
	"CastleEbyE_m500_v2_05_fv2_2/CastleEbyE.root",    // 49
	"CastleEbyE_m2000_v2_15_fv2_2/CastleEbyE.root",   // 50
	"CastleEbyE_m2000_v2_05_fv2_2/CastleEbyE.root",   // 51
	"CastleEbyE_m500_fm_50_v2_15_fv2_2/CastleEbyE.root",    // 52
	"CastleEbyE_m500_fm_50_v2_05_fv2_2/CastleEbyE.root",    // 53
	"CastleEbyE_m2000_fm_500_v2_15_fv2_2/CastleEbyE.root",   // 54
	"CastleEbyE_m2000_fm_500_v2_05_fv2_2/CastleEbyE.root",   // 55
	"CastleEbyE_Glauber/CastleEbyE.root",                    //56
	"CastleEbyE_PbPb/CastleEbyE.root",                        //57
	"CastleEbyE_pPb/CastleEbyE.root",                        //58
	"PbPb_2011/CastleEbyE.root",                             //59
	"PbPb_2015/CastleEbyE.root"                             //60
};


extern int s1;

using namespace std;
void ffwrite(TH1D *h) {
//	cout << "!!!!!!!" << endl;
//	h->Print();
//	cout << Form("txt/%s/%s.txt", fname[s1], h->GetName()) << endl;
	ofstream of(Form("txt/%s/%s.txt", txtfname[s1],h->GetName()));
	of << h->GetNbinsX() << endl;
	of.precision(10);
	for ( int i = 1; i <= h->GetNbinsX(); i++ ) {
		of << h->GetBinLowEdge(i) << "\t" << h->GetBinLowEdge(i+1) << "\t" << h->GetBinContent(i) << "\t" << h->GetBinError(i) << endl;
	}
	return;
}

char* txtfname[] = {
	"tmp/",		// 0 reserved
	"",	// 1
	"",	// 2
	"",	// 3
	"",	// 4
	"",	// 5
	"",	// 6
	"",	// 7
	"",	// 8
	"EbyE_pix_eff_subevt1_nt/",	// 9
	"EbyE_pix_acc_eff_nt/",	// 10
	"EbyE_pix_eff_pt05_pt15_nt/",	// 11
	"EbyE_pix_eff_pt05_pt15_subevt1_nt/",	// 12
	"EbyE_pix_eff_pt15_nt/",	// 13
	"EbyE_pix_eff_pt15_subevt1_nt/",	// 14
	"EbyE_pPb_eff_nt/",	// 15
	"EbyE_pPb_eff_subevt1_nt/",	// 16
	"EbyE_pPb_eff_nt/weighted",	// 17
	"EbyE_pix_acc_eff_nt/weighted",	// 18
	"EbyE_pix_eff_subevt1_nt/weighted",	// 19
	"EbyE_pPb_eff_subevt1_nt/weighted",	// 20
	"EbyE_pix_acc_eff_MatchMult_nt/",	// 21
	"EbyE_pPb_eff_MatchMult_nt/",		// 22
	"EbyE_pix_ppreco_acc_eff_nt/",		// 23
	"EbyE_pPb_eff_nt_align/",		// 24
	"EbyE_Pbp_eff_nt_align/",		// 25
	"EbyE_pPb_MB_eff_nt_align/",		// 26
	"EbyE_Pbp_MB_eff_nt_align/",		// 27
	"EbyE_pPb_merge_eff_nt_align/",		// 28
	"EbyE_pPb_MB_merge_eff_nt_align/",	// 29
	"EbyE_pPb_MB_HM_eff_nt_align/",		// 30
	"EbyE_pPb_eff_MatchMult_nt_align/",	// 31
	"EbyE_Pbp_eff_MatchMult_nt_align/",	// 32
	"EbyE_pPb_MB_eff_MatchMult_nt_align/",	// 33
	"EbyE_Pbp_MB_eff_MatchMult_nt_align/",	// 34
	"EbyE_pix_ppreco_MatchMult_acc_eff_nt/",	// 35
	"EbyE_PbPb_Hijing/",	// 36
	"Hijing_GEN_v2_15_fv2_2/",	// 37
	"hijing_nonflow_v2_15_fv2_2/",	// 38
	"hijing_nonflow_v2_5_fv2_2/",	// 39
	"Toy_m500_v2_15_fv2_2/",	// 40
	"Toy_m500_v2_05_fv2_2/",	// 41
	"Toy_m500_f50_v2_15_fv2_2/",	// 42
	"Toy_m500_f50_v2_05_fv2_2/",	// 43
	"Toy_m2000_f500_v2_15_fv2_2/",	// 44
	"Toy_m2000_f500_v2_05_fv2_2/",	// 45
	"Hijing_mod1_v2_15_fv2_2/",	// 46
	"Hijing_mod1_v2_05_fv2_2/",	// 47
	"CastleEbyE_m500_v2_15_fv2_2",    // 48, Castle Attempts
	"CastleEbyE_m500_v2_05_fv2_2",    // 49
	"CastleEbyE_m2000_v2_15_fv2_2",   // 50
	"CastleEbyE_m2000_v2_05_fv2_2",   // 51
	"CastleEbyE_m500_fm_50_v2_15_fv2_2",    // 52
	"CastleEbyE_m500_fm_50_v2_05_fv2_2",    // 53
	"CastleEbyE_m2000_fm_500_v2_15_fv2_2",   // 54
	"CastleEbyE_m2000_fm_500_v2_05_fv2_2",   // 55
	"CastleEbyE_Glauber",                    //56
	"CastleEbyE_PbPb",                       //57
	"CastleEbyE_pPb",                        //58
	"PbPb_2011",                             //59
	"PbPb_2015"                              //60
};

	string cent_pPb[15] = {
		string("N_{off}>=350"),     // 0
		string("300<=N_{off}<350"), // 1
		string("260<=N_{off}<300"), // 2
		string("220<=N_{off}<260"), // 3
		string("185<=N_{off}<220"), // 4
		string("150<=N_{off}<185"), // 5
		string("120<=N_{off}<150"), // 6
		string("100<=N_{off}<120"), // 7
		string("80<=N_{off}<100"), // 8
		string("60<=N_{off}<80"), // 9
		string("50<=N_{off}<60"), // 10
		string("40<=N_{off}<50"),  // 11
		string("30<=N_{off}<40"),  // 12
		string("20<=N_{off}<30"),  // 13
		string("0<=N_{off}<20"),  // 14
	};

