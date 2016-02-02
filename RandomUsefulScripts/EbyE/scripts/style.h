/*
 * =====================================================================================
 *
 *       Filename:  style.h
 *
 *    Description:  useful script
 *
 *        Version:  1.0
 *        Created:  11/12/2012 10:23:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan Wang, 
 *   Organization:  
 *
 * =====================================================================================
 */
 
 void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns,
                          const Int_t rows, const Float_t leftOffset=0.,
                          const Float_t bottomOffset=0.,
                          const Float_t leftMargin=0.2,
                          const Float_t bottomMargin=0.2,
                          const Float_t edge=0.05) {
   if (canv==0) {
      Error("makeMultiPanelCanvas","Got null canvas.");
      return;
   }
   canv->Clear();

   TPad* pad[columns][rows];
   
   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth =
      (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
			   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
      (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			  (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
   
   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
   
   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
			      Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);
	 
         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);
	 
         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}


void SetStyle() {

  TH1::SetDefaultSumw2();
  gROOT->SetStyle("Plain");

  TStyle *MITStyle = new TStyle("MIT-Style","The Perfect Style for Plots ;-)");
  gStyle = MITStyle;

  gStyle->SetTitleSize(0.08,"t");//  gStyle->SetTitleX(0.8);
  gStyle->SetTitleW(0.3);

  // Canvas

  MITStyle->SetCanvasColor     (0);
  MITStyle->SetCanvasBorderSize(10);
  MITStyle->SetCanvasBorderMode(0);

  MITStyle->SetCanvasDefH      (500);
  MITStyle->SetCanvasDefW      (550);
  //  MITStyle->SetCanvasDefX      (100);
  //  MITStyle->SetCanvasDefY      (100);

  // Pads

  MITStyle->SetPadColor       (0);
  //  MITStyle->SetPadBorderSize  (10);
  MITStyle->SetPadBorderMode  (0);
  MITStyle->SetPadBottomMargin(0.16);
  MITStyle->SetPadTopMargin   (0.06);
  MITStyle->SetPadLeftMargin  (0.1);
  MITStyle->SetPadRightMargin (0.14);
  MITStyle->SetPadGridX       (0);
  MITStyle->SetPadGridY       (0);
  MITStyle->SetPadTickX       (1);
  MITStyle->SetPadTickY       (1);

  // Frames

  MITStyle->SetFrameFillStyle ( 0);
  MITStyle->SetFrameFillColor ( 0);
  MITStyle->SetFrameLineColor ( 1);
  MITStyle->SetFrameLineStyle ( 0);
  MITStyle->SetFrameLineWidth ( 1);
  //  MITStyle->SetFrameBorderSize(10);
  MITStyle->SetFrameBorderMode( 0);

  // Histograms

  MITStyle->SetHistLineColor(1);
  MITStyle->SetHistLineStyle(0);
  MITStyle->SetHistLineWidth(3);
  MITStyle->SetNdivisions(505,"X");
  MITStyle->SetNdivisions(505,"Y");
  MITStyle->SetNdivisions(505,"Z");

  // Functions

  MITStyle->SetFuncColor(1);
  MITStyle->SetFuncStyle(0);
  MITStyle->SetFuncWidth(2);

  // Various

//  MITStyle->SetMarkerStyle(20);
//  MITStyle->SetMarkerColor(kBlack);
  MITStyle->SetMarkerSize (1.2);

//  MITStyle->SetTitleSize  (0.1,"t");

  MITStyle->SetTitleSize  (0.070,"X");
  MITStyle->SetTitleOffset(1.0,"X");
  MITStyle->SetTitleFont  (42,"X");
  MITStyle->SetLabelOffset(0.005,"X");
  MITStyle->SetLabelSize  (0.070,"X");
  MITStyle->SetLabelFont  (42   ,"X");

  MITStyle->SetTitleSize  (0.070,"Y");
  MITStyle->SetTitleOffset(1.0,"Y");
  MITStyle->SetTitleFont  (42,"Y");
  MITStyle->SetLabelOffset(0.005,"Y");
  MITStyle->SetLabelSize  (0.070,"Y");
  MITStyle->SetLabelFont  (42   ,"Y");

  MITStyle->SetTitleSize  (0.06,"Z");
  MITStyle->SetTitleOffset(1.640,"Z");
  MITStyle->SetTitleFont  (42,"Z");
  MITStyle->SetLabelOffset(0.005,"Z");
  MITStyle->SetLabelSize  (0.070,"Z");
  MITStyle->SetLabelFont  (42   ,"Z");

  MITStyle->SetStripDecimals(kFALSE);
  MITStyle->SetTitleBorderSize(0);
  MITStyle->SetTitleFillColor(0);
  //  MITStyle->SetTitleAlign(6);

  MITStyle->SetTextSize   (20);
  MITStyle->SetTextFont   (43);

  MITStyle->SetStatFont   (42);
  MITStyle->SetOptStat    (0);

  MITStyle->SetLegendBorderSize(0);
  MITStyle->SetEndErrorSize(0);
  MITStyle->SetErrorX(0);

  MITStyle->SetPalette    (1,0);
  return();
}
void InitHist(TH1        *hist,
	      const char *xtit,
	      const char *ytit  = "Number of Entries",
	      EColor      color = kBlack)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(color);
  hist->SetTitleSize  (0.06,"X");
  hist->SetTitleSize  (0.06,"Y");
  hist->SetTitleSize  (0.06,"Z");
  hist->SetTitleOffset(1.00,"Y");
  hist->SetTitleOffset(1.00,"X");
  hist->SetTitleOffset(1.00,"Z");
  hist->SetLabelOffset(0.008,"Y");
  hist->SetLabelSize  (0.050,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (0.6);
  // Strangely enough this cannot be set anywhere else??
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleFont(42);
  hist->SetTitle("");  
  return;
}

TCanvas* MakeCanvas(const char* name, const char *title)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.15);
  canvas->SetRightMargin    (0.06);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void Style()
{
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetFillColor(0);
	gStyle->SetTitleColor(kBlack);
	gStyle->SetStatColor(0);
	gStyle->SetHistFillColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit();
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetLineWidth(1);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetMarkerSize(2);
	gStyle->SetLabelColor(kBlack);
}

void myStyle(Int_t lStat=0)
{
	//Set gStyle
	int font = 42;
	// From plain
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(10);
	gStyle->SetCanvasColor(10);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetStatColor(10);
	gStyle->SetStatBorderSize(1);
	gStyle->SetLegendBorderSize(1);
	//
	gStyle->SetDrawBorder(0);
	gStyle->SetTextFont(font);
	gStyle->SetStatFont(font);
	gStyle->SetStatFontSize(0.05);
	gStyle->SetStatX(0.97);
	gStyle->SetStatY(0.98);
	gStyle->SetStatH(0.03);
	gStyle->SetStatW(0.3);
	gStyle->SetTickLength(0.02,"y");
	gStyle->SetEndErrorSize(3);
	gStyle->SetLabelSize(0.05,"xyz");
	gStyle->SetLabelFont(font,"xyz"); 
	gStyle->SetLabelOffset(0.01,"xyz");
	gStyle->SetTitleFont(font,"xyz");  
	gStyle->SetTitleOffset(1.0,"xyz");  
	gStyle->SetTitleSize(0.06,"xyz");  
	gStyle->SetMarkerSize(1); 
	gStyle->SetPalette(1,0); 
	if (lStat){
		gStyle->SetOptTitle(1);
		gStyle->SetOptStat(1111);
		gStyle->SetOptFit(1111);
	}
	else {
		gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
	}

}

void grStyle(TGraphErrors* gr, int mkstyle, int color, double mksize=1)
{
	gr->SetMarkerStyle(mkstyle);
	gr->SetMarkerColor(color);
	gr->SetLineColor(color);
	gr->SetMarkerSize(mksize);
}
