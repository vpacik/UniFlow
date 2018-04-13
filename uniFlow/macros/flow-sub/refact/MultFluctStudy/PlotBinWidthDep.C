/* PlotBinWidthDep
 *
 * Plot dependence of extracted v2 on N_{RFP} bin width produced by CompareMultFluctuationTest.C macro
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TColor.h"
#include "TLine.h"


void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
void SetCustomPalette();
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);


void PlotBinWidthDep()
{
  const Int_t iNumBinWidths = 6;
  Int_t iBinWidths[iNumBinWidths] = {1,5,10,30,50,75};

  const Int_t iNumMethods = 3;
  TString sMethods[iNumMethods] = {"GF_eventweighted","GF_noneventweighted","SP_nonscaled_noneventweighted"/*,"SP_scaled_noneventweighted"*/};
  TString sMethodLabels[iNumMethods] = {"scaled vect. with weights", "scaled vect. w/o weights", "not scaled vect. w/o weights"};


  // TLegend* leg = new TLegend(0.12,0.7,0.3,0.85);
  TLegend* leg = new TLegend(0.6,0.05,0.95,0.4);
  leg->SetBorderSize(0);
  leg->SetFillColorAlpha(0.,0.);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);

  TString sIntputFileName = "output.root";
  TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/multiplicity-fluctuations/CompareMultFluctuationTest/pp_gap08/";

  TString sGap = "gap08";
  TString sSpecies = "Charged";

  Double_t dXmin = 0.0;
  Double_t dXmax = 7.0;

  SetCustomPalette();
  Int_t nPnt  = iNumBinWidths+1;
  Int_t nnCol = gStyle->GetNumberOfColors();

  TLine* lineUnity = new TLine(dXmin,1.0,dXmax,1.0);
  lineUnity->SetLineColor(kGray+2);
  lineUnity->SetLineStyle(kDashed);

  // Preparring canvases / pads
  TCanvas* canMethods[iNumMethods];
  TPad* padMethods_Main[iNumMethods];
  TPad* padMethods_Ratio[iNumMethods];
  for(Int_t iMethod(0); iMethod < iNumMethods; ++iMethod)
  {
    TCanvas* canDiff = new TCanvas(Form("canDiff_method_%d",iMethod),Form("canDiff"),500,1200);
    canDiff->cd();
    TPad* padMain = new TPad("padMain","padMain", 0, 0.5, 1, 1.0);
    padMain->SetBottomMargin(0.0);
    padMain->SetRightMargin(0.03);
    padMain->SetLeftMargin(0.13);
    padMain->Draw();
    padMain->cd();
    TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(dXmin,-0.02,dXmax,0.35);
    frame_canDiff_1->SetTitle(Form("%s v_{2}{2,%s} (0-20%% V0A) ; ; v_{2}{2}",sSpecies.Data(),sGap.Data()));
    frame_canDiff_1->SetTitleFont(43,"X");
    frame_canDiff_1->SetTitleSize(18,"X");
    frame_canDiff_1->SetTitleOffset(4.3,"X");
    frame_canDiff_1->SetLabelFont(43,"X");
    frame_canDiff_1->SetLabelSize(18,"X");
    // frame_canDiff_1->SetTitleFont(43,"Y");
    // frame_canDiff_1->SetTitleOffset(2.2,"Y");

    canDiff->cd();
    TPad* padRatio = new TPad("padRatio","padRatio", 0, 0.0, 1, 0.5);
    padRatio->SetTopMargin(0.0);
    padRatio->SetBottomMargin(0.2);
    padRatio->SetRightMargin(0.03);
    padRatio->SetLeftMargin(0.13);
    padRatio->Draw();
    padRatio->cd();
    TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(dXmin,0.5,dXmax,1.5);
    frame_canDiff_2->SetTitle(Form("; p_{T} (GeV/c); wide / unit bin  "));
    frame_canDiff_2->SetNdivisions(510,"Y");
    frame_canDiff_2->SetTitleFont(43,"XY");
    frame_canDiff_2->SetTitleSize(18,"XY");
    frame_canDiff_2->SetTitleOffset(2.2,"X");
    frame_canDiff_2->SetTitleOffset(2.2,"Y");
    frame_canDiff_2->SetLabelFont(43,"XY");
    frame_canDiff_2->SetLabelSize(18,"XY");
    lineUnity->Draw("same");

    canMethods[iMethod] = canDiff;
    padMethods_Main[iMethod] = padMain;
    padMethods_Ratio[iMethod] = padRatio;
  }


  TH1D* hMerged[iNumMethods];
  TH1D* hUnit[iNumMethods];
  for(Int_t iWidth(0); iWidth < iNumBinWidths; ++iWidth)
  {
    Int_t iMultBinWidth = iBinWidths[iWidth];
    TString sInputFilePath = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/multiplicity-fluctuations/CompareMultFluctuationTest/pp_gap08/%d/",iMultBinWidth);

    TFile* fileIn = TFile::Open(Form("%s%s",sInputFilePath.Data(),sIntputFileName.Data()),"READ");
    if(!fileIn) { printf("ERROR: Input file not found!\n"); return; }

    // color index
    Int_t idx = iWidth * Float_t(nnCol-1) / (nPnt-1);

    TH1D* hWeight_GF_eventweighted = (TH1D*) fileIn->Get("hWeight_GF_eventweighted");
    if(!hWeight_GF_eventweighted) { printf("ERROR : 'hWeight_GF_eventweighted' not found!\n"); return; }
    StyleHist(hWeight_GF_eventweighted, gStyle->GetColorPalette(idx), kOpenCircle);
    leg->AddEntry(hWeight_GF_eventweighted,Form("Bin width %d",iMultBinWidth), "p");


    TH1D* hWeight_GF_noneventweighted = (TH1D*) fileIn->Get("hWeight_GF_noneventweighted");
    if(!hWeight_GF_noneventweighted) { printf("ERROR : 'hWeight_GF_noneventweighted' not found!\n"); return; }
    StyleHist(hWeight_GF_noneventweighted, gStyle->GetColorPalette(idx), kOpenCircle);

    TH1D* hWeight_SP_nonscaled_noneventweighted = (TH1D*) fileIn->Get("hWeight_SP_nonscaled_noneventweighted");
    if(!hWeight_SP_nonscaled_noneventweighted) { printf("ERROR : 'hWeight_SP_nonscaled_noneventweighted' not found!\n"); return; }
    StyleHist(hWeight_SP_nonscaled_noneventweighted, gStyle->GetColorPalette(idx), kOpenCircle);


    if(iWidth == 0)
    {
      // saving unit bin results

      StyleHist(hWeight_GF_eventweighted, gStyle->GetColorPalette(idx), kOpenSquare);
      StyleHist(hWeight_GF_noneventweighted, gStyle->GetColorPalette(idx), kOpenSquare);
      StyleHist(hWeight_SP_nonscaled_noneventweighted, gStyle->GetColorPalette(idx), kOpenSquare);
      hUnit[0] = hWeight_GF_eventweighted;
      hUnit[1] = hWeight_GF_noneventweighted;
      hUnit[2] = hWeight_SP_nonscaled_noneventweighted;

      // loading integrated / merged resulst (once)
      TH1D* hMerged_GF_eventweighted = (TH1D*) fileIn->Get("hMerged_GF_eventweighted");
      if(!hMerged_GF_eventweighted) { printf("ERROR : 'hMerged_GF_eventweighted' not found!\n"); return; }
      StyleHist(hMerged_GF_eventweighted, gStyle->GetColorPalette(nnCol-1), kOpenCircle);
      hMerged[0] = hMerged_GF_eventweighted;

      TH1D* hMerged_GF_noneventweighted = (TH1D*) fileIn->Get("hMerged_GF_noneventweighted");
      if(!hMerged_GF_noneventweighted) { printf("ERROR : 'hMerged_GF_noneventweighted' not found!\n"); return; }
      StyleHist(hMerged_GF_noneventweighted, gStyle->GetColorPalette(nnCol-1), kOpenCircle);
      hMerged[1] = hMerged_GF_noneventweighted;

      TH1D* hMerged_SP_nonscaled_noneventweighted = (TH1D*) fileIn->Get("hMerged_SP_nonscaled_noneventweighted");
      if(!hMerged_SP_nonscaled_noneventweighted) { printf("ERROR : 'hMerged_SP_nonscaled_noneventweighted' not found!\n"); return; }
      StyleHist(hMerged_SP_nonscaled_noneventweighted, gStyle->GetColorPalette(nnCol-1), kOpenCircle);
      hMerged[2] = hMerged_SP_nonscaled_noneventweighted;
    }

    padMethods_Main[0]->cd();
    hWeight_GF_eventweighted->DrawCopy("same");
    if(iWidth == 0) { hMerged[0]->DrawCopy("same"); }
    padMethods_Ratio[0]->cd();
    if(iWidth != 0)
    {
      TH1D* hRatio_GF_eventweighted = DivideHistos(hWeight_GF_eventweighted,hUnit[0],kTRUE);
      hRatio_GF_eventweighted->DrawCopy("hist p e1 same");
    }

    padMethods_Main[1]->cd();
    hWeight_GF_noneventweighted->DrawCopy("same");
    if(iWidth == 0) { hMerged[1]->DrawCopy("same"); }
    padMethods_Ratio[1]->cd();
    if(iWidth != 0)
    {
      TH1D* hRatio_GF_noneventweighted = DivideHistos(hWeight_GF_noneventweighted,hUnit[0],kTRUE);
      hRatio_GF_noneventweighted->DrawCopy("hist p e1 same");
    }

    padMethods_Main[2]->cd();
    hWeight_SP_nonscaled_noneventweighted->DrawCopy("same");
    if(iWidth == 0) { hMerged[2]->DrawCopy("same"); }
    padMethods_Ratio[2]->cd();
    if(iWidth != 0)
    {
      TH1D* hRatio_SP_nonscaled_noneventweighted = DivideHistos(hWeight_SP_nonscaled_noneventweighted,hUnit[0],kTRUE);
      hRatio_SP_nonscaled_noneventweighted->DrawCopy("hist p e1 same");
    }
  }

  // making ratios wrt to intefrated


  leg->AddEntry(hMerged[0],"Integrated","p");

  for(Int_t iMethod(0); iMethod < iNumMethods; ++iMethod)
  {
    padMethods_Main[iMethod]->cd();
    leg->SetHeader(Form("%s",sMethodLabels[iMethod].Data()));
    leg->Draw();

    TH1D* hRatio = DivideHistos(hMerged[iMethod],hUnit[iMethod],kTRUE);
    padMethods_Ratio[iMethod]->cd();
    hRatio->DrawCopy("hist p e1 same ");

    canMethods[iMethod]->SaveAs(Form("%s/BinWidthComp_method_%d.pdf",sOutputPath.Data(),iMethod),"pdf");
  }

  return;
}
// ==================================================================================================================
void StyleHist(TH1* hist, Color_t color, Style_t markerStyle)
{
  if(!hist) { printf("ERROR-DrawHist: Hist does not found.\n"); return; }
  hist->SetLineColor(color);
  // hist->SetLineStyle(color);
  // hist->SetLineStyle(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  return;
};
// ==================================================================================================================
void SetCustomPalette()
{
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,	0.6250, 0.7500, 0.8750, 1.0000};


  // Rain Bow
  // case 55:
  Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
  Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
  Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};

  // Bird
  //case 57:
  // Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  // Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  // Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};

  // Blue Green Yellow
  //case 71:
  // Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
  // Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
  // Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};

  // Solar
  // case 100:
  // Double_t red[9]   = { 99./255., 116./255., 154./255., 174./255., 200./255., 196./255., 201./255., 201./255., 230./255.};
  // Double_t green[9] = {  0./255.,   0./255.,   8./255.,  32./255.,  58./255.,  83./255., 119./255., 136./255., 173./255.};
  // Double_t blue[9]  = {  5./255.,   6./255.,   7./255.,   9./255.,   9./255.,  14./255.,  17./255.,  19./255.,  24./255.};

  // Viridis
  // case 112:
  // Double_t red[9]   = { 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
  // Double_t green[9] = {  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
  // Double_t blue[9]  = { 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};

  // Cividis
  // case 113:
  // Double_t red[9]   = {  0./255.,   5./255.,  65./255.,  97./255., 124./255., 156./255., 189./255., 224./255., 255./255.};
  // Double_t green[9] = { 32./255.,  54./255.,  77./255., 100./255., 123./255., 148./255., 175./255., 203./255., 234./255.};
  // Double_t blue[9]  = { 77./255., 110./255., 107./255., 111./255., 120./255., 119./255., 111./255.,  94./255.,  70./255.};

  Int_t pal = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
  const Int_t nCol = 255;
  Int_t colors[nCol];
  for (int i=0; i<nCol; i++) colors[i] =pal+i;

  gStyle->SetPalette(nCol,colors);
}
// ==================================================================================================================
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* ratio = (TH1D*) nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
  ratio->Reset();

  Double_t dContNom = 0, dErrNom = 0;
  Double_t dContDenom = 0, dErrDenom = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < binsDenom+1; iBin++)
  {
    if(iBin > binsNom) break;

    dContNom = nom->GetBinContent(iBin);
    dErrNom = nom->GetBinError(iBin);
    dContDenom = denom->GetBinContent(iBin);
    dErrDenom = denom->GetBinError(iBin);

    if(dContDenom == 0.0) continue;

    dContRatio =  dContNom / dContDenom;
    dErrRatio = TMath::Power(dErrNom/dContDenom, 2) + TMath::Power( dErrDenom*dContNom/(dContDenom*dContDenom), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCor) dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBin,dContRatio);
    ratio->SetBinError(iBin,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
