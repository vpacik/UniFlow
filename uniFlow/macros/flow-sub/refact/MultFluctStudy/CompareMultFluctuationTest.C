/* CompareMultFluctuationTest
 *
 * Compare raw (unsubtracted) result obtained from various flow-calculations methods.
 *
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

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);
void SetCustomPalette();

// Color_t colors[] = {kGreen+2, kSpring+2, kCyan+2, kBlue+2, kMagenta+2, kRed-2, kRed+2, kOrange+2, kYellow+2, kOrange-4};
// Color_t colors[] = {kGreen+2, kGreen, kSpring+2, kSpring, kCyan+2, kCyan, kBlue+2, kBlue, kMagenta+2, kRed-2, kRed+2, kOrange+2, kYellow+2, kOrange-4,};
Color_t colors[] = {kRed, kBlue+1, kGreen+2, kOrange-3};


void CompareMultFluctuationTest()
{
  const Int_t iNumMethods = 3;
  TString sMethods[iNumMethods] = {"GF_eventweighted","GF_noneventweighted","SP_nonscaled_noneventweighted"/*,"SP_scaled_noneventweighted"*/};
  TString sMethodLabels[iNumMethods] = {"scaled vect. with weights", "scaled vect. w/o weights", "not scaled vect. w/o weights"};

  TString sSpecies = "Charged";
  TString sGap = "gap00";

  Int_t iMultMin = 0;
  Int_t iMultBinWidth = 1;
  Int_t iMultMax = 150/iMultBinWidth;

  // file with correlations in 0-20% (V0A) vs. N_RFP
  TString sInputPathUnit = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-NchRFP-gap00/output_vn_%d/",iMultBinWidth);
  // file with correlations in 0-20% (V0A)  (i.e vs centrality percentile)
  TString sInputPathMerged = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-NchRFP-gap00/output_vn_merged/";

  TString sOutputPath = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/CompareMultFluctuationTest/pPb_gap00/%d/",iMultBinWidth);
  TString sOutputFileName = "output.root";
  // ##########################################################################################################
  gSystem->mkdir(sOutputPath.Data(),kTRUE);

  TFile* fileOut = TFile::Open(Form("%s%s",sOutputPath.Data(), sOutputFileName.Data()),"RECREATE");
  if(!fileOut) { printf("ERROR : Output file not found!\n"); return; }

  // gStyle->SetPalette(1);
  SetCustomPalette();

  Double_t dYmin = -10.0;
  Double_t dYmax = 20.0;
  Double_t dXmin = 0.0;
  Double_t dXmax = 7.0;

  TLine* lineUnity = new TLine(dXmin,1.0,dXmax,1.0);
  lineUnity->SetLineColor(kGray+2);
  lineUnity->SetLineStyle(kDashed);

  TLegend* legRatio = new TLegend(0.12,0.75,0.4,0.89);
  legRatio->SetBorderSize(0);
  legRatio->SetFillColorAlpha(0.,0.);
  legRatio->SetTextFont(42);
  legRatio->SetTextSize(0.03);
  legRatio->SetHeader(Form("Bin width %d",iMultBinWidth));

  TCanvas* canRatio = new TCanvas("canRatio","canRatio",800,800);
  TH1* frame_canRatio = (TH1*) gPad->DrawFrame(dXmin,0.8,dXmax,1.2);
  legRatio->Draw();
  lineUnity->Draw("same");

  TCanvas* canCompMethod = new TCanvas("canCompMethod","canCompMethod",1600,800);
  canCompMethod->Divide(2,1);
  canCompMethod->cd(1);
  TH1* frame_canCompMethod_1 = (TH1*) gPad->DrawFrame(dXmin,0.0,dXmax,0.5);
  frame_canCompMethod_1->SetTitle(Form("%s v_{2}{2,%s} (0-20%% V0A) ; p_{T} (GeV/c); v_{2}{2}",sSpecies.Data(),sGap.Data()));
  frame_canCompMethod_1->SetTitleOffset(1.45,"Y");
  canCompMethod->cd(2);
  TH1* frame_canCompMethod_2 = (TH1*) gPad->DrawFrame(dXmin,0.0,dXmax,0.5);
  frame_canCompMethod_2->SetTitle(Form("%s v_{2}{2,%s} (0-20%% V0A) ; p_{T} (GeV/c); v_{2}{2}",sSpecies.Data(),sGap.Data()));
  frame_canCompMethod_2->SetTitleOffset(1.45,"Y");

  TLegend* legMethods = new TLegend(0.12,0.75,0.4,0.89);
  legMethods->SetBorderSize(0);
  legMethods->SetFillColorAlpha(0.,0.);
  legMethods->SetTextFont(42);
  legMethods->SetTextSize(0.03);
  legMethods->SetHeader(Form("Integrated 0-20%%"));

  TLegend* legMethodsWeight = new TLegend(0.12,0.75,0.4,0.89);
  legMethodsWeight->SetBorderSize(0);
  legMethodsWeight->SetFillColorAlpha(0.,0.);
  legMethodsWeight->SetTextFont(42);
  legMethodsWeight->SetTextSize(0.03);
  legMethodsWeight->SetHeader(Form("Weighted mean (bin width %d)",iMultBinWidth));

  for(Int_t iMethod(0); iMethod < iNumMethods; ++iMethod)
  {
    TLegend* leg = new TLegend(0.15,0.70,0.43,0.89);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0.,0.);
    leg->SetHeader(Form("%s",sMethodLabels[iMethod].Data()));
    leg->SetTextFont(42);
    leg->SetTextSize(0.036);

    TLegend* legBins[3];
    for(Int_t iLeg(0); iLeg < 3; ++iLeg)
    {
      legBins[iLeg] = new TLegend(0.0+iLeg*0.33 ,0.0,(iLeg+1)*0.32, 1.0);
      legBins[iLeg]->SetBorderSize(0);
      legBins[iLeg]->SetFillColorAlpha(0.,0.);
      legBins[iLeg]->SetTextFont(42);
      legBins[iLeg]->SetTextSize(0.055);
    }

    TCanvas* canDiff = new TCanvas(Form("canDiff"),Form("canDiff"),500,1200);
    canDiff->cd();
    TPad* padMain = new TPad("padMain","padMain", 0, 0.3, 1, 1.0);
    padMain->SetBottomMargin(0.0);
    padMain->SetRightMargin(0.03);
    padMain->SetLeftMargin(0.13);
    padMain->Draw();
    padMain->cd();
    TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(dXmin,-0.45,dXmax,1);
    frame_canDiff_1->SetTitle(Form("%s v_{2}{2,%s} (0-20%% V0A) ; ; v_{2}{2}",sSpecies.Data(),sGap.Data()));
    frame_canDiff_1->SetTitleFont(43,"X");
    frame_canDiff_1->SetTitleSize(18,"X");
    frame_canDiff_1->SetTitleOffset(4.3,"X");
    frame_canDiff_1->SetLabelFont(43,"X");
    frame_canDiff_1->SetLabelSize(18,"X");
    // frame_canDiff_1->SetTitleFont(43,"Y");
    // frame_canDiff_1->SetTitleOffset(2.2,"Y");

    canDiff->cd();
    TPad* padRatio = new TPad("padRatio","padRatio", 0, 0.0, 1, 0.3);
    padRatio->SetTopMargin(0.0);
    padRatio->SetBottomMargin(0.25);
    padRatio->SetRightMargin(0.03);
    padRatio->SetLeftMargin(0.13);
    padRatio->Draw();
    padRatio->cd();
    TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(dXmin,0.85,dXmax,1.07);
    frame_canDiff_2->SetTitle(Form("; p_{T} (GeV/c); mean / integrated   "));
    frame_canDiff_2->SetNdivisions(505,"Y");
    frame_canDiff_2->SetTitleFont(43,"XY");
    frame_canDiff_2->SetTitleSize(18,"XY");
    frame_canDiff_2->SetTitleOffset(4.3,"X");
    frame_canDiff_2->SetTitleOffset(2.2,"Y");
    frame_canDiff_2->SetLabelFont(43,"XY");
    frame_canDiff_2->SetLabelSize(18,"XY");
    lineUnity->Draw("same");

    TFile* fileUnit = OpenFile(sInputPathUnit+sMethods[iMethod]+"/Processed.root","READ"); if(!fileUnit) return;
    TFile* fileMerged = OpenFile(sInputPathMerged+sMethods[iMethod]+"/Processed.root","READ"); if(!fileMerged) return;

    TH1D* histMerged = LoadHisto(Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),0),fileMerged); if(!histMerged) return;
    StyleHist(histMerged,kRed,kFullCircle);
    histMerged->SetLineWidth(2);
    leg->AddEntry(histMerged,"integrated","p");

    TProfile* p_weight_sigma = new TProfile("p_weight_sigma","p_weight_sigma",histMerged->GetXaxis()->GetNbins(),histMerged->GetXaxis()->GetXbins()->GetArray());
    if(!p_weight_sigma) { printf("Profile does not exists!\n"); return; }
    p_weight_sigma->Sumw2();
    StyleHist(p_weight_sigma,kViolet-1,kOpenSquare);

    TProfile* p_weight_average = new TProfile("p_weight_average","p_weight_average",histMerged->GetXaxis()->GetNbins(),histMerged->GetXaxis()->GetXbins()->GetArray());
    if(!p_weight_average) { printf("Profile does not exists!\n"); return; }
    p_weight_average->Sumw2();
    StyleHist(p_weight_average,kMagenta+2,kOpenSquare);

    // hist for single pt bin distribution of v2 from differrent N_RFP bins
    const Int_t iNumPtBins = histMerged->GetNbinsX()+2;
    TH1D* histVnDist[iNumPtBins];
    TH1D* histVnDistWeighted[iNumPtBins];
    for(Int_t pt(0); pt < iNumPtBins; ++pt)
    {
      histVnDist[pt] = new TH1D(Form("histVnDist_pt%d",pt),Form("v_{2}{2,%s} %g < p_{T} < %g GeV/c; v_{2}{2,%s}; Counts",sGap.Data(),histMerged->GetBinLowEdge(pt+1),histMerged->GetBinLowEdge(pt+2),sGap.Data()), 140, -0.4,1.);
      histVnDistWeighted[pt] = new TH1D(Form("histVnDistWeighted_pt%d",pt),Form("v_{2}{2,%s} (1/#sigma^{2} weight) %g < p_{T} < %g GeV/c; v_{2}{2,%s}; Counts (weighted)",sGap.Data(),histMerged->GetBinLowEdge(pt+1),histMerged->GetBinLowEdge(pt+2),sGap.Data()), 140, -0.4,1.);
    }

    TCanvas* canDist = new TCanvas("canDist","canDist",800,800);


    Int_t nPnt  = iMultMax-iMultMin;
    Int_t nnCol = gStyle->GetNumberOfColors();

    Int_t iLegInx = 0;
    for(Int_t iMult(iMultMin); iMult < iMultMax; ++iMult)
    {
      TString sHistoName = Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),iMult);

      TH1D* histUnit = LoadHisto(sHistoName,fileUnit); if(!histUnit) return;

      Int_t idx = iMult * Float_t(nnCol-1) / (nPnt-1);
      StyleHist(histUnit, gStyle->GetColorPalette(idx), kOpenCircle);

      padMain->cd();
      histUnit->DrawCopy("hist p x0 same");
      if(iMult == iMultMin) { leg->AddEntry(histUnit,Form("N_{RFP} bins (width %d)",iMultBinWidth),"p"); }
      legBins[iLegInx]->AddEntry(histUnit,Form("%d",iMult*iMultBinWidth),"p");
      if((iMult - iMultMin) % 50 == 49) { iLegInx++; }

      for(Int_t pt(1); pt < histUnit->GetNbinsX()+1; ++pt)
      {
        Double_t bin = histUnit->GetBinCenter(pt);
        Double_t dValue = histUnit->GetBinContent(pt);
        Double_t dError = histUnit->GetBinError(pt);

        if(dError <= 0.0001) continue;
        if(dValue >= 8.0 && dError >= 8.0) continue;

        Double_t dWeight = TMath::Power(dError,-2.0);
        p_weight_average->Fill(bin,dValue,1.0);
        p_weight_sigma->Fill(bin,dValue,dWeight);
        if(iMethod == 0 && pt == histUnit->GetNbinsX() )
        {
          // printf("pT %d: Val %f | weight %f\n",pt,dValue,dWeight);
          // printf("%f +- %f \n", dValue, dError);
          if(dError <= 0.0000001) printf("ZERO \n");
        }

        // filling vn dist
        histVnDist[pt-1]->Fill(dValue);
        histVnDistWeighted[pt-1]->Fill(dValue, TMath::Power(dError,-2.0));
      }
    }

    leg->AddEntry(p_weight_sigma,"1/#sigma^{2} weighted mean","p");
    // leg->AddEntry(p_weight_average,"arithmeric mean","p");

    // making histos out of profiles (for "safety" reason)
    TH1D* h_weight_average = (TH1D*) p_weight_average->ProjectionX();
    TH1D* h_weight_sigma = (TH1D*) p_weight_sigma->ProjectionX();
    StyleHist(h_weight_average,kBlack,kOpenSquare);
    StyleHist(h_weight_sigma,kViolet-1,kOpenSquare);
    h_weight_sigma->SetLineWidth(2);

    padMain->cd();
    histMerged->DrawCopy("same");
    // h_weight_average->DrawCopy("same");
    h_weight_sigma->DrawCopy("same");
    leg->Draw();

    padRatio->cd();
    TH1D* hRatioSigma = DivideHistos(h_weight_sigma,histMerged);
    TH1D* hRatioAverage = DivideHistos(h_weight_average,histMerged);
    hRatioSigma->DrawCopy("hist p e1 same");
    // hRatioAverage->DrawCopy("hist p e1 same");

    canDiff->SaveAs(Form("%s/%s_%s.pdf",sOutputPath.Data(),sSpecies.Data(),sMethods[iMethod].Data()),"pdf");

    canRatio->cd();
    StyleHist(hRatioSigma,colors[iMethod],kOpenSquare);
    hRatioSigma->DrawCopy("hist p e1 same");
    legRatio->AddEntry(hRatioSigma,sMethodLabels[iMethod].Data(),"p");


    TCanvas* canLegend = new TCanvas("canLegend","canLegend",200,1000);
    canLegend->cd();
    for(Int_t iLeg(0); iLeg < 3; ++iLeg) legBins[iLeg]->Draw();
    canLegend->SaveAs(Form("%s/%s_%s.pdf",sOutputPath.Data(),sSpecies.Data(),"BinsLegend"),"pdf");

    canCompMethod->cd(2);
    StyleHist(histMerged,colors[iMethod],kOpenCircle);
    histMerged->DrawCopy("hist p e1 same");
    legMethods->AddEntry(histMerged,sMethodLabels[iMethod].Data(),"p");
    legMethods->Draw();
    canCompMethod->cd(1);
    StyleHist(h_weight_sigma,colors[iMethod],kOpenSquare);
    h_weight_sigma->DrawCopy("hist p e1 same");
    legMethodsWeight->AddEntry(histMerged,sMethodLabels[iMethod].Data(),"p");
    legMethodsWeight->Draw();

    gSystem->mkdir(Form("%s/vn_dist_%s_%s/",sOutputPath.Data(),sSpecies.Data(),sMethods[iMethod].Data()),kTRUE);

    TLine* lineMerged = new TLine();
    lineMerged->SetLineColor(kRed);
    lineMerged->SetLineWidth(3);
    lineMerged->SetLineStyle(10);
    TLine* lineWeight = new TLine();
    lineWeight->SetLineColor(kViolet-1);
    lineWeight->SetLineWidth(3);
    for(Int_t pt(0); pt < iNumPtBins; ++pt)
    {

      gPad->SetLogy(kFALSE);
      canDist->cd();
      histVnDist[pt]->SetStats(kFALSE);
      histVnDist[pt]->Draw("hist");
      lineWeight->DrawLine(h_weight_sigma->GetBinContent(pt+1), 0.0, h_weight_sigma->GetBinContent(pt+1), histVnDist[pt]->GetMaximum());
      lineMerged->DrawLine(histMerged->GetBinContent(pt+1), 0.0, histMerged->GetBinContent(pt+1), histVnDist[pt]->GetMaximum());
      canDist->SaveAs(Form("%s/vn_dist_%s_%s/Dist_pt%d.pdf",sOutputPath.Data(),sSpecies.Data(),sMethods[iMethod].Data(),pt),"pdf");

      gPad->SetLogy(kTRUE);
      histVnDistWeighted[pt]->SetStats(kFALSE);
      histVnDistWeighted[pt]->Draw("hist");
      lineWeight->DrawLine(h_weight_sigma->GetBinContent(pt+1), 0.0, h_weight_sigma->GetBinContent(pt+1), histVnDistWeighted[pt]->GetMaximum());
      lineMerged->DrawLine(histMerged->GetBinContent(pt+1), 0.0, histMerged->GetBinContent(pt+1), histVnDistWeighted[pt]->GetMaximum());
      canDist->SaveAs(Form("%s/vn_dist_%s_%s/DistWeight_pt%d.pdf",sOutputPath.Data(),sSpecies.Data(),sMethods[iMethod].Data(),pt),"pdf");
    }


    // saving to output file
    fileOut->cd();
    h_weight_sigma->Write(Form("hWeight_%s",sMethods[iMethod].Data()));
    histMerged->Write(Form("hMerged_%s",sMethods[iMethod].Data()));


  }
  canCompMethod->SaveAs(Form("%s/%s_%s.pdf",sOutputPath.Data(),sSpecies.Data(),"CompMethods"),"pdf");
  canRatio->SaveAs(Form("%s/%s_%s.pdf",sOutputPath.Data(),sSpecies.Data(),"RATIO"),"pdf");


  return;
}
// ==================================================================================================================
TFile* OpenFile(TString sFileName, TString sMode)
{
  TFile* file = TFile::Open(sFileName.Data(),sMode.Data());
  if(!file) { printf("ERROR: Input file '%s' not found.\n",sFileName.Data()); return 0x0; }

  return file;
}
// ==================================================================================================================
TH1D* LoadHisto(TString sHistName, TFile* file)
{
  if(!file) { printf("ERROR-LoadHisto: File does not found.\n"); return 0x0; }

  TH1D* hist = (TH1D*) file->Get(sHistName.Data());
  if(!hist) { printf("ERROR-LoadHisto: Histo '%s' not found\n",sHistName.Data()); file->ls(); return 0x0; }

  return hist;
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
// ==================================================================================================================
void SetCustomPalette()
{
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,	0.6250, 0.7500, 0.8750, 1.0000};


  // Rain Bow
  // case 55:
  // Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
  // Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
  // Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};

  // Bird
  //case 57:
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};

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
