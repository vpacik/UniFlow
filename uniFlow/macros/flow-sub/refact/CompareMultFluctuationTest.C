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

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);
void SetCustomPalette();

// Color_t colors[] = {kGreen+2, kSpring+2, kCyan+2, kBlue+2, kMagenta+2, kRed-2, kRed+2, kOrange+2, kYellow+2, kOrange-4};
Color_t colors[] = {kGreen+2, kGreen, kSpring+2, kSpring, kCyan+2, kCyan, kBlue+2, kBlue, kMagenta+2, kRed-2, kRed+2, kOrange+2, kYellow+2, kOrange-4,};


void CompareMultFluctuationTest()
{
  const Int_t iNumMethods = 4;
  TString sMethods[iNumMethods] = {"GF_eventweighted","GF_noneventweighted","SP_nonscaled_noneventweighted","SP_scaled_noneventweighted"};

  TString sSpecies = "Charged";
  TString sGap = "gap00";

  Int_t iMultMin = 0;
  Int_t iMultMax = 30;
  Int_t iMultBinWidth = 5;

  // file with correlations in 0-20% (V0A) vs. N_RFP
  TString sInputPathUnit = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-NchRFP-gap00/output_vn_5/";
  // file with correlations in 0-20% (V0A)  (i.e vs centrality percentile)
  TString sInputPathMerged = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-NchRFP-gap00/output_vn_merged/";

  TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/CompareMultFluctuationTest/pPb_gap00/5/";
  // ##########################################################################################################
  gSystem->mkdir(sOutputPath.Data(),kTRUE);

  // gStyle->SetPalette(1);
  SetCustomPalette();

  Double_t dYmin = -100.0;
  Double_t dYmax = 100.0;
  Double_t dXmin = 0.0;
  Double_t dXmax = 10.0;

  TCanvas* canAll = new TCanvas("canAll","canAll", 1200,1200);
  canAll->Divide(2,2);
  canAll->cd(1);
  TH1* frame_canAll_1 = (TH1*) gPad->DrawFrame(dXmin,dYmin,dXmax,dYmax);
  canAll->cd(2);
  TH1* frame_canAll_2 = (TH1*) gPad->DrawFrame(dXmin,dYmin,dXmax,dYmax);
  canAll->cd(3);
  TH1* frame_canAll_3 = (TH1*) gPad->DrawFrame(dXmin,dYmin,dXmax,dYmax);
  canAll->cd(4);
  TH1* frame_canAll_4 = (TH1*) gPad->DrawFrame(dXmin,dYmin,dXmax,dYmax);

  TCanvas* canAllRatio = new TCanvas("canAllRatio","canAllRatio",800,800);
  canAllRatio->cd(1);
  TH1* frame_canAllRatio_1 = (TH1*) gPad->DrawFrame(dXmin,0.0,dXmax,2.0);

    // Int_t iMult = 0;
  for(Int_t iMethod(0); iMethod < iNumMethods; ++iMethod)
  {
    TLegend* leg = new TLegend(0.0,0.0,0.32,1.0);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0.,0.);
    leg->SetHeader("N_{RFP} bin");
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);

    TLegend* leg2 = new TLegend(0.33,0.0,0.65,1.0);
    leg2->SetBorderSize(0);
    leg2->SetFillColorAlpha(0.,0.);
    leg2->SetHeader("N_{RFP} bin");
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);

    TLegend* leg3 = new TLegend(0.66,0.0,1.0,1.0);
    leg3->SetBorderSize(0);
    leg3->SetFillColorAlpha(0.,0.);
    leg3->SetHeader("N_{RFP} bin");
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.03);

    TCanvas* canDiff = new TCanvas(Form("canDiff"),Form("canDiff"),2400,800);
    canDiff->Divide(3,1);
    canDiff->cd(1);
    TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(0.0,-1,10.0,1);
    frame_canDiff_1->SetTitle(Form("%s v_{2}{2,%s} (0-20%% V0A) ; p_{T} (GeV/c); v_{2}{2}",sSpecies.Data(),sGap.Data()));

    canDiff->cd(2);
    TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(0.0,0.5,10.0,1.5);
    frame_canDiff_2->SetTitle(Form("; p_{T} (GeV/c); ratio"));

    TFile* fileUnit = OpenFile(sInputPathUnit+sMethods[iMethod]+"/Processed.root","READ"); if(!fileUnit) return;
    TFile* fileMerged = OpenFile(sInputPathMerged+sMethods[iMethod]+"/Processed.root","READ"); if(!fileMerged) return;

    TH1D* histMerged = LoadHisto(Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),0),fileMerged); if(!histMerged) return;
    StyleHist(histMerged,kRed,kFullCircle);
    leg->AddEntry(histMerged,"integrated","p");

    Int_t nPnt  = iMultMax-iMultMin;
    Int_t nnCol = gStyle->GetNumberOfColors();


    // const Int_t iBins = histMerged->GetXaxis()->GetNbins();
    // const Double_t* bins = histMerged->GetXaxis()->GetXbins()->GetArray();

    TProfile* p_weight_average = new TProfile("p_weight_average","p_weight_average",histMerged->GetXaxis()->GetNbins(),histMerged->GetXaxis()->GetXbins()->GetArray());
    if(!p_weight_average) { printf("Profile does not exists!\n"); return; }
    p_weight_average->Sumw2();
    StyleHist(p_weight_average,kMagenta+2,kOpenSquare);
    leg->AddEntry(p_weight_average,"simple average","p");

    TProfile* p_weight_sigma = new TProfile("p_weight_sigma","p_weight_sigma",histMerged->GetXaxis()->GetNbins(),histMerged->GetXaxis()->GetXbins()->GetArray());
    if(!p_weight_sigma) { printf("Profile does not exists!\n"); return; }
    p_weight_sigma->Sumw2();
    StyleHist(p_weight_sigma,kBlack,kOpenSquare);
    leg->AddEntry(p_weight_sigma,"1/#sigma^{2}","p");

    for(Int_t iMult(iMultMin); iMult < iMultMax; ++iMult)
    {
      TString sHistoName = Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),iMult);

      TH1D* histUnit = LoadHisto(sHistoName,fileUnit); if(!histUnit) return;


      TLegend* curLeg = 0x0;
      if(iMult < 0.33*nPnt) { curLeg = leg; }
      else if(iMult < 0.66*nPnt) { curLeg = leg2; }
      else { curLeg = leg3; }
      curLeg->AddEntry(histUnit,Form("%d",iMult*iMultBinWidth),"p");

      Int_t idx = iMult * Float_t(nnCol-1) / (nPnt-1);
      StyleHist(histUnit, gStyle->GetColorPalette(idx), kOpenCircle);

      canDiff->cd(1);
      histUnit->Draw("hist p x0 same");
      canAll->cd(iMethod+1);
      histUnit->Draw("hist p e1 same");


      for(Int_t pt(0); pt < histUnit->GetNbinsX()+1; ++pt)
      {
        Double_t bin = histUnit->GetBinCenter(pt);
        Double_t dValue = histUnit->GetBinContent(pt);
        Double_t dError = histUnit->GetBinError(pt);

        if(dError <= 0.0) continue;
        if(dValue == 9.0 && dError == 9.0) continue;

        Double_t dWeight = TMath::Power(dError,-2.0);
        p_weight_average->Fill(bin,dValue,1.0);
        p_weight_sigma->Fill(bin,dValue,dWeight);
        printf("Val %f | weight %f\n",dValue,dWeight);
      }

    }
    TH1D* h_weight_average = (TH1D*) p_weight_average->ProjectionX();
    TH1D* h_weight_sigma = (TH1D*) p_weight_sigma->ProjectionX();
    StyleHist(h_weight_average,kMagenta+2,kOpenSquare);
    StyleHist(h_weight_sigma,kBlack,kOpenSquare);

    canDiff->cd(1);
    histMerged->Draw("same");
    h_weight_average->Draw("same");
    h_weight_sigma->Draw("same");



    canDiff->cd(2);
    TH1D* hRatioSigma = DivideHistos(h_weight_sigma,histMerged);
    TH1D* hRatioAverage = DivideHistos(h_weight_average,histMerged);
    hRatioSigma->Draw("hist p e1 same");
    hRatioAverage->Draw("hist p e1 same");

    canDiff->cd(3);
    leg->Draw();
    leg2->Draw();
    leg3->Draw();

    canDiff->SaveAs(Form("%s/%s_%s.pdf",sOutputPath.Data(),sSpecies.Data(),sMethods[iMethod].Data()),"pdf");

    // TCanvas* canProf = new TCanvas();
    // canProf->cd();
    // prof->Draw();
    canAll->cd(iMethod+1);
    histMerged->Draw("same");
    p_weight_average->Draw("same");
    p_weight_sigma->Draw("same");

    // canAll->cd(4);
    // leg->Draw();

    canAllRatio->cd();
    hRatioSigma->Draw("same");
    hRatioAverage->Draw("same");
  }
  canAll->SaveAs(Form("%s/%s_%s.pdf",sOutputPath.Data(),sSpecies.Data(),"compall"),"pdf");
  canAllRatio->SaveAs(Form("%s/%s_%s.pdf",sOutputPath.Data(),sSpecies.Data(),"compallratio"),"pdf");

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
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,
			0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  Int_t pal = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
  const Int_t nCol = 255;
  Int_t colors[nCol];
  for (int i=0; i<nCol; i++) colors[i] =pal+i;


  gStyle->SetPalette(nCol,colors);
}
