#include "TH1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
TH1D* PrepareSyst(TH1D* res, TH1D* syst);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);
void SetStyle(Bool_t graypalette=kFALSE);

// plotting order (mainly for syst)
Int_t plotOrder[] = {0,2,5,4,1,3}; // indexes in sSpecies_list

Bool_t bDrawSyst = kTRUE;

// regular order
TString sSpecies_labels[] = {" #pi^{#pm}"," K^{#pm}"," K^{0}_{S}"," p(#bar{p})"," #phi"," #Lambda(#bar{#Lambda})"};
TString sSpecies_list[] = {"Pion","Kaon","K0s","Proton","Phi","Lambda"};

Color_t colors[] = {kBlack, kRed, kOrange-1, kBlue, kGreen+2, kMagenta+1, kOrange-1};
Color_t colorsFill[] = {kBlack, kRed, kOrange-1, kBlue, kGreen+2, kMagenta+1, kOrange-1};

Int_t markers[] = {kOpenCircle, kOpenSquare, kOpenStar, kOpenCross, kOpenDiamond, kOpenCircle};
Double_t markersSize[] = {1.0,1.0,1.2,1.2,1.2,1.0,1.2};
Double_t dAlpha = 0.35;


TString sGapVal= "0.4"; TString sGapName= "gap04";


TString sCent[] = {"0-20","20-40","40-60"};

Double_t dPtBinErr = 0.0; // vertical stat err
Double_t dPtBinSyst = 0.05; // width of syst err box

std::vector<TGraphErrors*> vecGraphPoints;
std::vector<TGraphErrors*> vecGraphSyst;

Int_t iNumCent = sizeof(sCent) / sizeof(sCent[0]);
Int_t iNumSpecies = sizeof(sSpecies_list) / sizeof(sSpecies_list[0]);

void PlotPIDvnAfterSubWithSyst_0204060100()
{
  SetStyle();

  TString sInputPathTop = "./";

  TString sOutputPath = Form("%s/plots/",sInputPathTop.Data());
  gSystem->mkdir(sOutputPath.Data(),kTRUE);

  TFile* fileIn = OpenFile(Form("%s/v2-syst-final.root",sInputPathTop.Data()));
  if(!fileIn) { return; }

  // Int_t iCent(0);
  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    vecGraphPoints.clear();
    vecGraphSyst.clear();

    TLegend* leg = new TLegend(0.18,0.46,0.4,0.73);
    leg->SetFillColorAlpha(0,0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);

    // Int_t spec(0);
    for(Int_t spec(0); spec < iNumSpecies; ++spec)
    {
      TString sSpecies = sSpecies_list[spec];

      // loading central points
      TGraphErrors* graphPoints = (TGraphErrors*) fileIn->Get(Form("graphPoints_%s_cent%d",sSpecies.Data(),iCent));
      if(!graphPoints) { printf("No histPoints\n"); fileIn->ls(); return; }

      TGraphErrors* graphSyst = (TGraphErrors*) fileIn->Get(Form("graphSyst_%s_cent%d",sSpecies.Data(),iCent));
      if(!graphSyst) { printf("No histSyst \n"); fileIn->ls(); return; }


      graphPoints->SetMarkerStyle(markers[spec]);
      graphPoints->SetMarkerColor(colors[spec]);
      graphPoints->SetMarkerSize(markersSize[spec]);
      graphPoints->SetLineColor(colors[spec]);
      graphPoints->SetLineWidth(2);
      graphPoints->SetLineStyle(1);

      graphSyst->SetMarkerStyle(markers[spec]);
      graphSyst->SetMarkerColor(colors[spec]);
      graphSyst->SetMarkerSize(markersSize[spec]);
      graphSyst->SetLineColor(colors[spec]);
      graphSyst->SetLineWidth(2);
      graphSyst->SetLineStyle(1);
      graphSyst->SetFillStyle(1001);
      graphSyst->SetFillColorAlpha(colorsFill[spec],dAlpha);

      if(iCent < iNumCent-1 || spec != 4)
      {
        leg->AddEntry(graphSyst,sSpecies_labels[spec].Data(),"pf");
      }

      vecGraphPoints.push_back(graphPoints);
      vecGraphSyst.push_back(graphSyst);

      printf("vectors size %lu (points) | %lu (syst)\n", vecGraphPoints.size(), vecGraphSyst.size());
    }

    TLatex * text = new TLatex();
    TLatex * text2 = new TLatex();
    text->SetTextFont(42);
    text->SetTextSize(0.04);
    text2->SetTextFont(42);
    text2->SetTextSize(0.04);

    TCanvas* canPlot = new TCanvas("canPlot","canPlot",800,800);
    canPlot->cd();
    TH1* frame_canPlot = (TH1*) gPad->DrawFrame(0.0,0.0,7.0,0.25);
    frame_canPlot->SetTitle(Form("; #it{p}_{T} (GeV/#it{c}); v_{2}^{sub} {2, |#Delta#eta| > %s}",sGapVal.Data()));
    text2->DrawLatexNDC(0.186,0.84,Form("ALICE Preliminary"));
    text->DrawLatexNDC(0.186,0.79,Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
    text->DrawLatexNDC(0.186,0.74,Form("|#eta| < 0.8"));
    text->DrawLatexNDC(0.62,0.84,Form("%s%% (V0A)", sCent[iCent].Data()));
    // leg->SetHeader(Form("|#eta| < 0.8 (%s%% V0A)", sCent[iCent].Data()));
    leg->Draw();

    // Drawing

    if(bDrawSyst)
    {
      // systematics first
      for(Int_t i(0); i < iNumSpecies; ++i)
      {
        Int_t index = plotOrder[i];
        TGraphErrors* graphSyst = vecGraphSyst.at(index);
        canPlot->cd();
         // Phi in peripheral bin
        if(iCent == iNumCent-1 && index == 4) { continue; }
        graphSyst->Draw("p2 same");
      }
    }

    // central points
    for(Int_t i(0); i < iNumSpecies; ++i)
    {
      Int_t index = plotOrder[i];
      TGraphErrors* graphPoints = vecGraphPoints.at(index);
      canPlot->cd();
      // Phi in peripheral bin
      if(iCent == iNumCent-1 && index == 4) { continue; }
      graphPoints->Draw("p Z same");
    }

    canPlot->SaveAs(Form("%s/PIDvn_cent%d.pdf",sOutputPath.Data(),iCent),"pdf");
    canPlot->SaveAs(Form("%s/PIDvn_cent%d.eps",sOutputPath.Data(),iCent),"eps");
  }

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
//_____________________________________________________________________________
void SetStyle(Bool_t graypalette) {
  // cout << "Setting style!" << std::end;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetHistLineWidth(1);
  // gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.52,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  // gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  // gStyle->SetLegendTextSize(0.04);
}
//_____________________________________________________________________________
TH1D* PrepareSyst(TH1D* res, TH1D* syst)
{
  if(!res) {printf("nores\n"); return 0x0;}
  if(!syst) {printf("nosyst\n"); return 0x0;}

  TH1D* hResSyst = (TH1D*) res->Clone(Form("%s_sy",res->GetName()));

  for(Int_t bin(1); bin < hResSyst->GetNbinsX()+1; ++bin)
  {

    Double_t dCont = res->GetBinContent(bin);
    Double_t dSyst = syst->GetBinContent(bin);

    hResSyst->SetBinContent(bin,dCont);
    hResSyst->SetBinError(bin,dCont*dSyst);
  }

  return hResSyst;
}
//_____________________________________________________________________________
