#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TColor.h"
#include "TStyle.h"
#include "TMath.h"

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
void SetCustomPalette();
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);

// colors for centrality
Color_t colors[] = {kRed, kGreen+2, kBlue, kMagenta+1};


void Plot_etaGaps()
{
  TString sCentLabel[] = {"0-20%", "20-40%", "40-60%", "60-100%"};
  TString sSpecies = "Charged";
  TString sMethod = "GF_eventweighted";
  TString sOutputTag = "output_vn";
  TString sFileName = "Subt_results.root";


  // TString sListName = "list_SubtPP_vn_ppcent"; TString sHistName = "hSubPP_vn_cent"; TString sListHeader = "pPb - pp (cent)"; Int_t iNumCent = 4;
  // TString sListName = "list_SubtPP_vn_ppint"; TString sHistName = "hSubPP_vn_int_cent"; TString sListHeader = "pPb - pp (MB)"; Int_t iNumCent = 4;
  TString sListName = "list_SubtPPb_vn"; TString sHistName = "hSubPPb_vn_cent"; TString sListHeader = "pPb - pPb (60-100%)"; Int_t iNumCent = 3;


  TString sOutputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etagap_3/";

  const Int_t iNumFiles = 7;
  TString sFiles[iNumFiles] = {
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap00/" + sOutputTag + "/" + sMethod + "/pPb00_pp00_subt/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap02/" + sOutputTag + "/" + sMethod + "/pPb02_pp02_subt/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gaps-04-06-10-12/" + sOutputTag + "/gap04/" + sMethod + "/pPb_pp_subt_gap04/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gaps-04-06-10-12/" + sOutputTag + "/gap06/" + sMethod + "/pPb_pp_subt_gap06/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap08/" + sOutputTag + "/" + sMethod + "/pPb08_pp08_subt/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gaps-04-06-10-12/" + sOutputTag + "/gap10/" + sMethod + "/pPb_pp_subt_gap10/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gaps-04-06-10-12/" + sOutputTag + "/gap12/" + sMethod + "/pPb_pp_subt_gap12/" + sSpecies + "/" + sFileName
  };

  TString sFileLabel[iNumFiles] = {
    "|#Delta#eta| > 0.0",
    "|#Delta#eta| > 0.2",
    "|#Delta#eta| > 0.4",
    "|#Delta#eta| > 0.6",
    "|#Delta#eta| > 0.8",
    "|#Delta#eta| > 1.0",
    "|#Delta#eta| > 1.2"
  };

  gSystem->mkdir(sOutputFolder.Data(),1);


  SetCustomPalette();
  Int_t nPnt  = iNumFiles+1;
  Int_t nnCol = gStyle->GetNumberOfColors();

  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TCanvas* can = new TCanvas("can","can",400,400);
    can->cd();
    TH1* frame_can = (TH1*) gPad->DrawFrame(0.,0.,10.,0.3);
    frame_can->SetTitle(Form("v_{2}{2}^{sub} (%s V0A); p_{T} (GeV/c)",sCentLabel[cent].Data()));

    TLegend* leg = new TLegend(0.12,0.65,0.6,0.88);
    leg->SetHeader(sListHeader.Data());
    leg->SetBorderSize(0.);
    leg->SetFillColor(0);



    for(Int_t file(0); file < iNumFiles; ++file)
    {

      Int_t idx = file * Float_t(nnCol-1) / (nPnt-1);
      Color_t color = gStyle->GetColorPalette(idx);

      TFile* fileIn = TFile::Open(sFiles[file].Data(),"READ");
      TList* list = (TList*) fileIn->Get(sListName.Data());
      if(!list) { printf("List '%s' not found\n",sListName.Data()); fileIn->ls(); return; }
      TH1D* hist = (TH1D*) list->FindObject(Form("%s%d",sHistName.Data(),cent));
      if(!hist) { printf("Hist '%s%d' not found\n",sHistName.Data(),cent); list->ls(); return; }
      StyleHist(hist,color,kOpenCircle);
      hist->Draw("same");
      leg->AddEntry(hist,sFileLabel[file].Data(),"p");
    }

    leg->Draw();
    can->SaveAs(Form("%s%s_cent%d.pdf",sOutputFolder.Data(),sListName.Data(),cent),"pdf");
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
