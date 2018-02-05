#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
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

  Int_t iNumCent = 4;
  TString sListName = "list_SubtPP_vn_ppcent";
  TString sHistName = "hSubPP_vn_cent";
  TString sListHeader = "pPb - pp (cent)";


  TString sOutputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etagap/";

  const Int_t iNumFiles = 3;
  TString sFiles[iNumFiles] = {
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap00/" + sOutputTag + "/" + sMethod + "/pPb00_pp00_subt/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap02/" + sOutputTag + "/" + sMethod + "/pPb02_pp02_subt/" + sSpecies + "/" + sFileName,
    "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap08/" + sOutputTag + "/" + sMethod + "/pPb08_pp08_subt/" + sSpecies + "/" + sFileName
  };

  TString sFileLabel[iNumFiles] = {
    "|#Delta#eta| > 0.0",
    "|#Delta#eta| > 0.2",
    "|#Delta#eta| > 0.8"
  };

  gSystem->mkdir(sOutputFolder.Data(),1);

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
      TFile* fileIn = TFile::Open(sFiles[file].Data(),"READ");
      TList* list = (TList*) fileIn->Get(sListName.Data());
      if(!list) { fileIn->ls(); return; }
      TH1D* hist = (TH1D*) list->FindObject(Form("%s%d",sHistName.Data(),cent));
      if(!hist) { list->ls(); return; }
      StyleHist(hist,colors[file],kOpenCircle);
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
