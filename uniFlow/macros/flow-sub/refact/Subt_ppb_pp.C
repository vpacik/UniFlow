#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

TFile* OpenFile(TString sFileName);
TH1D* LoadHisto(TString sHistName, TFile* file);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);

TH1D* Scale(TH1D* base, Double_t factor);
TH1D* Subtract(TH1D* raw, TH1D* base, Double_t factor = 1.0);

// colors for centrality
Color_t colors[] = {kGreen+2, kBlue, kBlack, kMagenta+1};


void Subt_ppb_pp()
{
  TString sMethod = "GF_eventweighted";
  TString sOutputTag = "output_vn";
  TString sOutputTagInt = sOutputTag + "_int";

  Int_t iCent = 1;
  TString sInFileRaw = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap08/" + sOutputTag;
  TString sInFileBase = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run3-2-gap08/" + sOutputTag;
  TString sInFileBaseInt = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run3-2-gap08/" + sOutputTagInt;
  TString sOutFolder = sInFileRaw+"/"+sMethod+"/plots";

  const Int_t iNumCent = 4;
  TString sCentLabel[iNumCent] = {"0-20%", "20-40%", "40-60%", "60-100%"};

  // ==================================================================================================================
  // Making output folder
  gSystem->mkdir(sOutFolder,kTRUE);

  // LOADING INPUT

  TFile* fileInRaw = OpenFile(sInFileRaw+"/"+sMethod+"/Processed.root"); if(!fileInRaw) { return; }
  TFile* fileInBase = OpenFile(sInFileBase+"/"+sMethod+"/Processed.root"); if(!fileInBase) { return; }
  TFile* fileInBaseInt = OpenFile(sInFileBaseInt+"/"+sMethod+"/Processed.root"); if(!fileInBaseInt) { return; }

  TH1D* hBase_Cum_Charged_int = LoadHisto("hCum2_Charged_harm2_gap08_cent0",fileInBaseInt); if(!hBase_Cum_Charged_int) { return; }
  StyleHist(hBase_Cum_Charged_int, kRed, kOpenCircle);

  TList* list_Raw_Cum_Charged = new TList();
  TList* list_Base_Cum_Charged = new TList();

  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TH1D* temp = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",cent),fileInRaw); if(!temp) { return; }
    StyleHist(temp, colors[cent], kFullCircle);
    list_Raw_Cum_Charged->Add(temp);

    temp = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",cent),fileInBase); if(!temp) { return; }
    StyleHist(temp, colors[cent], kOpenSquare);
    list_Base_Cum_Charged->Add(temp);
  }

  // multiplicities
  TFile* fileInRaw_Mult = OpenFile(sInFileRaw+"/Mult.root"); if(!fileInRaw_Mult) { return; }
  TProfile* hRaw_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInRaw_Mult); if(!hRaw_Mult) { return; }

  TFile* fileInBase_Mult = OpenFile(sInFileBase+"/Mult.root"); if(!fileInBase_Mult) { return; }
  TProfile* hBase_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_Mult); if(!hBase_Mult) { return; }

  TFile* fileInBase_MultInt = OpenFile(sInFileBaseInt+"/Mult.root"); if(!fileInBase_MultInt) { return; }
  TProfile* hBase_MultInt = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_MultInt); if(!hBase_MultInt) { return; }
  Double_t dMult_Base_Int = hBase_MultInt->GetBinContent(1);

  // scaling base with pre-factor
  TH1D* hBase_Cum_Charged_int_scaled = Scale(hBase_Cum_Charged_int, hBase_MultInt->GetBinContent(1));

  TList* list_Raw_Cum_Charged_scaled = new TList();
  TList* list_Base_Cum_Charged_scaled = new TList();

  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TH1D* raw_scaled = Scale((TH1D*) list_Raw_Cum_Charged->At(cent), hRaw_Mult->GetBinContent(cent+1)); if(!raw_scaled) { return; }
    list_Raw_Cum_Charged_scaled->Add(raw_scaled);

    TH1D* base_scaled = Scale((TH1D*) list_Base_Cum_Charged->At(cent), hBase_Mult->GetBinContent(cent+1)); if(!base_scaled) { return; }
    list_Base_Cum_Charged_scaled->Add(base_scaled);
  }

  // subtracting (scaled) histos
  TList* list_SubtPP_Cum_Charged = new TList();

  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* temp_raw_cent = (TH1D*) list_Raw_Cum_Charged_scaled->At(centRaw);
    TH1D* hSubPP_Cum_Charged_int = Subtract(temp_raw_cent, hBase_Cum_Charged_int_scaled);

    TLegend* leg = new TLegend(0.1,0.5,0.6,0.89);
    leg->AddEntry(temp_raw_cent,Form("Unsub pPb (%s)",sCentLabel[centRaw].Data()),"p");
    leg->AddEntry(hBase_Cum_Charged_int,"pp (0-100%)","p");

    TCanvas* can = new TCanvas("can","can",1200,400);
    can->Divide(3,1);
    can->cd(1);
    TH1* frame = gPad->DrawFrame(0,0,10,0.2);
    frame->SetTitle("raw <<2'>>");
    ((TH1D*) list_Raw_Cum_Charged->At(centRaw))->Draw("same");
    hBase_Cum_Charged_int->Draw("same");

    can->cd(2);
    TH1* frame2 = gPad->DrawFrame(0,-0.03,10,1.0);
    frame2->SetTitle("<M> * <<2'>>");
    temp_raw_cent->Draw("same");
    hBase_Cum_Charged_int_scaled->Draw("same");

    can->cd(3);
    TH1* frame3 = gPad->DrawFrame(0,-0.3,10,0.5);
    frame3->SetTitle("<M>^{pPb}<<2'>>^{pPb} - <M>^{pp}<<2'>>^{pp} ");
    hSubPP_Cum_Charged_int->Draw("same");

    for(Int_t cent(0); cent < iNumCent; ++cent)
    {
      TH1D* hSubPP_Cum_Charged = Subtract(temp_raw_cent, (TH1D*) list_Base_Cum_Charged_scaled->At(cent));
      list_SubtPP_Cum_Charged->Add(hSubPP_Cum_Charged);

      leg->AddEntry(hSubPP_Cum_Charged,Form("pp (%s)",sCentLabel[cent].Data()),"p");

      can->cd(1);
      ((TH1D*) list_Base_Cum_Charged->At(cent))->Draw("same");

      can->cd(2);
      ((TH1D*) list_Base_Cum_Charged_scaled->At(cent))->Draw("same");

      can->cd(3);
      hSubPP_Cum_Charged->Draw("same");
    }

    can->cd(1);
    leg->Draw();
    can->SaveAs(Form("%s/Subt_pp-pbp_cent%d.pdf",sOutFolder.Data(),centRaw),"pdf");

  }

  return;
}
// ==================================================================================================================
TFile* OpenFile(TString sFileName)
{
  TFile* file = TFile::Open(sFileName.Data(),"READ");
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
TH1D* Scale(TH1D* base, Double_t factor)
{
  if(!base) { printf("ERROR-Scale: Hist 'base' does not found.\n"); return 0x0; }

  TH1D* scale = (TH1D*) base->Clone(Form("%s_scaled",base->GetName()));
  if(!scale) { printf("ERROR-scale: Hist 'scale' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < base->GetNbinsX()+1; ++bin)
  {
    scale->SetBinContent(bin, factor * base->GetBinContent(bin));
    scale->SetBinError(bin, TMath::Abs(factor * base->GetBinError(bin)));
  }

  return scale;
}

// ==================================================================================================================
TH1D* Subtract(TH1D* raw, TH1D* base, Double_t factor)
{
  if(!raw) { printf("ERROR-Subtract: Hist 'raw' does not found.\n"); return 0x0; }
  if(!base) { printf("ERROR-Subtract: Hist 'base' does not found.\n"); return 0x0; }

  TH1D* sub = (TH1D*) base->Clone(Form("%s_sub",raw->GetName()));
  if(!sub) { printf("ERROR-Subtract: Hist 'sub' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < sub->GetNbinsX()+1; ++bin)
  {
    Double_t con_raw = raw->GetBinContent(bin);
    Double_t err_raw = raw->GetBinError(bin);
    Double_t con_base = base->GetBinContent(bin);
    Double_t err_base = base->GetBinError(bin);

    sub->SetBinContent(bin, con_raw - factor * con_base);
    sub->GetBinError(bin, 0.01*con_raw);

  }
  return sub;
}
// ==================================================================================================================
