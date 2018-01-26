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
  // raw pPb cent dep
  TFile* fileInRaw = OpenFile(sInFileRaw+"/"+sMethod+"/Processed.root"); if(!fileInRaw) { return; }
  TH1D* hRaw_Cum_Charged_cent0 = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",0),fileInRaw); if(!hRaw_Cum_Charged_cent0) { return; }
  TH1D* hRaw_Cum_Charged_cent1 = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",1),fileInRaw); if(!hRaw_Cum_Charged_cent1) { return; }
  TH1D* hRaw_Cum_Charged_cent2 = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",2),fileInRaw); if(!hRaw_Cum_Charged_cent2) { return; }
  TH1D* hRaw_Cum_Charged_cent3 = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",3),fileInRaw); if(!hRaw_Cum_Charged_cent3) { return; }

  TFile* fileInRaw_Mult = OpenFile(sInFileRaw+"/Mult.root"); if(!fileInRaw_Mult) { return; }
  TProfile* hRaw_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInRaw_Mult); if(!hRaw_Mult) { return; }

  // base pp cent int
  TFile* fileInBaseInt = OpenFile(sInFileBaseInt+"/"+sMethod+"/Processed.root"); if(!fileInBaseInt) { return; }
  TH1D* hBase_Cum_Charged_int = LoadHisto("hCum2_Charged_harm2_gap08_cent0",fileInBaseInt); if(!hBase_Cum_Charged_int) { return; }

  TFile* fileInBase_MultInt = OpenFile(sInFileBaseInt+"/Mult.root"); if(!fileInBase_MultInt) { return; }
  TProfile* hBase_MultInt = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_MultInt); if(!hBase_MultInt) { return; }
  Double_t dMult_Base_Int = hBase_MultInt->GetBinContent(1);

  // base pp cent dependetn
  TFile* fileInBase = OpenFile(sInFileBase+"/"+sMethod+"/Processed.root"); if(!fileInBase) { return; }
  TH1D* hBase_Cum_Charged_cent0 = LoadHisto("hCum2_Charged_harm2_gap08_cent0",fileInBase); if(!hBase_Cum_Charged_cent0) { return; }
  TH1D* hBase_Cum_Charged_cent1 = LoadHisto("hCum2_Charged_harm2_gap08_cent1",fileInBase); if(!hBase_Cum_Charged_cent1) { return; }
  TH1D* hBase_Cum_Charged_cent2 = LoadHisto("hCum2_Charged_harm2_gap08_cent2",fileInBase); if(!hBase_Cum_Charged_cent2) { return; }
  TH1D* hBase_Cum_Charged_cent3 = LoadHisto("hCum2_Charged_harm2_gap08_cent3",fileInBase); if(!hBase_Cum_Charged_cent3) { return; }

  TFile* fileInBase_Mult = OpenFile(sInFileBase+"/Mult.root"); if(!fileInBase_Mult) { return; }
  TProfile* hBase_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_Mult); if(!hBase_Mult) { return; }


  // settiing plotting setting
  StyleHist(hRaw_Cum_Charged_cent0, kRed, kFullCircle);
  StyleHist(hRaw_Cum_Charged_cent1, kRed, kFullCircle);
  StyleHist(hRaw_Cum_Charged_cent2, kRed, kFullCircle);
  StyleHist(hRaw_Cum_Charged_cent3, kRed, kFullCircle);

  StyleHist(hBase_Cum_Charged_int, kRed, kOpenCircle);
  StyleHist(hBase_Cum_Charged_cent0, colors[0], kOpenSquare);
  StyleHist(hBase_Cum_Charged_cent1, colors[1], kOpenSquare);
  StyleHist(hBase_Cum_Charged_cent2, colors[2], kOpenSquare);
  StyleHist(hBase_Cum_Charged_cent3, colors[3], kOpenSquare);

  // scaling base with pre-factor
  TH1D* hRaw_Cum_Charged_cent0_scaled = Scale(hRaw_Cum_Charged_cent0, hRaw_Mult->GetBinContent(0+1));
  TH1D* hRaw_Cum_Charged_cent1_scaled = Scale(hRaw_Cum_Charged_cent1, hRaw_Mult->GetBinContent(1+1));
  TH1D* hRaw_Cum_Charged_cent2_scaled = Scale(hRaw_Cum_Charged_cent2, hRaw_Mult->GetBinContent(2+1));
  TH1D* hRaw_Cum_Charged_cent3_scaled = Scale(hRaw_Cum_Charged_cent3, hRaw_Mult->GetBinContent(3+1));

  TH1D* hBase_Cum_Charged_int_scaled = Scale(hBase_Cum_Charged_int, hBase_MultInt->GetBinContent(0+1));
  TH1D* hBase_Cum_Charged_cent0_scaled = Scale(hBase_Cum_Charged_cent0, hBase_Mult->GetBinContent(0+1));
  TH1D* hBase_Cum_Charged_cent1_scaled = Scale(hBase_Cum_Charged_cent1, hBase_Mult->GetBinContent(1+1));
  TH1D* hBase_Cum_Charged_cent2_scaled = Scale(hBase_Cum_Charged_cent2, hBase_Mult->GetBinContent(2+1));
  TH1D* hBase_Cum_Charged_cent3_scaled = Scale(hBase_Cum_Charged_cent3, hBase_Mult->GetBinContent(3+1));

  TH1D* hSubPP_Cum_Charged_int = Subtract(hRaw_Cum_Charged_cent0_scaled, hBase_Cum_Charged_int_scaled);
  TH1D* hSubPP_Cum_Charged_cent0 = Subtract(hRaw_Cum_Charged_cent0_scaled, hBase_Cum_Charged_cent0_scaled);
  TH1D* hSubPP_Cum_Charged_cent1 = Subtract(hRaw_Cum_Charged_cent0_scaled, hBase_Cum_Charged_cent1_scaled);
  TH1D* hSubPP_Cum_Charged_cent2 = Subtract(hRaw_Cum_Charged_cent0_scaled, hBase_Cum_Charged_cent2_scaled);
  TH1D* hSubPP_Cum_Charged_cent3 = Subtract(hRaw_Cum_Charged_cent0_scaled, hBase_Cum_Charged_cent3_scaled);

  TLegend* leg = new TLegend(0.1,0.5,0.6,0.89);
  // leg->SetBorderWidth(0);
  // leg->SetFillColorAplha(0,0);
  leg->AddEntry(hRaw_Cum_Charged_cent0,"Unsub pPb","p");
  leg->AddEntry(hBase_Cum_Charged_int,"pp (0-100%)","p");
  leg->AddEntry(hBase_Cum_Charged_cent0,Form("pp (%s)",sCentLabel[0].Data()),"p");
  leg->AddEntry(hBase_Cum_Charged_cent1,Form("pp (%s)",sCentLabel[1].Data()),"p");
  leg->AddEntry(hBase_Cum_Charged_cent2,Form("pp (%s)",sCentLabel[2].Data()),"p");
  leg->AddEntry(hBase_Cum_Charged_cent3,Form("pp (%s)",sCentLabel[3].Data()),"p");

  TCanvas* can = new TCanvas();
  can->Divide(3,1);
  can->cd(1);
  TH1* frame = gPad->DrawFrame(0,0,10,0.2);
  frame->SetTitle("raw <<2'>>");
  hRaw_Cum_Charged_cent0->Draw("same");
  hBase_Cum_Charged_int->Draw("same");
  hBase_Cum_Charged_cent0->Draw("same");
  hBase_Cum_Charged_cent1->Draw("same");
  hBase_Cum_Charged_cent2->Draw("same");
  hBase_Cum_Charged_cent3->Draw("same");
  leg->Draw();
  can->cd(2);
  TH1* frame2 = gPad->DrawFrame(0,-0.03,10,1.0);
  frame2->SetTitle("<M> * <<2'>>");
  hRaw_Cum_Charged_cent0_scaled->Draw("same");
  hBase_Cum_Charged_int_scaled->Draw("same");
  hBase_Cum_Charged_cent0_scaled->Draw("same");
  hBase_Cum_Charged_cent1_scaled->Draw("same");
  hBase_Cum_Charged_cent2_scaled->Draw("same");
  hBase_Cum_Charged_cent3_scaled->Draw("same");
  can->cd(3);
  TH1* frame3 = gPad->DrawFrame(0,-0.03,10,1.0);
  frame3->SetTitle("<M>^{pPb}<<2'>>^{pPb} - <M>^{pp}<<2'>>^{pp} ");
  hSubPP_Cum_Charged_int->Draw("same");
  hSubPP_Cum_Charged_cent0->Draw("same");
  hSubPP_Cum_Charged_cent1->Draw("same");
  hSubPP_Cum_Charged_cent2->Draw("same");
  hSubPP_Cum_Charged_cent3->Draw("same");

  can->SaveAs(Form("%s/Subt_pp-pbp_cent%d.pdf",sOutFolder.Data(),iCent),"pdf");

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
