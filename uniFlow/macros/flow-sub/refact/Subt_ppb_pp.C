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
  TString sOutFolder = sInFileRaw+"/"+sMethod+"/pPb-pp";
  TString sOutFile = sOutFolder+"/Subt_results.root";

  const Int_t iNumCent = 4;
  TString sCentLabel[iNumCent] = {"0-20%", "20-40%", "40-60%", "60-100%"};

  // ==================================================================================================================
  // Making output folder
  gSystem->mkdir(sOutFolder,kTRUE);

  // LOADING INPUT

  TFile* fileInRaw = OpenFile(sInFileRaw+"/"+sMethod+"/Processed.root"); if(!fileInRaw) { return; }
  TFile* fileInBase = OpenFile(sInFileBase+"/"+sMethod+"/Processed.root"); if(!fileInBase) { return; }
  TFile* fileInBaseInt = OpenFile(sInFileBaseInt+"/"+sMethod+"/Processed.root"); if(!fileInBaseInt) { return; }
  TFile* fileOut = OpenFile(sOutFile,"RECREATE"); if(!fileOut) { return; }

  TH1D* hBase_Cum_Charged_int = LoadHisto("hCum2_Charged_harm2_gap08_cent0",fileInBaseInt); if(!hBase_Cum_Charged_int) { return; }
  StyleHist(hBase_Cum_Charged_int, kBlue, kFullCircle);

  TList* list_Raw_Cum_Charged = new TList();
  TList* list_Base_Cum_Charged = new TList();

  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TH1D* temp = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",cent),fileInRaw); if(!temp) { return; }
    StyleHist(temp, kRed, kOpenSquare);
    list_Raw_Cum_Charged->Add(temp);

    temp = LoadHisto(Form("hCum2_Charged_harm2_gap08_cent%d",cent),fileInBase); if(!temp) { return; }
    StyleHist(temp, kGreen+2, kFullCircle);
    list_Base_Cum_Charged->Add(temp);
  }

  TH1D* hRaw_Cum_Refs = LoadHisto("hCum2_Refs_harm2_gap08", fileInRaw); if(!hRaw_Cum_Refs) { return; }
  StyleHist(hRaw_Cum_Refs, kRed, kOpenSquare);

  TH1D* hBase_Cum_Refs = LoadHisto("hCum2_Refs_harm2_gap08",fileInBase); if(!hBase_Cum_Refs) { return; }
  StyleHist(hBase_Cum_Refs, kGreen+2, kFullCircle);

  TH1D* hBase_Cum_Refs_int = LoadHisto("hCum2_Refs_harm2_gap08",fileInBaseInt); if(!hBase_Cum_Refs_int) { return; }
  StyleHist(hBase_Cum_Refs_int, kBlue, kFullCircle);

  // multiplicities
  TFile* fileInRaw_Mult = OpenFile(sInFileRaw+"/Mult.root"); if(!fileInRaw_Mult) { return; }
  TProfile* hRaw_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInRaw_Mult); if(!hRaw_Mult) { return; }

  TFile* fileInBase_Mult = OpenFile(sInFileBase+"/Mult.root"); if(!fileInBase_Mult) { return; }
  TProfile* hBase_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_Mult); if(!hBase_Mult) { return; }

  TFile* fileInBase_MultInt = OpenFile(sInFileBaseInt+"/Mult.root"); if(!fileInBase_MultInt) { return; }
  TProfile* hBase_MultInt = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_MultInt); if(!hBase_MultInt) { return; }
  Double_t dMult_Base_Int = hBase_MultInt->GetBinContent(1);

  // SUBTRACTION

  // Working on the c_n{2} Subtraction
  // cn{2}^sub = <M>^raw * cn{2}^raw - <M>^base * cn{2}^base

  // scaling REFS by <M>^2
  TH1D* hRaw_Cum_Refs_scaled = (TH1D*) hRaw_Cum_Refs->Clone(Form("%s_scaled",hRaw_Cum_Refs->GetName())); if(!hRaw_Cum_Refs_scaled) { return; }
  for(Int_t bin(1); bin < hRaw_Cum_Refs_scaled->GetNbinsX()+1; ++bin)
  {
    hRaw_Cum_Refs_scaled->SetBinContent(bin, hRaw_Cum_Refs->GetBinContent(bin) * TMath::Power(hRaw_Mult->GetBinContent(bin),2.0) );
  }

  TH1D* hBase_Cum_Refs_scaled = (TH1D*) hBase_Cum_Refs->Clone(Form("%s_scaled",hBase_Cum_Refs->GetName())); if(!hBase_Cum_Refs_scaled) { return; }
  for(Int_t bin(1); bin < hBase_Cum_Refs_scaled->GetNbinsX()+1; ++bin)
  {
    hBase_Cum_Refs_scaled->SetBinContent(bin, hBase_Cum_Refs->GetBinContent(bin) * TMath::Power(hBase_Mult->GetBinContent(bin),2.0) );
  }

  TH1D* hBase_Cum_Refs_int_scaled = (TH1D*) hBase_Cum_Refs_int->Clone(Form("%s_scaled",hBase_Cum_Refs_int->GetName())); if(!hBase_Cum_Refs_int_scaled) { return; }
  for(Int_t bin(1); bin < hBase_Cum_Refs_int_scaled->GetNbinsX()+1; ++bin)
  {
    hBase_Cum_Refs_int_scaled->SetBinContent(bin, hBase_Cum_Refs_int->GetBinContent(bin) * TMath::Power(hBase_MultInt->GetBinContent(bin),2.0) );
  }

  TH1D* hSub_Cum_Refs = Subtract(hRaw_Cum_Refs_scaled, hBase_Cum_Refs_scaled,1.0); if(!hSub_Cum_Refs) { return; }
  TH1D* hSub_Cum_Refs_int = Subtract(hRaw_Cum_Refs_scaled, hBase_Cum_Refs_int_scaled,1.0); if(!hSub_Cum_Refs_int) { return; }

  fileOut->cd();
  hRaw_Cum_Refs->Write("hRaw_Cum_Refs");
  hBase_Cum_Refs->Write("hBase_Cum_Refs");
  hBase_Cum_Refs_int->Write("hBase_Cum_Refs_int");
  hRaw_Cum_Refs_scaled->Write("hRaw_Cum_Refs_scaled");
  hBase_Cum_Refs_scaled->Write("hBase_Cum_Refs_scaled");
  hBase_Cum_Refs_int_scaled->Write("hBase_Cum_Refs_int_scaled");
  hSub_Cum_Refs->Write("hSub_Cum_Refs");
  hSub_Cum_Refs_int->Write("hSub_Cum_Refs_int");

  // dn{2}^sub = <M>dn{2}^raw - <M>dn{2}^base
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

  fileOut->cd();
  list_Raw_Cum_Charged->Write("list_Raw_Cum_Charged",TObject::kSingleKey);
  list_Base_Cum_Charged->Write("list_Base_Cum_Charged",TObject::kSingleKey);
  hBase_Cum_Charged_int->Write();
  list_Raw_Cum_Charged_scaled->Write("list_Raw_Cum_Charged_scaled",TObject::kSingleKey);
  list_Base_Cum_Charged_scaled->Write("list_Base_Cum_Charged_scaled",TObject::kSingleKey);
  hBase_Cum_Charged_int_scaled->Write();


  // subtracting (scaled) histos
  TList* list_SubtPP_Cum_Charged = new TList();
  TList* list_SubtPP_Cum_Charged_int = new TList();
  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* temp_raw_cent = (TH1D*) list_Raw_Cum_Charged_scaled->At(centRaw);
    // subtracting (cent) pp
    TH1D* hSubPP_Cum_Charged = Subtract(temp_raw_cent, (TH1D*) list_Base_Cum_Charged_scaled->At(centRaw));
    list_SubtPP_Cum_Charged->Add(hSubPP_Cum_Charged);

    // subtracting (int) pp
    TH1D* hSubPP_Cum_Charged_int = Subtract(temp_raw_cent, hBase_Cum_Charged_int_scaled);
    list_SubtPP_Cum_Charged_int->Add(hSubPP_Cum_Charged_int);

    TLegend* leg = new TLegend(0.1,0.5,0.6,0.89);
    leg->SetBorderSize(0.);
    leg->SetFillColor(0);
    leg->AddEntry(temp_raw_cent,Form("Unsub pPb (%s)",sCentLabel[centRaw].Data()),"p");
    leg->AddEntry(hBase_Cum_Charged_int,"pp (0-100%)","p");
    leg->AddEntry(hSubPP_Cum_Charged,Form("pp (%s)",sCentLabel[centRaw].Data()),"p");

    TCanvas* can = new TCanvas("can","can",1200,400);
    can->Divide(3,1);
    can->cd(1);
    TH1* frame = gPad->DrawFrame(0,0,10,0.2);
    frame->SetTitle("raw <<2'>>; p_{T} (GeV/c)");
    ((TH1D*) list_Raw_Cum_Charged->At(centRaw))->Draw("same");
    hBase_Cum_Charged_int->Draw("same");
    ((TH1D*) list_Base_Cum_Charged->At(centRaw))->Draw("same");
    leg->Draw();

    can->cd(2);
    TH1* frame2 = gPad->DrawFrame(0,-0.03,10,1.0);
    frame2->SetTitle("<M> * <<2'>>; p_{T} (GeV/c)");
    temp_raw_cent->Draw("same");
    hBase_Cum_Charged_int_scaled->Draw("same");
    ((TH1D*) list_Base_Cum_Charged_scaled->At(centRaw))->Draw("same");

    can->cd(3);
    TH1* frame3 = gPad->DrawFrame(0,-0.3,10,0.5);
    frame3->SetTitle("<M>^{pPb}<<2'>>^{pPb} - <M>^{pp}<<2'>>^{pp}; p_{T} (GeV/c)");
    hSubPP_Cum_Charged_int->Draw("same");
    hSubPP_Cum_Charged->Draw("same");

    can->SaveAs(Form("%s/Subt_pp-pbp_cent%d.pdf",sOutFolder.Data(),centRaw),"pdf");
  }

  // Saving to output ROOT file
  fileOut->cd();
  list_SubtPP_Cum_Charged->Write("list_SubtPP_Charged_cum_ppcent",TObject::kSingleKey);
  list_SubtPP_Cum_Charged_int->Write("list_SubtPP_Charged_cum_ppint",TObject::kSingleKey);

  // Making vn^sub out of dn^sub / sqrt(cn^sub)
  // making final vn{2}^{sub} = dn{2}^{sub} / sqrt(cn{2}^{sub})

  TList* list_SubtPP_vn_Charged = new TList();
  TList* list_SubtPP_vn_Charged_int = new TList();
  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* hSub_dn = (TH1D*) list_SubtPP_Cum_Charged->At(centRaw); if(!hSub_dn) { return; }

    // dividing by cn{2}^sub (cent)
    TH1D* hSubPP_vn = (TH1D*) hSub_dn->Clone(Form("%s_vn", hSub_dn->GetName())); if(!hSubPP_vn) { return; }
    StyleHist(hSubPP_vn, kGreen+2, kFullCircle);
    hSubPP_vn->Scale(1.0/TMath::Sqrt(hSub_Cum_Refs->GetBinContent(centRaw+1)));
    list_SubtPP_vn_Charged->Add(hSubPP_vn);

    // dividing by cn{2}^sub (int)
    TH1D* hSubPP_vn_int = (TH1D*) hSub_dn->Clone(Form("%s_vn_int", hSub_dn->GetName())); if(!hSubPP_vn_int) { return; }
    StyleHist(hSubPP_vn_int, kBlue, kFullCircle);
    hSubPP_vn_int->Scale(1.0/TMath::Sqrt(hSub_Cum_Refs_int->GetBinContent(1)));
    list_SubtPP_vn_Charged_int->Add(hSubPP_vn_int);

    TCanvas* canVn = new TCanvas("canVn","canVn",400,400);
    // canVn->Divide(,1);
    canVn->cd(1);
    TH1* frame_vn = (TH1*) gPad->DrawFrame(0.,0.,10.,0.5);
    frame_vn->SetTitle("sub vn{2}; p_{T} (GeV/c)");
    hSubPP_vn->Draw("same");
    hSubPP_vn_int->Draw("same");
    // canVn->cd(2);
    // canVn->cd(3);

    canVn->SaveAs(Form("%s/Subt_ppb_pp_vn_cent%d.pdf",sOutFolder.Data(),centRaw),"pdf");

  }

  fileOut->cd();
  list_SubtPP_vn_Charged->Write("list_SubtPP_Charged_vn_ppcent",TObject::kSingleKey);
  list_SubtPP_vn_Charged_int->Write("list_SubtPP_Charged_vn_ppint",TObject::kSingleKey);






  // UNIVERSALL PLOTTING
  // RFPs
  TLegend* leg_Refs = new TLegend(0.12,0.12,0.6,0.3);
  leg_Refs->SetBorderSize(0.);
  leg_Refs->SetFillColor(0);
  leg_Refs->AddEntry(hRaw_Cum_Refs,"Unsub pPb","p");
  leg_Refs->AddEntry(hBase_Cum_Refs,"pp ","p");
  leg_Refs->AddEntry(hBase_Cum_Refs_int,"pp (0-100%)","p");

  TCanvas* canRefs = new TCanvas("canRefs","canRefs",1200,400);
  canRefs->Divide(3,1);
  canRefs->cd(1);
  TH1* frame_Ref = (TH1*) gPad->DrawFrame(0,0,100,0.01);
  frame_Ref->SetTitle("raw <<2>>; cent %");
  hRaw_Cum_Refs->Draw("same");
  hBase_Cum_Refs->Draw("same");
  hBase_Cum_Refs_int->Draw("same");
  leg_Refs->Draw();
  canRefs->cd(2);
  TH1* frame_Ref_2 = (TH1*) gPad->DrawFrame(0,-1.0,100,10.0);
  frame_Ref_2->SetTitle("<M>^{2} * <<2>>; cent %");
  hRaw_Cum_Refs_scaled->Draw("same");
  hBase_Cum_Refs_scaled->Draw("same");
  hBase_Cum_Refs_int_scaled->Draw("same");
  canRefs->cd(3);
  TH1* frame_Ref_3 = (TH1*) gPad->DrawFrame(0,0,100,10.0);
  frame_Ref_3->SetTitle("<M>^{raw,2} <<2>> - <M>^{base,2}<<2>>; cent %");
  hSub_Cum_Refs->Draw("same");
  hSub_Cum_Refs_int->Draw("same");
  canRefs->SaveAs(Form("%s/cn_subt.pdf",sOutFolder.Data()),"pdf");

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
