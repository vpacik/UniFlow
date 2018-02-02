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
TH1D* ScaleCn(TH1D* base, TH1D* mult);
TH1D* Subtract(TH1D* raw, TH1D* base, Double_t factor = 1.0);

// colors for centrality
Color_t colors[] = {kGreen+2, kBlue, kBlack, kMagenta+1};


void Subt_ppb_pp()
{
  TString sSpecies = "Charged";
  TString sMethod = "GF_eventweighted";
  TString sOutputTag = "output_vn";
  TString sOutputTagInt = sOutputTag + "_int";

  TString sGapRaw = "gap08";
  TString sGapBase = "gap08";

  TString sInFileRaw = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-" + sGapRaw + "/" + sOutputTag;
  TString sInFileBase = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run3-2-" + sGapBase + "/" + sOutputTag;
  TString sInFileBaseInt = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run3-2-" + sGapBase + "/" + sOutputTagInt;
  TString sOutFolder = sInFileRaw+"/"+sMethod+"/pPb-pp_test/"+sSpecies;
  TString sOutFile = sOutFolder+"/Subt_results.root";

  const Int_t iNumCent = 4;
  TString sCentLabel[iNumCent] = {"0-20%", "20-40%", "40-60%", "60-100%"};

  // ==================================================================================================================

  // === LOADING INPUT ===
  // multiplicities
  TFile* fileInRaw_Mult = OpenFile(sInFileRaw+"/Mult.root"); if(!fileInRaw_Mult) { return; }
  TProfile* hRaw_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInRaw_Mult); if(!hRaw_Mult) { return; }

  TFile* fileInBase_Mult = OpenFile(sInFileBase+"/Mult.root"); if(!fileInBase_Mult) { return; }
  TProfile* hBase_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_Mult); if(!hBase_Mult) { return; }

  TFile* fileInBase_MultInt = OpenFile(sInFileBaseInt+"/Mult.root"); if(!fileInBase_MultInt) { return; }
  TProfile* hBase_MultInt = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_MultInt); if(!hBase_MultInt) { return; }
  Double_t dMult_Base_Int = hBase_MultInt->GetBinContent(1);

  // output files
  gSystem->mkdir(sOutFolder,kTRUE);  // Making output folder
  TFile* fileOut = OpenFile(sOutFile,"RECREATE"); if(!fileOut) { return; }

  // input files
  TFile* fileInRaw = OpenFile(sInFileRaw+"/"+sMethod+"/Processed.root"); if(!fileInRaw) { return; }
  TFile* fileInBase = OpenFile(sInFileBase+"/"+sMethod+"/Processed.root"); if(!fileInBase) { return; }
  TFile* fileInBaseInt = OpenFile(sInFileBaseInt+"/"+sMethod+"/Processed.root"); if(!fileInBaseInt) { return; }

  // cn
  TH1D* hRaw_cn = LoadHisto(Form("hCum2_Refs_harm2_%s",sGapRaw.Data()),fileInRaw); if(!hRaw_cn) { return; }
  StyleHist(hRaw_cn, kRed, kOpenSquare);
  TH1D* hBase_cn = LoadHisto(Form("hCum2_Refs_harm2_%s",sGapBase.Data()),fileInBase); if(!hBase_cn) { return; }
  StyleHist(hBase_cn, kGreen+2, kFullCircle);
  TH1D* hBase_cn_int = LoadHisto(Form("hCum2_Refs_harm2_%s",sGapBase.Data()),fileInBaseInt); if(!hBase_cn_int) { return; }
  StyleHist(hBase_cn_int, kBlue, kFullCircle);

  TH1D* hRaw_cn_scaled = ScaleCn(hRaw_cn,hRaw_Mult);
  TH1D* hBase_cn_scaled = ScaleCn(hBase_cn,hBase_Mult);
  TH1D* hBase_cn_int_scaled = ScaleCn(hBase_cn_int,hBase_MultInt);

  // dn
  TList* list_Raw_dn = new TList();
  TList* list_Raw_dn_scaled = new TList();
  TList* list_Base_dn = new TList();
  TList* list_Base_dn_scaled = new TList();

  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TH1D* temp = LoadHisto(Form("hCum2_%s_harm2_%s_cent%d",sSpecies.Data(),sGapRaw.Data(),cent),fileInRaw); if(!temp) { return; }
    StyleHist(temp, kRed, kOpenSquare);
    list_Raw_dn->Add(temp);

    temp = LoadHisto(Form("hCum2_%s_harm2_%s_cent%d",sSpecies.Data(),sGapBase.Data(),cent),fileInBase); if(!temp) { return; }
    StyleHist(temp, kGreen+2, kFullCircle);
    list_Base_dn->Add(temp);

    TH1D* raw_scaled = Scale((TH1D*) list_Raw_dn->At(cent), hRaw_Mult->GetBinContent(cent+1)); if(!raw_scaled) { return; }
    list_Raw_dn_scaled->Add(raw_scaled);

    TH1D* base_scaled = Scale((TH1D*) list_Base_dn->At(cent), hBase_Mult->GetBinContent(cent+1)); if(!base_scaled) { return; }
    list_Base_dn_scaled->Add(base_scaled);
  }

  TH1D* hBase_dn_int = LoadHisto(Form("hCum2_%s_harm2_%s_cent0",sSpecies.Data(),sGapBase.Data()),fileInBaseInt); if(!hBase_dn_int) { return; }
  StyleHist(hBase_dn_int, kBlue, kFullCircle);
  TH1D* hBase_dn_int_scaled = Scale(hBase_dn_int, hBase_MultInt->GetBinContent(1));

  fileOut->cd();
  hRaw_cn->Write("hRaw_cn");
  hBase_cn->Write("hBase_cn");
  hBase_cn_int->Write("hBase_cn_int");
  hRaw_cn_scaled->Write("hRaw_cn_scaled");
  hBase_cn_scaled->Write("hBase_cn_scaled");
  hBase_cn_int_scaled->Write("hBase_cn_int_scaled");

  list_Raw_dn->Write("list_Raw_dn",TObject::kSingleKey);
  list_Base_dn->Write("list_Base_dn",TObject::kSingleKey);
  hBase_dn_int->Write("hBase_dn_int");
  list_Raw_dn_scaled->Write("list_Raw_dn_scaled",TObject::kSingleKey);
  list_Base_dn_scaled->Write("list_Base_dn_scaled",TObject::kSingleKey);
  hBase_dn_int_scaled->Write("hBase_dn_int_scaled");

  // === SUBTRACTING pPb(cent) - pp (MB / cent) ===

  // cn{2}^sub = <M>^raw * cn{2}^raw - <M>^base * cn{2}^base
  TH1D* hSub_cn = Subtract(hRaw_cn_scaled, hBase_cn_scaled,1.0); if(!hSub_cn) { return; }
  TH1D* hSub_cn_int = Subtract(hRaw_cn_scaled, hBase_cn_int_scaled,1.0); if(!hSub_cn_int) { return; }

  // dn{2}^sub = <M>dn{2}^raw - <M>dn{2}^base
  TList* list_SubPP_dn = new TList();
  TList* list_SubPP_dn_int = new TList();
  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* temp_raw_cent = (TH1D*) list_Raw_dn_scaled->At(centRaw);
    // subtracting (cent) pp
    TH1D* hSubPP_dn = Subtract(temp_raw_cent, (TH1D*) list_Base_dn_scaled->At(centRaw));
    list_SubPP_dn->Add(hSubPP_dn);

    // subtracting (int) pp
    TH1D* hSubPP_dn_int = Subtract(temp_raw_cent, hBase_dn_int_scaled);
    list_SubPP_dn_int->Add(hSubPP_dn_int);

    TLegend* leg = new TLegend(0.12,0.5,0.6,0.89);
    leg->SetBorderSize(0.);
    leg->SetFillColor(0);
    leg->AddEntry(temp_raw_cent,Form("Unsub pPb (%s)",sCentLabel[centRaw].Data()),"p");
    leg->AddEntry(hBase_dn_int,"pp (0-100%)","p");
    leg->AddEntry(hSubPP_dn,Form("pp (%s)",sCentLabel[centRaw].Data()),"p");

    TCanvas* can = new TCanvas("can","can",1200,400);
    can->Divide(3,1);
    can->cd(1);
    TH1* frame = gPad->DrawFrame(0,0,10,0.2);
    frame->SetTitle("raw <<2'>>; p_{T} (GeV/c)");
    ((TH1D*) list_Raw_dn->At(centRaw))->Draw("same");
    hBase_dn_int->Draw("same");
    ((TH1D*) list_Base_dn->At(centRaw))->Draw("same");
    leg->Draw();

    can->cd(2);
    TH1* frame2 = gPad->DrawFrame(0,-0.03,10,1.0);
    frame2->SetTitle("<M> * <<2'>>; p_{T} (GeV/c)");
    temp_raw_cent->Draw("same");
    hBase_dn_int_scaled->Draw("same");
    ((TH1D*) list_Base_dn_scaled->At(centRaw))->Draw("same");

    can->cd(3);
    TH1* frame3 = gPad->DrawFrame(0,-0.3,10,0.5);
    frame3->SetTitle("<M>^{pPb}<<2'>>^{pPb} - <M>^{pp}<<2'>>^{pp}; p_{T} (GeV/c)");
    hSubPP_dn_int->Draw("same");
    hSubPP_dn->Draw("same");

    can->SaveAs(Form("%s/Subt_pp-pbp_cent%d.pdf",sOutFolder.Data(),centRaw),"pdf");
  }

  // Making vn^sub out of dn^sub / sqrt(cn^sub)
  TList* list_SubPP_vn = new TList();
  TList* list_SubPP_vn_int = new TList();
  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* hSub_dn = (TH1D*) list_SubPP_dn->At(centRaw); if(!hSub_dn) { return; }

    // dividing by cn{2}^sub (cent)
    TH1D* hSubPP_vn = (TH1D*) hSub_dn->Clone(Form("%s_vn", hSub_dn->GetName())); if(!hSubPP_vn) { return; }
    StyleHist(hSubPP_vn, kGreen+2, kFullCircle);
    hSubPP_vn->Scale(1.0/TMath::Sqrt(hSub_cn->GetBinContent(centRaw+1)));
    list_SubPP_vn->Add(hSubPP_vn);

    // dividing by cn{2}^sub (int)
    TH1D* hSubPP_vn_int = (TH1D*) hSub_dn->Clone(Form("%s_vn_int", hSub_dn->GetName())); if(!hSubPP_vn_int) { return; }
    StyleHist(hSubPP_vn_int, kBlue, kFullCircle);
    hSubPP_vn_int->Scale(1.0/TMath::Sqrt(hSub_cn_int->GetBinContent(1)));
    list_SubPP_vn_int->Add(hSubPP_vn_int);

    TCanvas* canVn = new TCanvas("canVn","canVn",400,400);
    // canVn->Divide(,1);
    canVn->cd(1);
    TH1* frame_vn = (TH1*) gPad->DrawFrame(0.,0.,10.,0.2);
    frame_vn->SetTitle("sub vn{2}; p_{T} (GeV/c)");
    hSubPP_vn->Draw("same");
    hSubPP_vn_int->Draw("same");
    // canVn->cd(2);
    // canVn->cd(3);

    canVn->SaveAs(Form("%s/Subt_ppb_pp_vn_cent%d.pdf",sOutFolder.Data(),centRaw),"pdf");

  }

  fileOut->cd();
  hSub_cn->Write("hSub_cn");
  hSub_cn_int->Write("hSub_cn_int");
  list_SubPP_dn->Write("list_SubtPP_dn_ppcent",TObject::kSingleKey);
  list_SubPP_dn_int->Write("list_SubtPP_dn_ppint",TObject::kSingleKey);
  list_SubPP_vn->Write("list_SubtPP_vn_ppcent",TObject::kSingleKey);
  list_SubPP_vn_int->Write("list_SubtPP_vn_ppint",TObject::kSingleKey);


  // === SUBTRACTING pPb(cent) - pPb (peripheral) ===

  // cn
  TH1D* hSubPPb_cn = (TH1D*) hRaw_cn_scaled->Clone(Form("%s_sub",hRaw_cn_scaled->GetName()));
  for(Int_t bin(1); bin < hRaw_cn_scaled->GetNbinsX()+1; ++bin)
  {
    Double_t base_con = hRaw_cn_scaled->GetBinContent(4);
    Double_t base_err = hRaw_cn_scaled->GetBinError(4);
    Double_t con = hRaw_cn_scaled->GetBinContent(bin);
    Double_t err = hRaw_cn_scaled->GetBinError(bin);
    hSubPPb_cn->SetBinContent(bin, con - base_con);
    hSubPPb_cn->SetBinError(bin,TMath::Sqrt(err*err + base_err*base_err));
  }

  TCanvas* can_pPb_ref = new TCanvas("can_pPb_ref","can_pPb_ref",1200,400);
  can_pPb_ref->Divide(3,1);
  can_pPb_ref->cd(1);
  TH1* frame_pPb_ref_1 = (TH1*) gPad->DrawFrame(0.,0.,100.,0.01);
  frame_pPb_ref_1->SetTitle("<<2>>; cent %");
  hRaw_cn->Draw("same");
  can_pPb_ref->cd(2);
  TH1* frame_pPb_ref_2 = (TH1*) gPad->DrawFrame(0.,0.,100.,10.0);
  frame_pPb_ref_2->SetTitle("<M>^{2}<<2>>; cent %");
  hRaw_cn_scaled->Draw("same");
  can_pPb_ref->cd(3);
  TH1* frame_pPb_ref_3 = (TH1*) gPad->DrawFrame(0.,0.,100.,10.);
  frame_pPb_ref_3->SetTitle("<M>^{raw,2}<<2>>^{raw} - <M>^{base,2}<<2>>^{base}; cent %");
  hSubPPb_cn->Draw("same");
  can_pPb_ref->SaveAs(Form("%s/Subt_ppb_ppb_cn.pdf",sOutFolder.Data()),"pdf");

  // dn
  TList* list_SubtPPb_dn = new TList();
  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TH1D* hSubPPb = Subtract((TH1D*)list_Raw_dn_scaled->At(cent),(TH1D*)list_Raw_dn_scaled->At(3));
    list_SubtPPb_dn->Add(hSubPPb);

    TCanvas* can_pPb_dn = new TCanvas("can_pPb_dn","can_pPb_dn",1200,400);
    can_pPb_dn->Divide(3,1);
    can_pPb_dn->cd(1);
    TH1* frame_pPb = (TH1*) gPad->DrawFrame(0.,0.0,10.,0.2);
    frame_pPb->SetTitle("<<2'>>; p_{T} (GeV/c)");
    ((TH1D*) list_Raw_dn->At(cent))->Draw("same");
    ((TH1D*) list_Raw_dn->At(3))->Draw("same");
    can_pPb_dn->cd(2);
    TH1* frame_pPb_2 = (TH1*) gPad->DrawFrame(0.,0.0,10.,1.);
    frame_pPb_2->SetTitle("<M><<2'>>; p_{T} (GeV/c)");
    ((TH1D*) list_Raw_dn_scaled->At(cent))->Draw("same");
    ((TH1D*) list_Raw_dn_scaled->At(3))->Draw("same");
    can_pPb_dn->cd(3);
    TH1* frame_pPb_3 = (TH1*) gPad->DrawFrame(0.,0.0,10.,1.);
    frame_pPb_3->SetTitle("<M>^{raw}<<2'>>^{raw} - <M>^{base}<<2'>>^{base} ; p_{T} (GeV/c)");
    hSubPPb->Draw("same");
    can_pPb_dn->SaveAs(Form("%s/Subt_ppb_ppb_dn_cent%d.pdf",sOutFolder.Data(),cent),"pdf");
  }

  // vn
  TList* list_SubtPPb_vn = new TList();
  for(Int_t cent(0); cent < iNumCent-1; ++cent)
  {
    TH1D* hSubPPb_vn = (TH1D*) ((TH1D*) list_SubtPPb_dn->At(cent))->Clone(Form("hSubPPb_vn_cent%d",cent));
    hSubPPb_vn->Scale(1.0/TMath::Sqrt(hSubPPb_cn->GetBinContent(cent+1)));
    list_SubtPPb_vn->Add(hSubPPb_vn);

    TCanvas* can_pPb_vn = new TCanvas("can_pPb_vn","can_pPb_vn",400,400);
    can_pPb_vn->cd();
    TH1* frame_pPb_vn = (TH1*) gPad->DrawFrame(0.,0.,10.,0.2);
    frame_pPb_vn->SetTitle("v_{2}{2}^{sub}; p_{T} (GeV/c)");
    hSubPPb_vn->Draw("same");
    can_pPb_vn->SaveAs(Form("%s/Subt_ppb_ppb_vn_cent%d.pdf",sOutFolder.Data(),cent),"pdf");
  }
  fileOut->cd();
  hSubPPb_cn->Write("list_SubtPPb_cn");
  list_SubtPPb_dn->Write("list_SubtPPb_dn",TObject::kSingleKey);
  list_SubtPPb_vn->Write("list_SubtPPb_vn",TObject::kSingleKey);

  // UNIVERSALL PLOTTING
  // RFPs
  TLegend* leg_Refs = new TLegend(0.12,0.12,0.6,0.3);
  leg_Refs->SetBorderSize(0.);
  leg_Refs->SetFillColor(0);
  leg_Refs->AddEntry(hRaw_cn,"Unsub pPb","p");
  leg_Refs->AddEntry(hBase_cn,"pp ","p");
  leg_Refs->AddEntry(hBase_cn_int,"pp (0-100%)","p");

  TCanvas* canRefs = new TCanvas("canRefs","canRefs",1200,400);
  canRefs->Divide(3,1);
  canRefs->cd(1);
  TH1* frame_Ref = (TH1*) gPad->DrawFrame(0,0,100,0.01);
  frame_Ref->SetTitle("raw <<2>>; cent %");
  hRaw_cn->Draw("same");
  hBase_cn->Draw("same");
  hBase_cn_int->Draw("same");
  leg_Refs->Draw();
  canRefs->cd(2);
  TH1* frame_Ref_2 = (TH1*) gPad->DrawFrame(0,-1.0,100,10.0);
  frame_Ref_2->SetTitle("<M>^{2} * <<2>>; cent %");
  hRaw_cn_scaled->Draw("same");
  hBase_cn_scaled->Draw("same");
  hBase_cn_int_scaled->Draw("same");
  canRefs->cd(3);
  TH1* frame_Ref_3 = (TH1*) gPad->DrawFrame(0,0,100,10.0);
  frame_Ref_3->SetTitle("<M>^{raw,2} <<2>> - <M>^{base,2}<<2>>; cent %");
  hSub_cn->Draw("same");
  hSub_cn_int->Draw("same");
  canRefs->SaveAs(Form("%s/cn_subt.pdf",sOutFolder.Data()),"pdf");


  // Comparison of various methods
  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TLegend* leg_comp = new TLegend(0.12,0.6,0.4,0.89);
    leg_comp->SetBorderSize(0);
    leg_comp->SetFillColor(0);
    leg_comp->SetFillStyle(0);
    leg_comp->SetHeader(Form("%s (V0A)", sCentLabel[cent].Data()));
    TCanvas* can_comp = new TCanvas("can_comp","can_comp",400,400);
    can_comp->cd();
    TH1* frame_comp = (TH1*) gPad->DrawFrame(0.,0.,10.,0.2);
    frame_comp->SetTitle("h^{#pm} v_{2}{2}^{sub}; p_{T} (Gev/c)");
    ((TH1D*)list_SubPP_vn->At(cent))->Draw("same");
    leg_comp->AddEntry((TH1D*)list_SubPP_vn->At(cent),Form("pp (%s)",sCentLabel[cent].Data()),"p");
    ((TH1D*)list_SubPP_vn_int->At(cent))->Draw("same");
    leg_comp->AddEntry((TH1D*)list_SubPP_vn_int->At(cent),"pp (MB)","p");
    if(cent < 3)
    {
      ((TH1D*)list_SubtPPb_vn->At(cent))->Draw("same");
      leg_comp->AddEntry((TH1D*)list_SubtPPb_vn->At(cent),"pPb (60-100%)","p");
    }

    leg_comp->Draw();
    can_comp->SaveAs(Form("%s/comp_cent%d.pdf",sOutFolder.Data(),cent),"pdf");
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
TH1D* ScaleCn(TH1D* base, TH1D* mult)
{
  if(!base) { printf("ERROR-ScaleCn: Hist 'base' does not found.\n"); return 0x0; }
  if(!mult) { printf("ERROR-ScaleCn: Hist 'mult' does not found.\n"); return 0x0; }

  TH1D* scaled = (TH1D*) base->Clone(Form("%s_scaled",base->GetName()));

  for(Int_t bin(1); bin < base->GetNbinsX()+1; ++bin)
  {
    Double_t con =  base->GetBinContent(bin);
    Double_t err =  base->GetBinError(bin);
    Double_t mult_con =  mult->GetBinContent(bin);
    Double_t mult_err =  mult->GetBinError(bin);

    scaled->SetBinContent(bin, con * mult_con*mult_con );
    scaled->SetBinError(bin, TMath::Power( mult_con*mult_con*err, 2.0) + TMath::Power(2.0*con*mult_con*mult_err, 2.0));
  }

  return scaled;
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
    // sub->GetBinError(bin, 0.01*con_raw);
    sub->GetBinError(bin, err_raw*err_raw + factor*factor*err_base*err_base );

  }
  return sub;
}
// ==================================================================================================================
