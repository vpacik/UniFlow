#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

TH1D* DoSubtraction(TH1D* hBase, TH1D* hPP, Double_t dNorm = 0.5);


Color_t gColorBase = kRed;
Color_t gColorPP = kGreen+2;

void FlowSubt()
{
  Int_t iNumCent = 4;
  TString sGap = "00";

  // TString sTag = "weighted";
  // TString sTag = "noweight";
  // TString sTag = "sub";
  // TString sTag = "sub_norm";
  TString sTag = "sub_norm_weighted";

  Int_t iNumSpecies = 7;
  TString sSpecies[] = {"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};

  // input files
  TString sInputDir_pPb = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap00/";
  TString sInputDir_pp = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap00/";

  // NOTE: potential loop over tags

  // output files files
  TString sOutputDir = Form("%s/output_%s/flow-subtraction",sInputDir_pPb.Data(),sTag.Data());
  TString sOutputFile = Form("Subtracted.root");
  TString sOutputFileRaw = Form("Raw.root");
  TString sOutputFilePPscaled = Form("ScaledPP.root");

  gSystem->mkdir(sOutputDir.Data(),kTRUE);
  TFile* fOutputFile = TFile::Open(Form("%s/%s",sOutputDir.Data(),sOutputFile.Data()),"RECREATE");
  if(!fOutputFile) { printf("ERROR: Output file not found.\n"); return; }
  TFile* fOutputFileRaw = TFile::Open(Form("%s/%s",sOutputDir.Data(),sOutputFileRaw.Data()),"RECREATE");
  if(!fOutputFileRaw) { printf("ERROR: Output file (Raw) not found.\n"); return; }
  TFile* fOutputFileScaledPP = TFile::Open(Form("%s/%s",sOutputDir.Data(),sOutputFilePPscaled.Data()),"RECREATE");
  if(!fOutputFileRaw) { printf("ERROR: Output file (ScaledPP) not found.\n"); return; }

  // processing normalisation
  // normalisation
  TString sMultFile_pp = Form("%s/Mult_int.root",sInputDir_pp.Data());
  TString sMultFile_pPb = Form("%s/Mult.root",sInputDir_pPb.Data());
  TFile* fMult_pp = TFile::Open(sMultFile_pp.Data(),"READ");
  if(!fMult_pp) { printf("ERROR: Mult file (pp) '%s' does not exists!\n",sMultFile_pp.Data()); return; }
  TFile* fMult_pPb = TFile::Open(sMultFile_pPb.Data(),"READ");
  if(!fMult_pPb) { printf("ERROR: Mult file (pPb) '%s' does not exists!\n",sMultFile_pPb.Data()); return; }

  TH1D* hMultPP = (TH1D*) fMult_pp->Get("fpRefsMult_rebin");
  if(!hMultPP) { printf("ERROR: Mult histo (pp) does not exists\n"); fMult_pp->ls(); return; }

  TH1D* hMultPPb = (TH1D*) fMult_pPb->Get("fpRefsMult_rebin");
  if(!hMultPPb) { printf("ERROR: Mult histo (Base) does not exists\n"); fMult_pPb->ls(); return; }

  TLegend* legMult = new TLegend(0.7,0.7,0.89,0.89);
  legMult->SetBorderSize(0);
  legMult->AddEntry(hMultPPb,"pPb","lp");
  legMult->AddEntry(hMultPP,"pp","lp");


  TCanvas* canNorm = new TCanvas("canNorm","canNorm");
  canNorm->cd();
  hMultPPb->SetStats(0);
  hMultPPb->SetTitle("Mean RFPs multiplicity");
  hMultPPb->GetYaxis()->SetTitle("<M>");
  hMultPPb->GetXaxis()->SetTitle("centrality %");
  hMultPPb->SetMinimum(0.);
  hMultPPb->SetLineColor(gColorBase);
  hMultPPb->SetMarkerColor(gColorBase);
  hMultPPb->SetMarkerStyle(kFullCircle);
  hMultPPb->Draw();
  hMultPP->SetLineColor(gColorPP);
  hMultPP->SetMarkerColor(gColorPP);
  hMultPP->SetMarkerStyle(kFullCircle);
  hMultPP->Draw("same");
  legMult->Draw();
  canNorm->SaveAs(Form("%s/normalisation.pdf",sOutputDir.Data()),"pdf");

  // loop over species
  for(Int_t iSpecies(0); iSpecies < iNumSpecies; ++iSpecies)
  {
    TString sInputFile_pp = Form("%s/output_%s/Processed.root",sInputDir_pp.Data(),sTag.Data());
    TString sInputFile_pPb = Form("%s/output_%s/Processed.root",sInputDir_pPb.Data(),sTag.Data());

    TFile* file_pp = TFile::Open(sInputFile_pp.Data(),"READ");
    if(!file_pp) { printf("ERROR: File (pp) '%s' does not exists!\n",sInputFile_pp.Data()); return; }
    TFile* file_pPb = TFile::Open(sInputFile_pPb.Data(),"READ");
    if(!file_pPb) { printf("ERROR: File (pPb) '%s' does not exists!\n",sInputFile_pPb.Data()); return; }


    for(Int_t iCent(0); iCent < iNumCent; ++iCent)
    {
      TString sHistoName_pp = Form("hFlow2_%s_harm2_gap%s",sSpecies[iSpecies].Data(),sGap.Data());
      TString sHistoName_pPb = Form("hFlow2_%s_harm2_gap%s",sSpecies[iSpecies].Data(),sGap.Data()); // wihtout centrality index

      if( sSpecies[iSpecies].EqualTo("Charged") ||
          sSpecies[iSpecies].EqualTo("Pion") ||
          sSpecies[iSpecies].EqualTo("Kaon") ||
          sSpecies[iSpecies].EqualTo("Proton") )
      {
         sHistoName_pp.Append("_cent0");
         sHistoName_pPb.Append(Form("_cent%d",iCent));
      }
      else
      {
         sHistoName_pp.Append("_mult0");
         sHistoName_pPb.Append(Form("_mult%d",iCent));
      }

      // loading histos
      TH1D* hPP = (TH1D*) file_pp->Get(sHistoName_pp.Data());
      if(!hPP) { printf("ERROR: Histo (pp) '%s' not found!\n",sHistoName_pp.Data()); file_pp->ls(); return; }

      TH1D* hPPb = (TH1D*) file_pPb->Get(sHistoName_pPb.Data());
      if(!hPPb) { printf("ERROR: Histo (pPb) '%s' not found!\n",sHistoName_pPb.Data()); file_pPb->ls(); return; }

      // estimating normalisation
      Double_t dMultPPb = hMultPPb->GetBinContent(iCent+1);
      Double_t dMultPP = hMultPP->GetBinContent(1);
      Double_t dNorm = dMultPP / dMultPPb;

      // do subtraction
      TH1D* hSubted = DoSubtraction(hPPb, hPP, dNorm);
      if(!hSubted) { printf("ERROR: Subtraction unsuccesfull!\n"); return; }

      TH1D* hPP_Scaled = (TH1D*) hPP->Clone(Form("%s_scaled_mult%d",hPP->GetName(),iCent));
      hPP_Scaled->Scale(dNorm);

      // plotting the histos

      TLegend* legFlow = new TLegend(0.15,0.7,0.4,0.88);
      legFlow->SetBorderSize(0);
      legFlow->AddEntry(hPP,"pp (raw)","pl");
      legFlow->AddEntry(hPPb,"pPb (raw)","pl");
      legFlow->AddEntry(hSubted,"pPb (sub)","pl");

      // flow
      TCanvas* can = new TCanvas("can","can");
      // can->Divide(2,1);
      can->cd(1);

      hPP->SetStats(0);
      hPP->SetTitle(Form("%s v_{2}(p_{T})",sSpecies[iSpecies].Data()));
      hPP->GetYaxis()->SetTitle("v_{2}");
      hPP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hPP->SetMinimum(-0.1);
      hPP->SetMaximum(1.0);
      hPP->SetLineColor(gColorPP);
      hPP->SetMarkerColor(gColorPP);
      hPP->SetMarkerStyle(kFullCircle);
      hPP->Draw("same");
      hPPb->SetLineColor(gColorBase);
      hPPb->SetMarkerColor(gColorBase);
      hPPb->SetMarkerStyle(kFullCircle);
      hPPb->Draw("same");

      hSubted->SetLineColor(kBlue);
      hSubted->SetMarkerColor(kBlue);
      hSubted->SetMarkerStyle(kFullSquare);

      hSubted->Draw("same");

      legFlow->Draw();

      // can->cd(2);
      // hSubted->SetMinimum(0.0);
      // hSubted->Draw();

      can->SaveAs(Form("%s/flow_subtraction_%s_gap%s_cent%d.pdf",sOutputDir.Data(),sSpecies[iSpecies].Data(),sGap.Data(),iCent),"pdf");

      // output file
      fOutputFile->cd();
      hSubted->Write();

      fOutputFileRaw->cd();
      hPPb->Write();

      fOutputFileScaledPP->cd();
      hPP_Scaled->Write();
    }
  }


  return;
}

TH1D* DoSubtraction(TH1D* hBase, TH1D* hPP, Double_t dNorm)
{

  if(!hBase) { printf("ERROR-DoSubtraction: Histo (hBase) not found!\n"); return 0x0; }
  if(!hPP) { printf("ERROR-DoSubtraction: Histo (hPP) not found!\n"); return 0x0; }

  TH1D* hSubt = (TH1D*) hBase->Clone(Form("%s_subt",hBase->GetName()));
  if(!hSubt) { printf("ERROR-DoSubtraction: Cloning of hBase not successfull!\n"); return 0x0; }

  for(Int_t ibin(1); ibin < hSubt->GetNbinsX()+1; ++ibin)
  {
    Double_t dConBase = hBase->GetBinContent(ibin);
    Double_t dConPP = hPP->GetBinContent(ibin);

    Double_t dErrBase = hBase->GetBinError(ibin);
    Double_t dErrPP = hPP->GetBinError(ibin);

    Double_t dErr2 = dErrBase*dErrBase + TMath::Power(dNorm*dErrPP,2.);

    hSubt->SetBinContent(ibin, dConBase - dNorm * dConPP);
    hSubt->SetBinError(ibin, TMath::Sqrt(dErr2) );
  }

  return hSubt;
}
