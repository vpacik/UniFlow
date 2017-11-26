#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"

TH1D* DoSubtraction(TH1D* hBase, TH1D* hPP);

Color_t gColorBase = kRed;
Color_t gColorPP = kGreen+2;

void FlowSubt()
{
  const char* sInputFile_pp = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0/output_multint/Processed.root";
  const char* sInputFile_pPb = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output/Processed.root";
  const char* sHistoName_pp = "hFlow2_Charged_harm2_gap00_cent0";
  const char* sHistoName_pPb = "hFlow2_Charged_harm2_gap00_cent2";

  // loading files
  TFile* file_pp = TFile::Open(sInputFile_pp,"READ");
  if(!file_pp) { printf("ERROR: File (pp) '%s' does not exists!\n",sInputFile_pp); return; }
  TFile* file_pPb = TFile::Open(sInputFile_pPb,"READ");
  if(!file_pPb) { printf("ERROR: File (pPb) '%s' does not exists!\n",sInputFile_pPb); return; }

  // loading histos
  TH1D* hPP = (TH1D*) file_pp->Get(sHistoName_pp);
  if(!hPP) { printf("ERROR: Histo (pp) '%s' not found!\n",sHistoName_pp); file_pp->ls(); return; }

  TH1D* hPPb = (TH1D*) file_pPb->Get(sHistoName_pPb);
  if(!hPPb) { printf("ERROR: Histo (pPb) '%s' not found!\n",sHistoName_pPb); file_pPb->ls(); return; }

  // do subtraction
  TH1D* hSubted = DoSubtraction(hPPb, hPP);
  if(!hSubted) { printf("ERROR: Subtraction unsuccesfull!\n"); return; }

  // plotting the histos
  TCanvas* can = new TCanvas("can","can");
  can->Divide(2,1);
  can->cd(1);

  hPPb->SetLineColor(gColorBase);
  hPPb->SetMarkerColor(gColorBase);
  hPPb->Draw();

  hPP->SetLineColor(gColorPP);
  hPP->SetMarkerColor(gColorPP);
  hPP->Draw("same");

  can->cd(2);
  hSubted->Draw();

  return;
}

TH1D* DoSubtraction(TH1D* hBase, TH1D* hPP)
{
  if(!hBase) { printf("ERROR-DoSubtraction: Histo (hBase) not found!\n"); return 0x0; }
  if(!hPP) { printf("ERROR-DoSubtraction: Histo (hPP) not found!\n"); return 0x0; }

  TH1D* hSubt = (TH1D*) hBase->Clone(Form("%s_subt",hBase->GetName()));
  if(!hSubt) { printf("ERROR-DoSubtraction: Cloning of hBase not successfull!\n"); return 0x0; }

  return hSubt;
}
