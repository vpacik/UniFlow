// macro for comparing list of histogram

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"

Color_t gColors[] = { kRed, kGreen+2, kBlue, kMagenta+2, kBlack};
Color_t gMarkerStyles[] = { kFullCircle, kOpenCircle, kOpenSquare,kOpenSquare,kOpenSquare };

// unsubtracted
const char* gsFileList[] = {
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0/output_multint_weighted/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0/output_multint_noweight/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0/output_multint_sub_norm/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0/output_multint_sub_norm_weighted/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0/output_multint_sub/Processed.root"
};
// const char* gsHistoNames = "hFlow2_Charged_harm2_gap00_cent3";
const char* gsOutputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0/plots/";
Double_t gdMinimum = 0.;
Double_t gdMaximum = 1.;

// //subtracted
// const char* gsFileList[] = {
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_weighted/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_noweight/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub_norm/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub_norm_weighted/Subtracted_Charged.root"
// };
// const char* gsOutputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/plots_sub/";
// const char* gsHistoNames = "hFlow2_Charged_harm2_gap00_cent3_subt";
// Double_t gdMinimum = -0.15;
// Double_t gdMaximum = 0.2;


const char* gsLabels[] = {
  "GF [A]",
  "GF wo event weights [B]",
  "SP [C]",
  "SP w event weights [D]",
  "SP wo M scaling [E]"
};

Int_t giNumFiles = sizeof(gsFileList)/sizeof(gsFileList[0]);
Int_t iNumCent = 1;

TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);

void CompareRatio(const char** sFileList = gsFileList, Int_t iNumFiles = giNumFiles)
{
  printf("Num files: %d\n",iNumFiles);

  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    const char* sHistoName = Form("hFlow2_Charged_harm2_gap00_cent%d",iCent);
    TList* list = new TList();
    list->SetOwner(kTRUE);

    // finding maximum
    Double_t dMinimum = 0.0;
    Double_t dMaximum = 0.0;

    for(Int_t i(0); i < iNumFiles; ++i)
    {
      printf("Loading form : %s\n",sFileList[i]);

      // loading from root file
      TFile* file = TFile::Open(sFileList[i],"READ");
      if(!file) { printf("ERROR : Cannot open file '%s'!\n",sFileList[i]); continue; }

      TH1D* hist = (TH1D*) file->Get(sHistoName);
      if(!hist) { printf("ERROR : Cannot find histo '%s'!\n",sHistoName); file->ls(); continue; }

      // finding maximum & minimum
      if(hist->GetMaximum() > dMaximum) dMaximum = 1.2 * hist->GetMaximum();
      if(hist->GetMinimum() < dMinimum) dMinimum = hist->GetMinimum();

      list->Add(hist);
    }

    if(dMaximum < gdMaximum) dMaximum = gdMaximum;
    if(dMinimum > 0.) dMinimum = 0.0;

    printf("Min %g \t Max %g \n",dMinimum,dMaximum);

    Int_t iEntries = list->GetEntries();
    printf("Number of entries in list : %d\n",iEntries);

    TLegend* leg = new TLegend(0.15,0.5,0.6,0.88);
    leg->SetTextSize(0.03);

    TCanvas* can = new TCanvas("can","can");
    can->Divide(2,1);
    can->cd(1);

    for(Int_t i(0); i < iEntries; ++i)
    {
      TH1D* hist = (TH1D*) list->At(i);

      hist->SetMaximum(1.*dMaximum);
      hist->SetMinimum(dMinimum);
      // hist->SetMaximum(gdMaximum);
      hist->SetTitle("");
      hist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      hist->GetYaxis()->SetTitle("v_{2}");
      hist->SetStats(kFALSE);
      hist->SetLineColor(gColors[i]);
      hist->SetMarkerColor(gColors[i]);
      hist->SetMarkerStyle(gMarkerStyles[i]);
      leg->AddEntry(hist,gsLabels[i],"pel");
      hist->Draw("same");
    }
    leg->SetFillColorAlpha(0,0);
    leg->SetBorderSize(0);
    leg->Draw();
    // reference histo (first in list)
    TH1D* hBaseline = (TH1D*) list->At(0);

    can->cd(2);
    for(Int_t i(0); i < iEntries; ++i)
    {
      TH1D* hist = (TH1D*) list->At(i);
      TH1D* ratio = (TH1D*) DivideHistos(hist,hBaseline);
      hist->SetTitle("");
      ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      ratio->GetYaxis()->SetTitle("Ratio wrt. standard v_{2}{2}");
      ratio->SetLineColor(gColors[i]);
      ratio->SetMarkerColor(gColors[i]);
      ratio->SetMinimum(0.);
      ratio->SetMaximum(2.);
      ratio->Draw("same");
    }

    // writing
    gSystem->mkdir(gsOutputFolder,kTRUE);

    can->SaveAs(Form("%s/comp_%s.pdf",gsOutputFolder,sHistoName),"pdf");

  }

  return;
}


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
