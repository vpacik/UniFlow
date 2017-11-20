// macro for comparing list of histogram

Color_t gColors[] = { kRed, kGreen+2, kBlue, kMagenta+2 };
Color_t gMarkerStyles[] = { kFullCircle, kOpenCircle, kOpenCircle };

const char* gsOutputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-test/plots";
const char* gsHistoNames = "hFlow2_Lambda_harm2_gap08_mult3";

const char* gsFileList[] = {
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-test/output_old/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-test/output_sub/Processed.root",
  "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-test/output_sub_norm/Processed.root"
};

const char* gsLabels[] = {
  "v2{2}",
  "SP",
  "SP norm"
};

Int_t giNumFiles = sizeof(gsFileList)/sizeof(gsFileList[0]);


TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);

void CompareRatio(const char* sHistoName = gsHistoNames, const char** sFileList = gsFileList, Int_t iNumFiles = giNumFiles)
{
  printf("Num files: %d\n",iNumFiles);

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

  printf("Min %g \t Max %g \n",dMinimum,dMaximum);

  Int_t iEntries = list->GetEntries();
  printf("Number of entries in list : %d\n",iEntries);

  TLegend* leg = new TLegend(0.6,0.2,0.89,0.4);

  TCanvas* can = new TCanvas("can","can");
  can->Divide(2,1);
  can->cd(1);


  for(Int_t i(0); i < iEntries; ++i)
  {
    TH1D* hist = (TH1D*) list->At(i);


    hist->SetMaximum(1.*dMaximum);

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
    ratio->SetLineColor(gColors[i]);
    ratio->SetMarkerColor(gColors[i]);
    ratio->SetMinimum(0.8);
    ratio->SetMaximum(1.2);
    ratio->Draw("same");
  }

  // writing
  gSystem->mkdir(gsOutputFolder,kTRUE);

  can->SaveAs(Form("%s/comp_%s.pdf",gsOutputFolder,sHistoName),"pdf");

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
