void TransformQAtoSparse()
{
  const Int_t iDim = 40;
  const Int_t iNumBins = 200;
  Int_t iNumEntries = 3*iDim*iNumBins;


  TFile* fileIn = TFile::Open("./AnalysisResults.root", "READ");
  if(!fileIn) { printf("E: fileIn not found!\n"); return; }

  TFile* fileOutOrig = TFile::Open("./outQA_orig.root", "RECREATE");
  if(!fileOutOrig) { printf("E: fileOutOrig not found!\n"); return; }

  TFile* fileOutSparse = TFile::Open("./outQA_sparse.root", "RECREATE");
  if(!fileOutSparse) { printf("E: fileOutSparse not found!\n"); return; }

  TRandom3 random = TRandom3();

  TString name = "sparse; ";
  Int_t bins[iDim] = {};
  Double_t min[iDim] = {};
  Double_t max[iDim] = {};
  Double_t mean[iDim] = {};
  Double_t sigma[iDim] = {};
  TH1D* histos[iDim] = { 0x0 };

  for(Int_t i(0); i < iDim; ++i)
  {
    name = TString(Form("hist_%d; ",i));
    bins[i] = iNumBins;
    min[i] = -20.0;
    max[i] = 20.0;
    mean[i] = random.Uniform(40) - 20.0;
    sigma[i] = random.Uniform(5);
    histos[i] =  new TH1D(name.Data(),name.Data(),bins[i], min[i], max[i]);


  }

  THnSparseD* sparse = new THnSparseD("sparse","sparse",iDim,bins,min,max);


  // filling
  for(Int_t e(0); e < iNumEntries; ++e)
  {
    Double_t numbers[iDim] = {0};

    for(Int_t i(0); i < iDim; ++i)
    {
      Double_t number = random.Gaus(mean[i], sigma[i]);
      numbers[i] = number;
      histos[i]->Fill(number);
    }

    sparse->Fill(numbers);
  }

  fileOutOrig->cd();
  for(Int_t i(0); i < iDim; ++i)
  {
    histos[i]->Write();
  }

  // saving sparse to different file
  fileOutSparse->cd();
  sparse->Write();


  printf("Mem %g\n",sparse->GetSparseFractionMem());
  printf("Bins %g\n",sparse->GetSparseFractionBins());

}
