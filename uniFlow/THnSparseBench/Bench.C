void Bench()
{
  TFile* file = TFile::Open("AnalysisResults.root","READ");
  file->cd("UniFlow");
  file->ls();
  TList* list = gDirectory->Get("Flow_K0s_UniFlow");

  TH3D* hist = (TH3D*) list->FindObject("fh3V0sEntriesK0s_gap08");

  hist->Draw();

  Int_t iBinsX = hist->GetNbinsX();
  Int_t iBinsY = hist->GetNbinsY();
  Int_t iBinsZ = hist->GetNbinsZ();
  Int_t iBinsTot = iBinsX*iBinsY*iBinsZ;
  Int_t iBinsEmpty = 0;

  printf("Bins %d | X %d | Y %d | Z %d\n",iBinsTot,iBinsX,iBinsY,iBinsZ);
  for(Int_t x(0); x < iBinsX+2; ++x)
    for(Int_t y(0); y < iBinsY+2; ++y)
      for(Int_t z(0); z < iBinsZ+2; ++z)
      {
        if(hist->GetBinContent(x,y,z) == 0) iBinsEmpty++;
      }


   THnSparseD* sparse = THnSparse::CreateSparse("sparse","sparse",hist);
   sparse->Draw();


   printf("Fraction of empty bins %g  (%d out of %d) | %g\n",(float) iBinsEmpty/iBinsTot,iBinsEmpty,iBinsTot,sparse->GetSparseFractionBins());
   printf("Fraction of memory %g\n",sparse->GetSparseFractionMem());
   printf("Size %g (sparse) | %g (hist)\n",sizeof(*sparse), sizeof(*hist));

   TFile* outSparse = TFile::Open("outSparse.root","RECREATE");
   outSparse->cd();
   sparse->Write("sparse");
   sparse->Write("sparse_2");
   sparse->Write("sparse_3");
   sparse->Write("sparse_3");
   sparse->Write("sparse_3");

   TFile* outHist= TFile::Open("outHist.root","RECREATE");
   hist->Write("hist");
   hist->Write("hist_2");
   hist->Write("hist_3");
   hist->Write("hist_3");
   hist->Write("hist_3");
}
