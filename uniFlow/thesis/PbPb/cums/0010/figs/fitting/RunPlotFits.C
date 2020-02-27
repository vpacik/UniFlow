void RunPlotFits()
{
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/6815/thesis/";
    TString sCorr = "<<2>>(2,-2)";
    // TString sCorr = "<<4>>(2,2,-2,-2)";

    Int_t iNumCent = 6;
    // TString sSpecies = "K0s"; Int_t iNumPt = 14;
    // TString sSpecies = "Lambda"; Int_t iNumPt = 13;
    TString sSpecies = "Phi"; Int_t iNumPt = 5;

    gROOT->Macro(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/macros/fits/PlotFits.C(\"%s\",\"%s\",\"%s\",%d,%d)",sPath.Data(),sCorr.Data(),sSpecies.Data(),iNumPt,iNumCent));

}
