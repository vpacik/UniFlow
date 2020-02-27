void RunPlotFitsSingle()
{
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/6815/thesis/";
    TString sCorr; TString sCorrTitle;

    Int_t iCent = 1;
    Bool_t bPhi = 0;
    // TString sSpecies = "K0s"; Int_t iPt = 10;
    // TString sSpecies = "Lambda"; Int_t iPt = 5;
    TString sSpecies = "Phi"; Int_t iPt = 3; bPhi = 1;

    sCorr = "<<2>>(2,-2)"; sCorrTitle = "#LT#LT2'#GT#GT_{2,-2}";
    // gROOT->Macro(Form("PlotFitsSingleCorr.C(\"%s\",\"%s\",\"%s\",\"%s\",%d,%d)",sPath.Data(),sCorr.Data(),sCorrTitle.Data(),sSpecies.Data(),iPt,iCent));
    //
    // sCorr = "<<4>>(2,2,-2,-2)"; sCorrTitle = "#LT#LT4'#GT#GT_{2,2,-2,-2}";
    // gROOT->Macro(Form("PlotFitsSingleCorr.C(\"%s\",\"%s\",\"%s\",\"%s\",%d,%d)",sPath.Data(),sCorr.Data(),sCorrTitle.Data(),sSpecies.Data(),iPt,iCent));

    gROOT->Macro(Form("PlotFitsSingleMass.C(\"%s\",\"%s\",\"%s\",%d,%d,%d)",sPath.Data(),sCorr.Data(),sSpecies.Data(),iPt,iCent,bPhi));
}
