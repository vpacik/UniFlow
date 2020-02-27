void FlowFluctAll(TString sPath = "~/Codes/ALICE/Flow/uniFlow/thesis/PbPb/cums/0010/v2/process/", Bool_t bCor = 1,  Int_t iNumCent = 6);

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>


void FlowFluct(
    TString sSpecies = "Charged",
    Int_t iCent = 0,
    TString sPath = "~/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/MC_2210/species_direct/",
    Bool_t bCor = 0,
    TString sHistoTwo = "hFlow2_harm2_gap08",
    TString sHistoFour = "hFlow4_harm2_gap-10",
    TString sOutFileName = "FlowFluct"
)
{
    if(bCor) { sOutFileName.Append("Cor"); }

    printf("%s | %d\n", sSpecies.Data(), iCent);

    TFile* fileIn = TFile::Open(Form("%s/Processed.root",sPath.Data()),"READ");

    TFile* fileOut = TFile::Open(Form("%s/%s.root",sPath.Data(), sOutFileName.Data()),"UPDATE");

    TString sNameTwo = Form("%s_%s_cent%d",sSpecies.Data(), sHistoTwo.Data(), iCent);
    TString sNameFour = Form("%s_%s_cent%d",sSpecies.Data(), sHistoFour.Data(), iCent);

    fileIn->ls();

    TH1D* hist2 = (TH1D*) fileIn->Get(sNameTwo.Data());
    if(!hist2) { printf("Hist2 '%s' not found!\n",sNameTwo.Data()); /*fileIn->ls();*/ return; }

    TH1D* hist4 = (TH1D*) fileIn->Get(sNameFour.Data());
    if(!hist4) { printf("Hist4 '%s' not found!\n",sNameFour.Data()); /*fileIn->ls();*/ return; }

    // make mean
    TH1D* histMean = (TH1D*) hist2->Clone(Form("%s_Mean_cent%d",sSpecies.Data(),iCent));
    if(!histMean) { printf("Mean not created!\n"); return; }
    histMean->Reset();

    // make std
    TH1D* histSTD = (TH1D*) hist2->Clone(Form("%s_Std_cent%d",sSpecies.Data(),iCent));
    if(!histSTD) { printf("STD not created!\n"); return; }
    histSTD->Reset();

    // make relative
    TH1D* histRel = (TH1D*) hist2->Clone(Form("%s_Rel_cent%d",sSpecies.Data(),iCent));
    if(!histRel) { printf("Rel not created!\n"); return; }
    histRel->Reset();

    Int_t iNumBins = hist2->GetNbinsX();
    for(Int_t iBin(1); iBin < iNumBins+1; ++iBin) {
        Double_t c2 = hist2->GetBinContent(iBin);
        Double_t e2 = hist2->GetBinError(iBin);
        Double_t c4 = hist4->GetBinContent(iBin);
        Double_t e4 = hist4->GetBinError(iBin);

        if(c2 == 0.0 || c4 == 0.0) { continue; }

        Double_t mean = (c2*c2 + c4*c4) * 0.5;

        Double_t dErrMean_c2 = TMath::Power(mean,-0.5)*c2*e2;
        Double_t dErrMean_c4 = TMath::Power(mean,-0.5)*c4*e4;

        Double_t err_mean = TMath::Power(dErrMean_c2,2.0) + TMath::Power(dErrMean_c4,2.0);
        if(bCor) { err_mean = TMath::Power(dErrMean_c2 + dErrMean_c4,2.0); }

        // Double_t err_mean = 2.0*(c2*c2*e2*e2 + c4*c4*e4*e4)/(c2*c2 + c4*c4);
        // Double_t err_mean = (c2*c2*e2*e2+c4*c4*e4*e4)*0.25/(c2*c2+c4*c4);

        histMean->SetBinContent(iBin, TMath::Sqrt(mean));
        histMean->SetBinError(iBin, TMath::Sqrt(err_mean));

        if(c2 < c4) { continue; }

        Double_t std = (c2*c2 - c4*c4) * 0.5;

        Double_t dErrStd_c2 = TMath::Power(std,-0.5)*c2*e2;
        Double_t dErrStd_c4 = -1.0*TMath::Power(std,-0.5)*c4*e4;

        Double_t err_std = TMath::Power(dErrStd_c2,2.0) + TMath::Power(dErrStd_c4,2.0);
        if(bCor) { err_std = TMath::Power(dErrStd_c2 + dErrStd_c4,2.0); }
        // Double_t err_std = 2.0*(c2*c2*e2*e2 + c4*c4*e4*e4)/(c2*c2 - c4*c4);
        // Double_t err_std = err_mean;

        histSTD->SetBinContent(iBin, TMath::Sqrt(std));
        histSTD->SetBinError(iBin, TMath::Sqrt(err_std));

        // if(mean == 0.0) { continue; }
        // Double_t rel = std / mean;
        // histRel->SetBinContent(iBin, rel);
    }

    histRel = Utils::DivideHistos(histSTD,histMean,1);
    histRel->Divide(histSTD,histMean);

    fileOut->cd();
    histMean->Write(Form("%s_mean_cent%d",sSpecies.Data(),iCent),TObject::kOverwrite);
    histSTD->Write(Form("%s_std_cent%d",sSpecies.Data(),iCent),TObject::kOverwrite);
    histRel->Write(Form("%s_rel_cent%d",sSpecies.Data(),iCent),TObject::kOverwrite);
    fileOut->Close();

    // TCanvas* can = new TCanvas();
    // can->Divide(3,1);
    // can->cd(1);
    // histMean->Draw();
    // can->cd(2);
    // histSTD->Draw();
    // can->cd(3);
    // histRel->Draw();


    return;
}

void FlowFluctAll(TString sPath, Bool_t bCor, Int_t iNumCent)
{
    std::vector<TString> sSpecies = {"Charged","Pion","Kaon","Proton","K0s","Lambda"};
    // std::vector<TString> sSpecies = {"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};

    Int_t iNumSpec = sSpecies.size();

    for(Int_t iS(0); iS < iNumSpec; ++iS) {
        TString sSpe = sSpecies.at(iS);

        for(Int_t iCent(0); iCent < iNumCent; ++iCent) {
            FlowFluct(sSpe,iCent,sPath,bCor);
        }
    }
}
