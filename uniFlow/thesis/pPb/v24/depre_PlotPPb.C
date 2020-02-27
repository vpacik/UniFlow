// Macro for comparing p-Pb results of v2{4} with v2{2}^{sub} (Preliminaries)

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>
#include <TString.h>
#include <TFile.h>


TGraphErrors* Transfer(TGraphErrors* histo, TH1D* syst);

void PlotPPb(
    TString sSpecies = "Pion",
    Int_t iCent = 0
)
{
    Bool_t bSyst = kFALSE;

    TString sPathOut = "./plots/";

    TString sFileSub = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/pPb/v22_subt_pPb_QM/v2-syst-final.root";
    TString sFileCums = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/pPb/v24/2464/V0A_020/gap_all/Graphs.root";
    TString sFileCumSyst = Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/pPb/v24/syst_2517/%s/final_diff/syst_total.root",sSpecies.Data());
    // TString sFileSub = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/pPb/cums/v22_subt_pPb_QM/v2-syst-final.root";
    // TString sFileCums = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/pPb/cums/2464/V0A_020/gap_all/Graphs.root";
    // TString sFileCumSyst = Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/pPb/cums/syst_2517/%s/final_diff/syst_total.root",sSpecies.Data());

    TString sHistoSub = Form("graphPoints_%s_cent%d", sSpecies.Data(), iCent);
    TString sHistoSubSyst = Form("graphSyst_%s_cent%d", sSpecies.Data(), iCent);
    TString sHistoCumTwo = Form("%s_hFlow2_harm2_gap-10_cent%d", sSpecies.Data(), iCent);
    TString sHistoCumFour = Form("%s_hFlow4_harm2_gap-10_cent%d", sSpecies.Data(), iCent);
    TString sHistoCumSyst = Form("%s_hFlow4_harm2_gap-10_cent%d", sSpecies.Data(), iCent);

    Color_t col[] = { kBlack, kRed, kBlue};
    Int_t markers[] = { kOpenSquare, kOpenCircle, kOpenCircle};

    /////////

    TFile* fileSub = TFile::Open(sFileSub.Data(), "READ");
    if(!fileSub) {  Utils::Error("fileSub not found!"); return; }
    TGraphErrors* grSub = (TGraphErrors*) fileSub->Get(sHistoSub.Data());
    if(!grSub) { Utils::Error(Form("grSub '%s' not found!",sFileSub.Data())); fileSub->ls(); return; }
    grSub->SetMarkerStyle(markers[0]);
    grSub->SetMarkerColor(col[0]);
    grSub->SetLineColor(col[0]);

    TGraphErrors* grSubSyst = (TGraphErrors*) fileSub->Get(sHistoSubSyst.Data());
    if(!grSubSyst) { Utils::Error(Form("grSubSyst '%s' not found!",sFileSub.Data())); fileSub->ls(); return; }
    grSubSyst->SetMarkerStyle(markers[0]);
    grSubSyst->SetMarkerColor(col[0]);
    grSubSyst->SetLineColor(col[0]);
    grSubSyst->SetFillColorAlpha(col[0],0.2);

    TFile* fileCums = TFile::Open(sFileCums.Data(), "READ");
    if(!fileCums) {  Utils::Error("fileCums not found!"); return; }

    TFile* fileCumsSyst = TFile::Open(sFileCumSyst.Data(), "READ");
    if(!fileCumsSyst) {  Utils::Error("sFileCumSyst not found!"); return; }

    TGraphErrors* grCumTwo = (TGraphErrors*) fileCums->Get(sHistoCumTwo);
    if(!grCumTwo) { Utils::Error(Form("grCumTwo '%s' not found!",sHistoCumTwo.Data())); fileCums->ls(); return; }
    grCumTwo->SetMarkerStyle(markers[1]);
    grCumTwo->SetMarkerColor(col[1]);
    grCumTwo->SetLineColor(col[1]);

    TGraphErrors* grCumFour = (TGraphErrors*) fileCums->Get(sHistoCumFour);
    if(!grCumFour) { Utils::Error(Form("grCumFour '%s' not found!",sHistoCumFour.Data())); fileCums->ls(); return; }
    grCumFour->SetMarkerStyle(markers[2]);
    grCumFour->SetMarkerColor(col[2]);
    grCumFour->SetLineColor(col[2]);

    TH1D* hCumSyst = (TH1D*) fileCumsSyst->Get(sHistoCumSyst);
    if(!hCumSyst) { Utils::Error(Form("hCumSyst '%s' not found!",sHistoCumSyst.Data())); fileCumsSyst->ls(); return; }
    TGraphErrors* grCumSyst = Transfer(grCumFour,hCumSyst);
    grCumSyst->SetMarkerStyle(markers[2]);
    grCumSyst->SetMarkerColor(col[2]);
    grCumSyst->SetLineColor(col[2]);
    grCumSyst->SetFillColorAlpha(col[2],0.2);

    TCanvas* can = new TCanvas("can","can", 600,600);
    can->cd();

    TLegend* leg = Utils::MakeLegend(Utils::kLegTopLeft);
    leg->SetHeader(Form("%s",sSpecies.Data()));
    leg->AddEntry(grCumTwo, "v_{2}{2}", "pl");
    leg->AddEntry(grSub, "v_{2}^{sub}{2,|#Delta#eta|>0.4}", "pl");
    leg->AddEntry(grCumFour, "v_{2}{4}", "pl");

    TH1* frame = (TH1*) gPad->DrawFrame(0.0,0.0,6.0,0.35);

    if(bSyst) grSubSyst->Draw("same p2");
    if(bSyst) grCumSyst->Draw("same p2");
    grSub->Draw("same pe");
    grCumTwo->Draw("same pe");
    grCumFour->Draw("same pe");
    leg->Draw();

    gSystem->mkdir(sPathOut.Data(),1);
    can->SaveAs(Form("%s/comp_%s_cent%d.pdf",sPathOut.Data(), sSpecies.Data(), iCent), "pdf");

    fileSub->Close();
    fileCums->Close();
}

void PlotAll(Int_t iCent = 0)
{
    PlotPPb("Charged",iCent);
    PlotPPb("Pion",iCent);
    PlotPPb("Kaon",iCent);
    PlotPPb("Proton",iCent);
    PlotPPb("K0s",iCent);
    PlotPPb("Lambda",iCent);
    // PlotPPb("Phi",iCent);
}

TGraphErrors* Transfer(TGraphErrors* points, TH1D* syst)
{
    if(!points) { printf("E::Transfer : Input points not found!\n"); return nullptr; }
    if(!syst) { printf("E::Transfer : Input syst histo not found!\n"); return nullptr; }

    Double_t dSystRel = 0.0;

    Double_t dX = 0.0;
    Double_t dY = 0.0;

    TGraphErrors* err = (TGraphErrors*) points->Clone("syst");

    // SetEx(err,0.05);
    Int_t np = err->GetN();
    for (Int_t i=0; i<np; i++) {
        err->GetPoint(i,dX,dY);
        dSystRel = syst->GetBinContent(i+1);
        err->SetPointError(i,0.05, dSystRel*dY);
    }

    return err;
}
