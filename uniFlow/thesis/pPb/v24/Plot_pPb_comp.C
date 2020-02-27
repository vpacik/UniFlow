// Macro for comparing p-Pb results of v2{4} with v2{2}^{sub} (Preliminaries)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStylePID.C"

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>
#include <TString.h>
#include <TFile.h>


TGraphErrors* Transfer(TGraphErrors* histo, TH1D* syst);

void Plot_pPb_comp(
    TString sSpecies = "Pion",
    Int_t iCent = 0
)
{
    Bool_t bSyst = kFALSE;

    TString sPathOut = "./plots/";
    TString sOutPlotName = "pPb_perSpecies_v2comp.pdf";

    TString sFileSub = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/pPb/v22_subt_pPb_QM/v2-syst-final.root";
    TString sFileCums = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/pPb/v24/2464/V0A_020/gap_all/Graphs.root";
    TString sFileCumSyst = Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/pPb/v24/syst_2517/%s/final_diff/syst_total.root",sSpecies.Data());
    // TString sFileSub = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/pPb/cums/v22_subt_pPb_QM/v2-syst-final.root";
    // TString sFileCums = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/pPb/cums/2464/V0A_020/gap_all/Graphs.root";
    // TString sFileCumSyst = Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/pPb/cums/syst_2517/%s/final_diff/syst_total.root",sSpecies.Data());

    Double_t dXmin = 0.0;
    Double_t dXmax = 6.0;
    Double_t dYmin = -0.1;
    Double_t dYmax = 0.4;

    TGaxis* gaX = new TGaxis(0.5, dYmin, 5.5, dYmin, 0.5, 5.5);
    gaX->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TGaxis* gaY = new TGaxis(dXmin, -0.05, dXmin, 0.35, -0.05,0.35);
    gaY->SetTitle("#it{v}_{2}{m}");
    // TString sCorrName = "v2{4}";

    TCanvas* can = new TCanvas("can","can",900,900);
    TList* frames = PrepareCanvas2x3(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);


    // can->SetLeftMargin(0.1);
    // gaY->SetTitleOffset(1.2);
    // gLatexPad->SetTextSize(22);

    can->cd(0);
    gLatexCan->DrawLatexNDC(0.13,0.97,Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

    const Int_t iNumPads = 6;
    TString sPadNames[iNumPads] = {"Charged", "Pion", "Kaon", "K0s", "Proton", "Lambda"};
    TString sPadLabels[iNumPads] = {"h^{#pm}","#pi^{#pm}","K^{#pm}","K_{S}^{0}","p/#bar{p}","#Lambda/#bar{#Lambda}"};

    // TString sPadNames[iNumPads] = {""};
    // TString sPadLabels[iNumPads] = {" 0-20% (V0A)"};

    // const Int_t iNumHistos = 3;
    // TString* sHistoNames = {sHistoCumTwo, sHistoCumFour, sHistoSub};
    // TString* sHistoLabels = {"a","b","c"};

    TString sDrawOpt = "pe";


    Color_t col[] = { kBlack, kBlue, kRed};
    Int_t markers[] = { kOpenSquare, kOpenCircle, kFullCircle};


    // Color_t cols[] = {kGray+2, kRed+1, kGreen+3, kGreen+1, kBlue, kOrange+2};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,kFullDiamond};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,47};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,43,kOpenSquare};
    // Double_t size[] = { 1.0,1.2,1.3,1.0,1.4,1.0};

    // gStyle->SetPalette(kRainBow);
    // TColor::InvertPalette();

    TLegend* leg1 = new TLegend(0.25,0.67,0.55,0.95);
    // leg1->SetHeader("");
    SetLegendStyle(leg1);

    TLegend* leg2 = new TLegend(0.22,0.12,0.45,0.33);
    // leg2->SetHeader("");
    SetLegendStyle(leg2);

    TFile* fileSub = TFile::Open(sFileSub.Data(), "READ");
    if(!fileSub) {  Utils::Error("fileSub not found!"); return; }
    // TGraphErrors* grSub = (TGraphErrors*) fileSub->Get(sHistoSub.Data());
    // if(!grSub) { Utils::Error(Form("grSub '%s' not found!",sFileSub.Data())); fileSub->ls(); return; }
    // grSub->SetMarkerStyle(markers[0]);
    // grSub->SetMarkerColor(col[0]);
    // grSub->SetLineColor(col[0]);
    //

    TFile* fileCums = TFile::Open(sFileCums.Data(), "READ");
    if(!fileCums) {  Utils::Error("fileCums not found!"); return; }

    // loop over pads
    for(Int_t p(0); p < iNumPads; ++p) {
        can->cd(p+1);

        TString sPadName = sPadNames[p];
        TString sPadLabel = sPadLabels[p];
        if(!sPadLabel.IsNull()) gLatexPad->DrawLatexNDC(0.94,0.96,sPadLabel.Data());

        Int_t iCent = 0;

        TString sHistoSub = Form("graphPoints_%s_cent%d", sPadName.Data(), iCent);
        TString sHistoSubSyst = Form("graphSyst_%s_cent%d", sPadName.Data(), iCent);
        TString sHistoCumTwo = Form("%s_hFlow2_harm2_gap-10_cent%d", sPadName.Data(), iCent);
        TString sHistoCumFour = Form("%s_hFlow4_harm2_gap-10_cent%d", sPadName.Data(), iCent);
        TString sHistoCumSyst = Form("%s_hFlow4_harm2_gap-10_cent%d", sPadName.Data(), iCent);

        TGraphErrors* grSub = (TGraphErrors*) fileSub->Get(sHistoSub.Data());
        if(!grSub) { Utils::Error(Form("grSub '%s' not found!",sHistoSub.Data())); fileSub->ls(); return; }

        TGraphErrors* grCumTwo = (TGraphErrors*) fileCums->Get(sHistoCumTwo);
        if(!grCumTwo) { Utils::Error(Form("grCumTwo '%s' not found!",sHistoCumTwo.Data())); fileCums->ls(); return; }

        TGraphErrors* grCumFour = (TGraphErrors*) fileCums->Get(sHistoCumFour);
        if(!grCumFour) { Utils::Error(Form("grCumFour '%s' not found!",sHistoCumFour.Data())); fileCums->ls(); return; }

        grSub->SetMarkerStyle(markers[0]);
        grSub->SetMarkerColor(col[0]);
        grSub->SetLineColor(col[0]);

        grCumTwo->SetMarkerStyle(markers[1]);
        grCumTwo->SetMarkerColor(col[1]);
        grCumTwo->SetLineColor(col[1]);

        grCumFour->SetMarkerStyle(markers[2]);
        grCumFour->SetMarkerColor(col[2]);
        grCumFour->SetLineColor(col[2]);

        grSub->Draw(sDrawOpt.Data());
        grCumFour->Draw(sDrawOpt.Data());
        grCumTwo->Draw(sDrawOpt.Data());

        if(p < 1) {
            leg1->AddEntry(grCumTwo,"#it{v}_{2}{2}","p");
            leg1->AddEntry(grSub,"#it{v}_{2}^{subt}{2,|#Delta#eta|>0.4}","p");
            leg1->AddEntry(grCumFour,"#it{v}_{2}{4}","p");
        }
    }


    can->cd(1);
    leg1->Draw();
    leg2->Draw();
    can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));
}
//
//
//
//     TFile* fileSub = TFile::Open(sFileSub.Data(), "READ");
//     if(!fileSub) {  Utils::Error("fileSub not found!"); return; }
//     TGraphErrors* grSub = (TGraphErrors*) fileSub->Get(sHistoSub.Data());
//     if(!grSub) { Utils::Error(Form("grSub '%s' not found!",sFileSub.Data())); fileSub->ls(); return; }
//     grSub->SetMarkerStyle(markers[0]);
//     grSub->SetMarkerColor(col[0]);
//     grSub->SetLineColor(col[0]);
//
//     TGraphErrors* grSubSyst = (TGraphErrors*) fileSub->Get(sHistoSubSyst.Data());
//     if(!grSubSyst) { Utils::Error(Form("grSubSyst '%s' not found!",sFileSub.Data())); fileSub->ls(); return; }
//     grSubSyst->SetMarkerStyle(markers[0]);
//     grSubSyst->SetMarkerColor(col[0]);
//     grSubSyst->SetLineColor(col[0]);
//     grSubSyst->SetFillColorAlpha(col[0],0.2);
//
//     TFile* fileCums = TFile::Open(sFileCums.Data(), "READ");
//     if(!fileCums) {  Utils::Error("fileCums not found!"); return; }
//
//     TFile* fileCumsSyst = TFile::Open(sFileCumSyst.Data(), "READ");
//     if(!fileCumsSyst) {  Utils::Error("sFileCumSyst not found!"); return; }
//
//     TGraphErrors* grCumTwo = (TGraphErrors*) fileCums->Get(sHistoCumTwo);
//     if(!grCumTwo) { Utils::Error(Form("grCumTwo '%s' not found!",sHistoCumTwo.Data())); fileCums->ls(); return; }
//     grCumTwo->SetMarkerStyle(markers[1]);
//     grCumTwo->SetMarkerColor(col[1]);
//     grCumTwo->SetLineColor(col[1]);
//
//     TGraphErrors* grCumFour = (TGraphErrors*) fileCums->Get(sHistoCumFour);
//     if(!grCumFour) { Utils::Error(Form("grCumFour '%s' not found!",sHistoCumFour.Data())); fileCums->ls(); return; }
//     grCumFour->SetMarkerStyle(markers[2]);
//     grCumFour->SetMarkerColor(col[2]);
//     grCumFour->SetLineColor(col[2]);
//
//     TH1D* hCumSyst = (TH1D*) fileCumsSyst->Get(sHistoCumSyst);
//     if(!hCumSyst) { Utils::Error(Form("hCumSyst '%s' not found!",sHistoCumSyst.Data())); fileCumsSyst->ls(); return; }
//     TGraphErrors* grCumSyst = Transfer(grCumFour,hCumSyst);
//     grCumSyst->SetMarkerStyle(markers[2]);
//     grCumSyst->SetMarkerColor(col[2]);
//     grCumSyst->SetLineColor(col[2]);
//     grCumSyst->SetFillColorAlpha(col[2],0.2);
//
//     can->cd();
//
//     TLegend* leg = Utils::MakeLegend(Utils::kLegTopLeft);
//     leg->SetHeader(Form("%s",sSpecies.Data()));
//     leg->AddEntry(grCumTwo, "v_{2}{2}", "pl");
//     leg->AddEntry(grSub, "v_{2}^{sub}{2,|#Delta#eta|>0.4}", "pl");
//     leg->AddEntry(grCumFour, "v_{2}{4}", "pl");
//
//     if(bSyst) grSubSyst->Draw("same p2");
//     if(bSyst) grCumSyst->Draw("same p2");
//     grSub->Draw("same pe");
//     grCumTwo->Draw("same pe");
//     grCumFour->Draw("same pe");
//     leg->Draw();
//
//     gSystem->mkdir(sPathOut.Data(),1);
//     can->SaveAs(Form("%s/comp_%s_cent%d.pdf",sPathOut.Data(), sSpecies.Data(), iCent), "pdf");
//
//     fileSub->Close();
//     fileCums->Close();
// }
//
// // void Plot_pPb_comp(Int_t iCent = 0)
// // {
// //     PlotPPb("Charged",iCent);
// //     PlotPPb("Pion",iCent);
// //     PlotPPb("Kaon",iCent);
// //     PlotPPb("Proton",iCent);
// //     PlotPPb("K0s",iCent);
// //     PlotPPb("Lambda",iCent);
// //     // PlotPPb("Phi",iCent);
// // }
//
// TGraphErrors* Transfer(TGraphErrors* points, TH1D* syst)
// {
//     if(!points) { printf("E::Transfer : Input points not found!\n"); return nullptr; }
//     if(!syst) { printf("E::Transfer : Input syst histo not found!\n"); return nullptr; }
//
//     Double_t dSystRel = 0.0;
//
//     Double_t dX = 0.0;
//     Double_t dY = 0.0;
//
//     TGraphErrors* err = (TGraphErrors*) points->Clone("syst");
//
//     // SetEx(err,0.05);
//     Int_t np = err->GetN();
//     for (Int_t i=0; i<np; i++) {
//         err->GetPoint(i,dX,dY);
//         dSystRel = syst->GetBinContent(i+1);
//         err->SetPointError(i,0.05, dSystRel*dY);
//     }
//
//     return err;
// }
