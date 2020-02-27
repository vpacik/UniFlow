// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libs/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libs/libCanvas.C"

void Plot_PbPb_v22gap08_perCent()
{
    // === Preparing canvas =================

    Double_t dXmin = 0.0;
    Double_t dXmax = 10.0;
    Double_t dYmin = -0.03;
    Double_t dYmax = 0.35;

    TCanvas* can = new TCanvas("can","can",900,600);

    TString sHistoName = "v2{2,gap08}";

    TGaxis* gaX = new TGaxis(1.0, dYmin, 9.0, dYmin, 1.0, 9.0);
    gaX->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TGaxis* gaY = new TGaxis(dXmin, 0.0, dXmin, 0.3, 0.0,0.3);
    gaY->SetTitle("v_{2}{2,|#Delta#eta|>0.8}");

    TList* frames = PrepareCanvas3x2(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    // Drawing lines at zero
    TLine lineZero = TLine();
    lineZero.SetLineColor(kGray+1);
    lineZero.SetLineStyle(kDashed);

    for(Int_t i(1); i < 7; ++i) {
        can->cd(i);
        lineZero.DrawLine(dXmin, 0.0, dXmax, 0.0);
    }

    // === Filling head =================

    TLatex * textCan = new TLatex();
    textCan->SetTextFont(43);
    textCan->SetTextSize(20);
    textCan->SetTextAlign(11);

    TLatex* textPad = new TLatex();
    textPad->SetTextFont(43);
    textPad->SetTextSize(20);
    textPad->SetTextAlign(33);

    TLegend* leg1 = new TLegend(0.27,0.61,0.52,0.96);
    leg1->SetFillColorAlpha(0,0);
    leg1->SetBorderSize(0);
    leg1->SetTextFont(43);
    leg1->SetTextSize(14);

    TLegend* leg2 = new TLegend(0.4,0.7,0.65,0.96);
    leg2->SetFillColorAlpha(0,0);
    leg2->SetBorderSize(0);
    leg2->SetTextFont(43);
    leg2->SetTextSize(14);

    // === Loading datapoints =================

    TFile* fileIn = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/results/PbPb/cums/v24_prel/graphs_PbPb_v22_v24.root","READ");
    if(!fileIn) { printf("Input file not openned!\n"); return; }

    // TString sSpecies[iNumPads] = {"Charged", "Pion", "Kaon", "K0s", "Proton", "Lambda", "Phi"};

    const Int_t iNumCent = 7;
    TString sNames[iNumCent] = {"Charged","Pion","Kaon","K0s","Proton","Lambda","Phi"};
    TString sLabels[iNumCent] = {"h^{#pm}","#pi^{#pm}","K^{#pm}","K^{0}_{S}","p/#bar{p}","#Lambda/#bar{#Lambda}","#phi"};

    Color_t cols[iNumCent] = {kGray+2, kRed+1, kGreen+3, kGreen+1, kBlue, kOrange+2, kMagenta};
    Int_t mark[iNumCent] = { kOpenCircle,kFullCircle,kFullTriangleUp,kFullTriangleDown,kFullSquare,kFullDiamond, kFullCross };
    Double_t size[iNumCent] = { 0.8,0.8,1.0,1.0,0.7,1.0 };
    gStyle->SetPalette(kRainBow);
    // TColor::InvertPalette();

    TGraphErrors* gr = nullptr;

    for(Int_t i(0); i < iNumCent; ++i) {
        can->cd(1);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s",sNames[i].Data(),sHistoName.Data(),"0510"));
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"5-10%");
        if(i < 4) leg1->AddEntry(gr,sLabels[i].Data(),"pl");
        else leg2->AddEntry(gr,sLabels[i].Data(),"pl");

        can->cd(2);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s",sNames[i].Data(),sHistoName.Data(),"1020"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"10-20%");

        can->cd(3);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s",sNames[i].Data(),sHistoName.Data(),"2030"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"20-30%");

        can->cd(4);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s",sNames[i].Data(),sHistoName.Data(),"3040"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"30-40%");

        can->cd(5);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s",sNames[i].Data(),sHistoName.Data(),"4050"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"40-50%");

        can->cd(6);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s",sNames[i].Data(),sHistoName.Data(),"5060"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"50-60%");
    }

    can->cd();
    textCan->DrawLatexNDC(0.09,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8"));

    can->cd(1);
    leg1->Draw();
    leg2->Draw();
    can->SaveAs("./plots/PbPb_v22gap08_perCent.pdf","pdf");
}
