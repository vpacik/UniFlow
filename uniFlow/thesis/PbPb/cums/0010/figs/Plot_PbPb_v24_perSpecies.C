// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libs/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libs/libCanvas.C"

void Plot_PbPb_v24_perSpecies()
{
    // === Preparing canvas =================

    Double_t dXmin = 0.0;
    Double_t dXmax = 10.0;
    Double_t dYmin = -0.02;
    Double_t dYmax = 0.28;

    TCanvas* can = new TCanvas("can","can",900,900);

    TGaxis* gaX = new TGaxis(1.0, dYmin, 9.0, dYmin, 1.0, 9.0);
    gaX->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TGaxis* gaY = new TGaxis(dXmin, 0.0, dXmin, 0.25, 0.0,0.25);
    gaY->SetTitle("v_{2}{4}");

    TString sHistoName = "v2{4}";

    TList* frames = PrepareCanvas3x3(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    // Drawing lines at zero
    TLine lineZero = TLine();
    lineZero.SetLineColor(kGray+1);
    lineZero.SetLineStyle(kDashed);

    for(Int_t i(1); i < 10; ++i) {
        can->cd(i);
        lineZero.DrawLine(dXmin, 0.0, dXmax, 0.0);
    }

    // Deleting TPads and frames
    delete frames->At(1);
    delete frames->At(2);

    can->cd(2);
    delete gPad;

    can->cd(3);
    delete gPad;

    // === Filling head =================

    TLatex * textCan = new TLatex();
    textCan->SetTextFont(43);
    textCan->SetTextSize(20);
    textCan->SetTextAlign(11);

    TLatex* textPad = new TLatex();
    textPad->SetTextFont(43);
    textPad->SetTextSize(20);
    textPad->SetTextAlign(33);

    TLegend* leg = new TLegend(0.4,0.7,0.6,0.95);
    leg->SetFillColorAlpha(0,0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(20);

    leg->SetHeader("Event class (V0M)");
    TLegendEntry *header = (TLegendEntry*) leg->GetListOfPrimitives()->First();
    header->SetTextFont(43);
    header->SetTextSize(20);

    // === Loading datapoints =================

    TFile* fileIn = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/results/PbPb/cums/v24_prel/graphs_PbPb_v22_v24.root","READ");
    if(!fileIn) { printf("Input file not openned!\n"); return; }

    const Int_t iNumCent = 6;
    TString sNames[iNumCent] = {"0510","1020","2030","3040","4050","5060"};
    TString sLabels[iNumCent] = {"  5-10%","10-20%","20-30%","30-40%","40-50%","50-60%"};

    Color_t cols[iNumCent] = {kBlack, kRed, kBlue, kGreen, kOrange, kMagenta};
    Int_t mark[iNumCent] = { kFullCross,kFullDiamond,kFullSquare,kFullCircle,kOpenSquare,kOpenCircle };
    Double_t size[iNumCent] = { 1.2,1.4,1.0,1.0,1.0,1.0 };
    gStyle->SetPalette(kRainBow);
    // TColor::InvertPalette();

    TGraphErrors* gr = nullptr;

    for(Int_t i(0); i < iNumCent; ++i) {
        can->cd(1);
        gr = (TGraphErrors*) fileIn->Get(Form("Charged_%s_cent%s",sHistoName.Data(),sNames[i].Data()));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe PLC PMC");
        textPad->DrawLatexNDC(0.94,0.94,"h^{#pm}");

        can->cd(4);
        gr = (TGraphErrors*) fileIn->Get(Form("Pion_%s_cent%s",sHistoName.Data(),sNames[i].Data()));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe PLC PMC");
        leg->AddEntry(gr,sLabels[i],"pl");
        textPad->DrawLatexNDC(0.94,0.94,"#pi^{#pm}");

        can->cd(5);
        gr = (TGraphErrors*) fileIn->Get(Form("Kaon_%s_cent%s",sHistoName.Data(),sNames[i].Data()));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe PLC PMC");
        textPad->DrawLatexNDC(0.94,0.94,"K^{#pm}");

        can->cd(6);
        gr = (TGraphErrors*) fileIn->Get(Form("K0s_%s_cent%s",sHistoName.Data(),sNames[i].Data()));
        // SetColor(gr,cols[i])
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe PLC PMC");
        textPad->DrawLatexNDC(0.94,0.94,"K^{0}_{S}");

        can->cd(7);
        gr = (TGraphErrors*) fileIn->Get(Form("Proton_%s_cent%s",sHistoName.Data(),sNames[i].Data()));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe PLC PMC");
        textPad->DrawLatexNDC(0.94,0.94,"p/#bar{p}");

        can->cd(8);
        gr = (TGraphErrors*) fileIn->Get(Form("Lambda_%s_cent%s",sHistoName.Data(),sNames[i].Data()));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe PLC PMC");
        textPad->DrawLatexNDC(0.92,0.92,"#Lambda/#bar{#Lambda}");

        can->cd(9);
        gr = (TGraphErrors*) fileIn->Get(Form("Phi_%s_cent%s",sHistoName.Data(),sNames[i].Data()));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe PLC PMC");
        textPad->DrawLatexNDC(0.92,0.92,"#phi");
    }

    can->cd();
    textCan->DrawLatexNDC(0.085,0.968,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8"));
    leg->Draw();

    can->SaveAs("./plots/PbPb_v24_perSpecies.pdf","pdf");

}
