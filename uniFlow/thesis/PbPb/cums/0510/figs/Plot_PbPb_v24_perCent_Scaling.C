// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libs/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libs/libCanvas.C"

void Plot_PbPb_v24_perCent_Scaling()
{
    // === Preparing canNCQvas =================

    Double_t dXmin = 0.0;
    Double_t dXmax = 5.0;
    Double_t dYmin = -0.01;
    Double_t dYmax = 0.12;

    TString sHistoName = "v2{4}";

    TCanvas* canNCQ = new TCanvas("canNCQ","canNCQ",900,600);
    TCanvas* canKEt = new TCanvas("canKEt","canKEt",900,600);

    TGaxis* gaY = new TGaxis(dXmin, 0.0, dXmin, 0.11, 0.0,0.11);
    gaY->SetTitle("v_{2}{4}/n_{q}");

    TGaxis* gaXNCQ = new TGaxis(0.5, dYmin, 4.5, dYmin, 0.5, 4.5);
    gaXNCQ->SetTitle("#it{p}_{T}/n_{q} (GeV/#it{c})");
    TGaxis* gaXKEt = new TGaxis(0.5, dYmin, 4.5, dYmin, 0.5, 4.5);
    gaXKEt->SetTitle("(#it{m}_{T} - #it{m}_{0})/n_{q} (GeV/#it{c}^{2})");

    TList* framesNCQ = PrepareCanvas3x2(canNCQ, gaXNCQ, gaY, dXmin, dYmin, dXmax, dYmax);
    TList* framesKEt = PrepareCanvas3x2(canKEt, gaXKEt, gaY, dXmin, dYmin, dXmax, dYmax);

    // Drawing lines at zero
    TLine lineZero = TLine();
    lineZero.SetLineColor(kGray+1);
    lineZero.SetLineStyle(kDashed);

    for(Int_t i(1); i < 7; ++i) {
        canNCQ->cd(i);
        lineZero.DrawLine(dXmin, 0.0, dXmax, 0.0);
        canKEt->cd(i);
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

    TLegend* leg1 = new TLegend(0.27,0.7,0.52,0.96);
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

    TFile* fileIn = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/results/PbPb/cums/v24_prel/graphs_PbPb_v22_v24_scaled.root","READ");
    if(!fileIn) { printf("Input file not openned!\n"); return; }

    // TString sSpecies[iNumPads] = {"Charged", "Pion", "Kaon", "K0s", "Proton", "Lambda", "Phi"};

    const Int_t iNumCent = 6;
    TString sNames[iNumCent] = {"Pion","Kaon","K0s","Proton","Lambda","Phi"};
    TString sLabels[iNumCent] = {"#pi^{#pm}","K^{#pm}","K^{0}_{S}","p/#bar{p}","#Lambda/#bar{#Lambda}","#phi"};

    Color_t cols[iNumCent] = { kRed+1, kGreen+3, kGreen+1, kBlue, kOrange+2, kMagenta};
    Int_t mark[iNumCent] = { kFullCircle,kFullTriangleUp,kFullTriangleDown,kFullSquare,kFullDiamond, kFullCross };
    Double_t size[iNumCent] = { 0.8,1.0,1.0,0.7,1.0 };
    gStyle->SetPalette(kRainBow);
    // TColor::InvertPalette();

    TGraphErrors* gr = nullptr;

    for(Int_t i(0); i < iNumCent; ++i) {

        // NCQ
        canNCQ->cd(1);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_NCQ",sNames[i].Data(),sHistoName.Data(),"0510"));
        if(!gr) { fileIn->ls(); printf("%s_%s_cent%s_KEt\n",sNames[i].Data(),sHistoName.Data(),"0510"); return; }
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"5-10%");
        if(i < 3) leg1->AddEntry(gr,sLabels[i].Data(),"pel");
        else leg2->AddEntry(gr,sLabels[i].Data(),"pel");

        canNCQ->cd(2);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_NCQ",sNames[i].Data(),sHistoName.Data(),"1020"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"10-20%");

        canNCQ->cd(3);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_NCQ",sNames[i].Data(),sHistoName.Data(),"2030"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"20-30%");

        canNCQ->cd(4);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_NCQ",sNames[i].Data(),sHistoName.Data(),"3040"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"30-40%");

        canNCQ->cd(5);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_NCQ",sNames[i].Data(),sHistoName.Data(),"4050"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"40-50%");

        canNCQ->cd(6);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_NCQ",sNames[i].Data(),sHistoName.Data(),"5060"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"50-60%");

        // KEt

        canKEt->cd(1);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_KEt",sNames[i].Data(),sHistoName.Data(),"0510"));
        if(!gr) { fileIn->ls(); printf("%s_%s_cent%s_KEt\n",sNames[i].Data(),sHistoName.Data(),"0510"); return; }
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"5-10%");
        // leg->AddEntry(gr,sLabels[i].Data(),"pel");

        canKEt->cd(2);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_KEt",sNames[i].Data(),sHistoName.Data(),"1020"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"10-20%");

        canKEt->cd(3);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_KEt",sNames[i].Data(),sHistoName.Data(),"2030"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"20-30%");

        canKEt->cd(4);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_KEt",sNames[i].Data(),sHistoName.Data(),"3040"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"30-40%");

        canKEt->cd(5);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_KEt",sNames[i].Data(),sHistoName.Data(),"4050"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"40-50%");

        canKEt->cd(6);
        gr = (TGraphErrors*) fileIn->Get(Form("%s_%s_cent%s_KEt",sNames[i].Data(),sHistoName.Data(),"5060"));
        // SetColor(gr,cols[i]);
        SetMarker(gr,mark[i],size[i]);
        gr->Draw("pe");
        textPad->DrawLatexNDC(0.94,0.94,"50-60%");
    }

    canNCQ->cd();
    textCan->DrawLatexNDC(0.09,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8"));

    canKEt->cd();
    textCan->DrawLatexNDC(0.09,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8"));

    canNCQ->cd(1);
    leg1->Draw();
    leg2->Draw();
    canNCQ->SaveAs("./plots/PbPb_v24_perCent_NCQ.pdf","pdf");

    canKEt->cd(1);
    leg1->Draw();
    leg2->Draw();
    canKEt->SaveAs("./plots/PbPb_v24_perCent_KEt.pdf","pdf");

}
