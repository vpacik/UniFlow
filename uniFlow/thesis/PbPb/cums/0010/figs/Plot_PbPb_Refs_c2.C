// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"

void Plot_PbPb_Refs_c2()
{
    // === Preparing canvas =================

    TString sInputFileName = "PbPb/cums/0010/v2/process/Processed.root";
    TString sOutPlotName =  "PbPb_Refs_c2.pdf";

    Double_t dXmin = 0.0;
    Double_t dXmax = 60.0;
    Double_t dYmin = -0.007;
    Double_t dYmax = 0.013;

    TGaxis* gaX = new TGaxis(dXmin, dYmin, dXmax, dYmin, dXmin, dXmax);
    gaX->SetTitle("V0M (%)");
    TGaxis* gaY = new TGaxis(dXmin, dYmin, dXmin, dYmax, dYmin,dYmax);
    gaY->SetTitle("#it{c}_{2}{m}");
    gaY->SetMaxDigits(2);

    TCanvas* can = new TCanvas("can","can",600,600);
    TList* frames = PrepareCanvas1x1(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    can->cd();
    gLatexCan->DrawLatexNDC(0.23,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

    const Int_t iNumPads = 1;
    TString sPadNames[iNumPads] = {""};
    TString sPadLabels[iNumPads] = {""};

    const Int_t iNumHistos = 5;
    TString sHistoNames[iNumHistos] = {"Refs_hCum2_harm2_gap-10","Refs_hCum2_harm2_gap00", "Refs_hCum2_harm2_gap08","Refs_hCum4_harm2_gap-10","Refs_hCum4_harm2_gap00"};
    TString sHistoLabels[iNumHistos] = {"#it{c}_{2}{2}","#it{c}_{2}{2,|#Delta#eta| > 0}","#it{c}_{2}{2,|#Delta#eta| > 0.8}","#it{c}_{2}{4} (#times10^{2})","#it{c}_{2}{4,|#Delta#eta| > 0} (#times10^{2})"};

    TString sDrawOpt = "pe same PMC PLC";

    Color_t cols[] = {kGray+2, kRed+1, kGreen+3, kGreen+1, kBlue, kOrange+2};
    Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,kFullDiamond};
    Double_t size[] = { 1.0,1.2,1.3,1.0,1.0,1.0};

    // Color_t cols[] = {kGreen+2, kRed, kBlue-3, kViolet-3, kOrange+2};
    // Int_t mark[] = { kFullSquare,kFullCircle,kFullDiamond,kOpenSquare,kOpenCircle,kOpenCircle };
    // // Int_t mark[] = { kOpenSquare,kOpenCircle,kOpenDiamond,kFullSquare,kFullCircle,kOpenCircle };
    // Double_t size[] = { 1.0,1.0,1.4,1.0,1.0,1.0 };

    gStyle->SetPalette(kRainBow);
    TColor::InvertPalette();

    // TLegend* leg1 = new TLegend(0.27,0.5,0.52,0.94);
    TLegend* leg1 = new TLegend(0.6,0.44,0.9,0.7);
    leg1->SetHeader("RFPs (0.2 < #it{p}_{T} < 3 GeV/#it{c})");
    SetLegendStyle(leg1);



 ////// old -down

    // Double_t dXmin = 0.0;
    // Double_t dXmax = 60.0;
    // Double_t dYmin = 0.0;
    // Double_t dYmax = 0.12;
    //
    // TString sOutPlotName = "PbPb_Refs_v2_etaGap";
    //
    // TGaxis* gaX = new TGaxis(dXmin, dYmin, dXmax, dYmin, dXmin, dXmax);
    // gaX->SetTitle("V0M (%)");
    // TGaxis* gaY = new TGaxis(dXmin, dYmin, dXmin, dYmax, dYmin,dYmax);
    // gaY->SetTitle("#it{v}_{2}{m}");q
    //
    // TCanvas* can = new TCanvas("can","can",600,600);
    // PrepareCanvas1x1(can,gaX,gaY,dXmin,dYmin,dXmax,dYmax);
    //
    // const Int_t iNumNames = 5;
    // TString sNames[iNumNames] = {"Refs_v2{2}","Refs_v2{2,gap00}", "Refs_v2{2,gap08}","Refs_v2{4}","Refs_v2{4,gap00}"};
    // TString sLabels[iNumNames] = {"#it{v}_{2}{2}","#it{v}_{2}{2,|#Delta#eta| > 0}","#it{v}_{2}{2,|#Delta#eta| > 0.8}","#it{v}_{2}{4}","#it{v}_{2}{4,|#Delta#eta| > 0}"};
    //
    // Color_t cols[] = {kGreen+2, kRed, kBlue-3, kViolet-3, kOrange+2};
    // Int_t mark[] = { kFullSquare,kFullCircle,kFullDiamond,kOpenSquare,kOpenCircle,kOpenCircle };
    // // Int_t mark[] = { kOpenSquare,kOpenCircle,kOpenDiamond,kFullSquare,kFullCircle,kOpenCircle };
    // Double_t size[] = { 1.0,1.0,1.4,1.0,1.0,1.0 };
    //
    // // === Filling head =================
    //
    // TLatex * textCan = new TLatex();
    // textCan->SetTextFont(43);
    // textCan->SetTextSize(20);
    // textCan->SetTextAlign(11);
    //
    // TLatex* textPad = new TLatex();
    // textPad->SetTextFont(43);
    // textPad->SetTextSize(20);
    // textPad->SetTextAlign(33);
    //
    // TLegend* leg = new TLegend(0.6,0.14,0.9,0.4);
    // leg->SetFillColorAlpha(0,0);
    // leg->SetBorderSize(0);
    // leg->SetTextFont(43);
    // leg->SetTextSize(20);
    //
    // // leg->SetHeader("Event class (V0M)");
    // // TLegendEntry *header = (TLegendEntry*) leg->GetListOfPrimitives()->First();
    // // header->SetTextFont(43);
    // // header->SetTextSize(20);

    // === Loading datapoints =================

    TFile* fileIn = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileName.Data()),"READ");
    if(!fileIn) { printf("Input file '%s' not openned!\n",sInputFileName.Data()); return; }

    TH1* gr = nullptr;

    // loop over pads
    for(Int_t p(0); p < iNumPads; ++p) {
        can->cd(p+1);

        TString sPadName = sPadNames[p];
        TString sPadLabel = sPadLabels[p];
        if(!sPadLabel.IsNull()) gLatexPad->DrawLatexNDC(0.94,0.94,sPadLabel.Data());


        // loop over histos within given pad
        for(Int_t h(0); h < iNumHistos; ++h) {

            TString sHistoName = sHistoNames[h];
            TString sHistoLabel = sHistoLabels[h];

            TString sGraphName = Form("%s", sHistoName.Data());
            // printf("%s\n",sGraphName.Data());

            gr = (TH1D*) fileIn->Get(sGraphName.Data());
            if(!gr) {
                printf("Graph '%s' not found!\n",sGraphName.Data());
                // fileIn->ls();
                return;
            }

            if(h > 2) {
                gr->Scale(100);

            }

            // SetMarker(gr,mark[h],size[h]);
            gr->SetMarkerStyle(mark[h]);
            gr->SetMarkerSize(size[h]);
            // gr->SetMarkerColor(size[h]);
            gr->Draw(sDrawOpt.Data());

            if(p < 1 && !sHistoLabel.IsNull()) {
                leg1->AddEntry(gr,sHistoLabel.Data(),"pl");
            }
        }
    }

    can->cd(1);
    leg1->Draw();
    can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));
}
