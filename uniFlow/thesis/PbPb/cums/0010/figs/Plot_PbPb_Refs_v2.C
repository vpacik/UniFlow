// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"

void Plot_PbPb_Refs_v2()
{
    // === Preparing canvas =================

    TString sInputFileName = "PbPb/cums/0010/v2/graphs_PbPb_v22_v24.root";
    TString sOutPlotName =  "PbPb_Refs_v2.pdf";

    Double_t dXmin = 0.0;
    Double_t dXmax = 60.0;
    Double_t dYmin = 0.0;
    Double_t dYmax = 0.12;

    TGaxis* gaX = new TGaxis(dXmin, dYmin, dXmax, dYmin, dXmin, dXmax);
    gaX->SetTitle("V0M (%)");
    TGaxis* gaY = new TGaxis(dXmin, dYmin, dXmin, dYmax, dYmin,dYmax);
    gaY->SetTitle("#it{v}_{2}{m}");

    TCanvas* can = new TCanvas("can","can",600,600);
    TList* frames = PrepareCanvas1x1(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    can->cd();
    gLatexCan->DrawLatexNDC(0.14,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

    const Int_t iNumPads = 1;
    TString sPadNames[iNumPads] = {""};
    TString sPadLabels[iNumPads] = {""};

    const Int_t iNumHistos = 5;
    TString sHistoNames[iNumHistos] = {"Refs_v2{2}","Refs_v2{2,gap00}", "Refs_v2{2,gap08}","Refs_v2{4}","Refs_v2{4,gap00}"};
    TString sHistoLabels[iNumHistos] = {"#it{v}_{2}{2}","#it{v}_{2}{2,|#Delta#eta| > 0}","#it{v}_{2}{2,|#Delta#eta| > 0.8}","#it{v}_{2}{4}","#it{v}_{2}{4,|#Delta#eta| > 0}"};

    TString sDrawOpt = "pe PLC PMC";

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
    TLegend* leg1 = new TLegend(0.6,0.14,0.9,0.4);
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

    TGraphErrors* gr = nullptr;

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

            gr = (TGraphErrors*) fileIn->Get(sGraphName.Data());
            if(!gr) {
                printf("Graph '%s' not found!\n",sGraphName.Data());
                // fileIn->ls();
                return;
            }

            SetMarker(gr,mark[h],size[h]);
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
