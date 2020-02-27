
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStylePID.C"


void Plot_pPb_Ref_v24()
{
    // === Preparing canvas =================


    TString sInputFileName = "pPb/v24/2464/V0A_norebin/gap_all_refs/Processed.root";
    TString sInputFileNameMerge = "pPb/v24/2464/V0A_020/gap_all/Processed.root";
    TString sOutPlotName = "pPb_Refs_c24.pdf";
    // TString sOutPlotForm = "PbPb_v24_perCent";

    Double_t dXmin = 0.0;
    Double_t dXmax = 60.0;
    Double_t dYmin = -0.00002;
    Double_t dYmax = 0.00008;

    // TGaxis* gaX = new TGaxis(dXmin, dYmin, dXmax, dYmin, dXmin, dXmax);
    TGaxis* gaX = new TGaxis(0.0, dYmin, 60.0, dYmin, 0., 60.0);
    gaX->SetTitle("V0A (%)");
    TGaxis* gaY = new TGaxis(dXmin, dYmin, dXmin, dYmax, dYmin, dYmax);
    // TGaxis* gaY = new TGaxis(dXmin, -0.1, dXmin, 0.25, -0.1,0.25);
    gaY->SetTitle("#it{c}_{2}{4}");
    TString sCorrName = "hCum4_harm2_gap-10";

    TCanvas* can = new TCanvas("can","can",600,600);
    TList* frames = PrepareCanvas1x1(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);


    // can->SetLeftMargin(0.14);
    // gaY->SetTitleOffset(1.5);
    gLatexPad->SetTextSize(22);

    // can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));

    // return;

    can->cd();
    gLatexCan->DrawLatexNDC(0.25,0.955,Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

// return;
    const Int_t iNumPads = 1;
    TString sPadNames[iNumPads] = {""};
    TString sPadLabels[iNumPads] = {""};

    SetStylePID();

    const Int_t iNumHistos = 1;
    TString sHistoNames[iNumHistos] = {"Refs"};
    TString sHistoLabels[iNumHistos] = {"#it{c}_{2}{4}"};

    TString sDrawOpt = "pe";

    // Color_t cols[] = {kGray+2, kRed+1, kGreen+3, kGreen+1, kBlue, kOrange+2};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,kFullDiamond};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,47};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,43,kOpenSquare};
    // Double_t size[] = { 1.0,1.2,1.3,1.0,1.4,1.0};

    // gStyle->SetPalette(kRainBow);
    // TColor::InvertPalette();

    TLegend* leg1 = new TLegend(0.12,0.12,0.35,0.33);
    // leg1->SetHeader("");
    SetLegendStyle(leg1);

    TLegend* leg2 = new TLegend(0.22,0.12,0.45,0.33);
    // leg2->SetHeader("");
    SetLegendStyle(leg2);

    // === Loading datapoints =================


    TFile* fileIn = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileName.Data()),"READ");
    if(!fileIn) { printf("Input file '%s' not openned!\n",sInputFileName.Data()); return; }

    TFile* fileInMerge = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileNameMerge.Data()),"READ");
    if(!fileInMerge) { printf("Input file '%s' not openned!\n",sInputFileNameMerge.Data()); return; }

    TGraphErrors* gr = nullptr;

    // loop over pads
    for(Int_t p(0); p < iNumPads; ++p) {
        can->cd(p+1);

        TString sPadName = sPadNames[p];
        TString sPadLabel = sPadLabels[p];
        if(!sPadLabel.IsNull()) gLatexPad->DrawLatexNDC(0.94,0.9,sPadLabel.Data());

        TString sGraphName = Form("%s_%s%s", "Refs", sCorrName.Data(), sPadName.Data());
        printf("%s\n",sGraphName.Data());

        TGraphErrors* gr = LoadHistoToGraph(fileIn,sGraphName);
        if(!gr) {
            printf("Graph '%s' not found!\n",sGraphName.Data());
            // fileIn->ls();
            return;
        }

        TGraphErrors* grMerge = LoadHistoToGraph(fileInMerge,sGraphName);
        if(!grMerge) {
            printf("Graph '%s' not found!\n",sGraphName.Data());
            // fileIn->ls();
            return;
        }


        SetGraphError(gr,0.0);

        SetMarkerPID(gr,kSpecies(kRefs));
        gr->Draw(sDrawOpt.Data());
        // grMerge->Draw(sDrawOpt.Data());

        if(p < 1) {
            // leg1->AddEntry(gr,"C","p");
        }

    }

    can->cd(4);
    leg1->Draw();
    leg2->Draw();
    can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));
}
