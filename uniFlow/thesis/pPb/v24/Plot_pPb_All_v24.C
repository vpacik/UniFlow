
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStylePID.C"


void Plot_pPb_All_v24()
{
    // === Preparing canvas =================

    TString sInputFileNameModel = "pPb/v24/2464/VOA_020/gap_all/Graphs.root";
    TString sInputFileNameTRENTo = "PbPb/models/graphs_PbPb_Fv2_TRENTo.root";

    TString sInputFileName = "pPb/v24/2464/V0A_020/gap_all/Graphs.root";
    TString sOutPlotName = "pPb_All_v24.pdf";
    // TString sOutPlotForm = "PbPb_v24_perCent";

    Double_t dXmin = 0.0;
    Double_t dXmax = 6.0;
    Double_t dYmin = -0.1;
    Double_t dYmax = 0.25;

    TGaxis* gaX = new TGaxis(0.0, dYmin, 6.0, dYmin, 0., 6.0);
    gaX->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TGaxis* gaY = new TGaxis(dXmin, -0.1, dXmin, 0.25, -0.1,0.25);
    gaY->SetTitle("#it{v}_{2}{4}");
    TString sCorrName = "hFlow4_harm2_gap-10";

    TCanvas* can = new TCanvas("can","can",600,600);
    TList* frames = PrepareCanvas1x1(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    // can->SetLeftMargin(0.14);
    // gaY->SetTitleOffset(1.5);
    gLatexPad->SetTextSize(22);

    can->cd();
    gLatexCan->DrawLatexNDC(0.15,0.955,Form("p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

    const Int_t iNumPads = 1;
    TString sPadNames[iNumPads] = {"cent0"};
    TString sPadLabels[iNumPads] = {" 0-20% (V0A)"};

    SetStylePID();

    const Int_t iNumHistos = 6;
    TString* sHistoNames = gSpecies;
    TString* sHistoLabels = gSpeciesLabels;
    // TString sHistoNames[iNumHistos] = {"Charged","Pion","Kaon","K0s","Proton","Lambda","Phi"};
    // TString sHistoLabels[iNumHistos] = {"h^{#pm}","#pi^{#pm}","K^{#pm}","K^{0}_{S}","p/#bar{p}","#Lambda/#bar{#Lambda}","#phi"};

    TString sDrawOpt = "pe";

    // Color_t cols[] = {kGray+2, kRed+1, kGreen+3, kGreen+1, kBlue, kOrange+2};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,kFullDiamond};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,47};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,43,kOpenSquare};
    // Double_t size[] = { 1.0,1.2,1.3,1.0,1.4,1.0};

    // gStyle->SetPalette(kRainBow);
    // TColor::InvertPalette();

    TLegend* leg1 = new TLegend(0.16,0.12,0.4,0.33);
    // leg1->SetHeader("");
    SetLegendStyle(leg1);

    TLegend* leg2 = new TLegend(0.26,0.12,0.5,0.33);
    // leg2->SetHeader("");
    SetLegendStyle(leg2);

    // === Loading datapoints =================


    TFile* fileIn = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileName.Data()),"READ");
    if(!fileIn) { printf("Input file '%s' not openned!\n",sInputFileName.Data()); return; }

    TGraphErrors* gr = nullptr;

    // loop over pads
    for(Int_t p(0); p < iNumPads; ++p) {
        can->cd(p+1);

        TString sPadName = sPadNames[p];
        TString sPadLabel = sPadLabels[p];
        if(!sPadLabel.IsNull()) gLatexPad->DrawLatexNDC(0.94,0.9,sPadLabel.Data());

        // loop over histos within given pad
        for(Int_t h(0); h < iNumHistos; ++h) {

            TString sHistoName = sHistoNames[h];
            TString sHistoLabel = sHistoLabels[h];

            TString sGraphName = Form("%s_%s_%s", sHistoName.Data(), sCorrName.Data(), sPadName.Data());
            printf("%s\n",sGraphName.Data());

            gr = (TGraphErrors*) fileIn->Get(sGraphName.Data());
            if(!gr) {
                printf("Graph '%s' not found!\n",sGraphName.Data());
                // fileIn->ls();
                return;
            }

            if(h < 1) { gr->SetMarkerStyle(kOpenSquare); }

            SetGraphError(gr,0.0);

            SetMarkerPID(gr,kSpecies(h));
            gr->Draw(sDrawOpt.Data());

            if(p < 1 && !sHistoLabel.IsNull()) {
                if(h < iNumHistos / 2 ) {
                    leg1->AddEntry(gr,sHistoLabel.Data(),"p");
                } else {
                    leg2->AddEntry(gr,sHistoLabel.Data(),"p");
                }

            }
        }

    }

    can->cd(4);
    leg1->Draw();
    leg2->Draw();
    can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));
}
