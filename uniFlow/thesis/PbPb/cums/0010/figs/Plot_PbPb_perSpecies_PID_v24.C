// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStyle.C"


void Plot_PbPb_perSpecies_PID_v24()
{
    // === Preparing canvas =================

    TString sInputFileName = "PbPb/cums/0010/v2/graphs_PbPb_v22_v24.root";
    TString sOutPlotName = "PbPb_perSpecies_PID_v24.pdf";
    // TString sOutPlotForm = "PbPb_v24_perCent";

    Double_t dXmin = 0.0;
    Double_t dXmax = 10.0;
    Double_t dYmin = -0.03;
    Double_t dYmax = 0.4;

    TGaxis* gaX = new TGaxis(1.0, dYmin, 9.0, dYmin, 1.0, 9.0);
    gaX->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TGaxis* gaY = new TGaxis(dXmin, 0.0, dXmin, 0.35, 0.0,0.35);
    gaY->SetTitle("#it{v}_{2}{4}");
    TString sCorrName = "v2{4}";

    TCanvas* can = new TCanvas("can","can",900,600);
    TList* frames = PrepareCanvas3x2(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    can->cd();
    gLatexCan->DrawLatexNDC(0.09,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

    const Int_t iNumPads = 6;
    TString sPadNames[iNumPads] = {"Pion","Kaon","K0s","Proton","Lambda","Phi"};
    TString sPadLabels[iNumPads] = {"#pi^{#pm}","K^{#pm}","K^{0}_{S}","p/#bar{p}","#Lambda/#bar{#Lambda}","#phi"};

    const Int_t iNumHistos = 6;
    TString sHistoNames[iNumHistos] = {"cent0010","cent1020","cent2030","cent3040","cent4050","cent5060"};
    TString sHistoLabels[iNumHistos] = {" 0-10%","10-20%","20-30%","30-40%","40-50%","50-60%"};

    TString sDrawOpt = "pe PLC PMC";

    Color_t cols[] = {kGray+2, kRed+1, kGreen+3, kGreen+1, kBlue, kOrange+2};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,kFullDiamond};
    // Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,kOpenSquare,47};
    Int_t mark[] = { kFullSquare,kFullCross,kFullDiamond,kFullCircle,43,kOpenSquare};
    Double_t size[] = { 1.0,1.2,1.3,1.0,1.4,1.0};

    gStyle->SetPalette(kRainBow);
    TColor::InvertPalette();

    TLegend* leg1 = new TLegend(0.27,0.65,0.52,0.96);
    leg1->SetHeader("Event class (V0M)");
    SetLegendStyle(leg1);

    TLegend* leg2 = new TLegend(0.55,0.65,0.82,0.96);
    leg2->SetHeader("");
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
        if(!sPadLabel.IsNull()) gLatexPad->DrawLatexNDC(0.94,0.94,sPadLabel.Data());


        // loop over histos within given pad
        for(Int_t h(0); h < iNumHistos; ++h) {

            TString sHistoName = sHistoNames[h];
            TString sHistoLabel = sHistoLabels[h];

            TString sGraphName = Form("%s_%s_%s", sPadName.Data(), sCorrName.Data(), sHistoName.Data());
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
                if(h < iNumPads / 2) {
                    leg1->AddEntry(gr,sHistoLabel.Data(),"pl");
                } else {
                    leg2->AddEntry(gr,sHistoLabel.Data(),"pl");
                }

            }
        }
    }

    can->cd(1);
    leg1->Draw();
    leg2->Draw();
    can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));
}
