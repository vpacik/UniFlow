// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStylePID.C"


void Plot_PbPb_perCent_All_v24_theory()
{
    // === Preparing canvas =================

    TString sInputFileNameModel = "PbPb/models/graphs_PbPb_v22_v24_AMPT.root";
    TString sInputFileNameTRENTo = "PbPb/models/graphs_PbPb_v22_v24_TRENTo.root";

    TString sInputFileName = "PbPb/cums/0010/v2/graphs_PbPb_v22_v24.root";
    TString sOutPlotName = "PbPb_perCent_All_v24_theory.pdf";
    // TString sOutPlotForm = "PbPb_v24_perCent";

    Double_t dXmin = 0.0;
    Double_t dXmax = 3.0;
    Double_t dYmin = -0.03;
    Double_t dYmax = 0.33;

    TGaxis* gaX = new TGaxis(0.2, dYmin, 2.8, dYmin, 0.2, 2.8);
    gaX->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TGaxis* gaY = new TGaxis(dXmin, 0.0, dXmin, 0.30, 0.0,0.30);
    gaY->SetTitle("#it{v}_{2}{4}");
    TString sCorrName = "v2{4}";

    TCanvas* can = new TCanvas("can","can",900,600);
    TList* frames = PrepareCanvas3x2(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    can->cd();
    gLatexCan->DrawLatexNDC(0.09,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

    const Int_t iNumPads = 6;
    TString sPadNames[iNumPads] = {"cent0010","cent1020","cent2030","cent3040","cent4050","cent5060"};
    TString sPadLabels[iNumPads] = {" 0-10%","10-20%","20-30%","30-40%","40-50%","50-60%"};

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

    TLegend* leg1 = new TLegend(0.27,0.65,0.5,0.96);
    // leg1->SetHeader("");
    SetLegendStyle(leg1);

    TLegend* leg2 = new TLegend(0.45,0.65,0.68,0.96);
    // leg2->SetHeader("");
    SetLegendStyle(leg2);

    TLegend* legAMPT = new TLegend(0.27,0.69,0.57,0.96);
    legAMPT->SetHeader("AMPT");
    SetLegendStyle(legAMPT);

    TLegend* legTRENTo = new TLegend(0.45,0.69,0.75,0.96);
    legTRENTo->SetHeader("T_{R}ENTo");
    SetLegendStyle(legTRENTo);

    // === Loading datapoints =================

    TFile* fileInModel = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileNameModel.Data()),"READ");
    if(!fileInModel) { printf("Input file '%s' not openned!\n",sInputFileNameModel.Data()); return; }

    TFile* fileInTRENTo = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileNameTRENTo.Data()),"READ");
    if(!fileInTRENTo) { printf("Input file '%s' not openned!\n",sInputFileNameTRENTo.Data()); return; }

    TFile* fileIn = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileName.Data()),"READ");
    if(!fileIn) { printf("Input file '%s' not openned!\n",sInputFileName.Data()); return; }

    TGraphErrors* gr = nullptr;

    // loop over pads
    for(Int_t p(0); p < iNumPads; ++p) {
        can->cd(p+1);

        TString sPadName = sPadNames[p];
        TString sPadLabel = sPadLabels[p];
        if(!sPadLabel.IsNull()) gLatexPad->DrawLatexNDC(0.94,0.94,sPadLabel.Data());

        if(p > 0) {

            // plotting models
            TGraphErrors* grModPion = (TGraphErrors*) fileInModel->Get(Form("pion_v24_%s",sPadName.Data()));
            if(!grModPion) { printf("Pion model not found!\n"); fileInModel->ls(); return; }

            TGraphErrors* grModKaon = (TGraphErrors*) fileInModel->Get(Form("kaon_v24_%s",sPadName.Data()));
            if(!grModKaon) { printf("Kaon model not found!\n"); fileInModel->ls(); return; }

            TGraphErrors* grModProton = (TGraphErrors*) fileInModel->Get(Form("proton_v24_%s",sPadName.Data()));
            if(!grModProton) { printf("Proton model not found!\n"); fileInModel->ls(); return; }

            Double_t dAlpha = 1.0;

            grModPion->SetFillColorAlpha(kRed,dAlpha);
            grModPion->SetLineColor(kRed);
            // grModPion->SetFillStyle(1001);
            grModPion->SetFillStyle(3345);
            grModPion->Draw("3");

            grModKaon->SetFillColorAlpha(kGreen-2,dAlpha);
            grModKaon->SetLineColor(kGreen-2);
            // grModKaon->SetFillStyle(1001);
            grModKaon->SetFillStyle(3357);
            grModKaon->Draw("3");

            grModProton->SetFillColorAlpha(kBlue-4,0.3);
            grModProton->SetLineColor(kBlue-4);
            // grModProton->SetFillStyle(322);
            grModProton->SetFillStyle(1001);
            // grModProton->SetFillStyle(3020);
            grModProton->Draw("3");

            // plotting TRENTO
            TGraphErrors* grModPionTRENTo = (TGraphErrors*) fileInTRENTo->Get(Form("pion_v24_%s",sPadName.Data()));
            if(!grModPionTRENTo) { printf("Pion model not found!\n"); fileInTRENTo->ls(); return; }

            TGraphErrors* grModKaonTRENTo = (TGraphErrors*) fileInTRENTo->Get(Form("kaon_v24_%s",sPadName.Data()));
            if(!grModKaonTRENTo) { printf("Kaon model not found!\n"); fileInTRENTo->ls(); return; }

            TGraphErrors* grModProtonTRENTo = (TGraphErrors*) fileInTRENTo->Get(Form("proton_v24_%s",sPadName.Data()));
            if(!grModProtonTRENTo) { printf("Proton model not found!\n"); fileInTRENTo->ls(); return; }


            grModPionTRENTo->SetFillColorAlpha(kRed,1.0);
            grModPionTRENTo->SetLineColor(kRed);
            grModPionTRENTo->SetFillStyle(322);
            grModPionTRENTo->Draw("3");

            grModKaonTRENTo->SetFillColorAlpha(kGreen-2,1.0);
            grModKaonTRENTo->SetLineColor(kGreen-2);
            grModKaonTRENTo->SetFillStyle(322);
            grModKaonTRENTo->Draw("3");

            grModProtonTRENTo->SetFillColorAlpha(kBlue,1.0);
            // grModProton->SetFillStyle(322);
            grModProtonTRENTo->SetLineColor(kBlue);
            grModProtonTRENTo->SetFillStyle(322);
            grModProtonTRENTo->Draw("3");

            if(p == 1) {
                // legAMPT->AddEntry(grModPion," ","f");
                // legAMPT->AddEntry(grModKaon," ","f");
                // legAMPT->AddEntry(grModProton," ","f");
                legAMPT->AddEntry(grModPion,gSpeciesLabels[kPion].Data(),"f");
                legAMPT->AddEntry(grModKaon,gSpeciesLabels[kKaon].Data(),"f");
                legAMPT->AddEntry(grModProton,gSpeciesLabels[kProton].Data(),"f");

                legTRENTo->AddEntry(grModPionTRENTo,gSpeciesLabels[kPion].Data(),"f");
                legTRENTo->AddEntry(grModKaonTRENTo,gSpeciesLabels[kKaon].Data(),"f");
                legTRENTo->AddEntry(grModProtonTRENTo,gSpeciesLabels[kProton].Data(),"f");

                // legTRENTo->AddEntry(grModPionTRENTo," ","f");
                // legTRENTo->AddEntry(grModKaonTRENTo," ","f");
                // legTRENTo->AddEntry(grModProtonTRENTo," ","f");
            }

        }

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
                if(h < iNumPads / 2) {
                    leg1->AddEntry(gr,sHistoLabel.Data(),"p");
                } else {
                    leg2->AddEntry(gr,sHistoLabel.Data(),"p");
                }

            }
        }

    }

    can->cd(1);
    leg1->Draw();
    leg2->Draw();
    can->cd(4);
    legAMPT->Draw();
    legTRENTo->Draw();
    can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));
}
