// centrality dependence (pads per species)

#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libCanvas.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStylePID.C"

TGraphErrors* TransferSyst(TGraphErrors* histo, TH1D* syst);

void Plot_PbPb_perCent_All_v24_KEt_syst()
{
    // === Preparing canvas =================

    TString sInputFileName = "PbPb/cums/0010/v2/graphs_PbPb_v22_v24_scaled.root";
    TString sOutPlotName = "PbPb_perCent_All_v24_KEt_syst.pdf";
    // TString sOutPlotForm = "PbPb_v24_perCent";
    TString sInputFileNameSyst = "PbPb/cums/0510/preliminary/syst/";


    Double_t dXmin = 0.0;
    Double_t dXmax = 5.0;
    Double_t dYmin = -0.01;
    Double_t dYmax = 0.13;

    TGaxis* gaX = new TGaxis(0.5, dYmin, 4.5, dYmin, 0.5, 4.5);
    gaX->SetTitle("(#it{m}_{T} - #it{m}_{0})/#it{n}_{q} (GeV/#it{c}^{2})");
    TGaxis* gaY = new TGaxis(dXmin, 0.0, dXmin, 0.12, 0.0,0.12);
    gaY->SetTitle("#it{v}_{2}{4}/#it{n}_{q}");
    TString sCorrName = "v2{4}_KEt";

    TCanvas* can = new TCanvas("can","can",900,600);
    TList* frames = PrepareCanvas3x2(can, gaX, gaY, dXmin, dYmin, dXmax, dYmax);

    can->cd();
    gLatexCan->DrawLatexNDC(0.09,0.955,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, |#eta| < 0.8 (This thesis)"));
    // gLatexCan->DrawLatexNDC(0.8,0.955,Form("This thesis"));

    const Int_t iNumPads = 6;
    TString sPadNames[iNumPads] = {"cent0010","cent1020","cent2030","cent3040","cent4050","cent5060"};
    TString sPadLabels[iNumPads] = {" 0-10%","10-20%","20-30%","30-40%","40-50%","50-60%"};

    SetStylePID();

    const Int_t iNumHistos = 7;
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

    // === Loading datapoints =================

    TFile* fileIn = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s",sInputFileName.Data()),"READ");
    if(!fileIn) { printf("Input file '%s' not openned!\n",sInputFileName.Data()); return; }


    TGraphErrors* gr = nullptr;
    TGraphErrors* grSyst = nullptr;
    Double_t dAlpha = 0.35;


    // loop over pads
    for(Int_t p(0); p < iNumPads; ++p) {
        can->cd(p+1);

        TList* list = new TList();
        TList* listSyst = new TList();

        TString sPadName = sPadNames[p];
        TString sPadLabel = sPadLabels[p];
        if(!sPadLabel.IsNull()) gLatexPad->DrawLatexNDC(0.94,0.94,sPadLabel.Data());


        // loop over histos within given pad
        for(Int_t h(0); h < iNumHistos; ++h) {

            if(h < 1) continue; //skipping Charged

            TString sHistoName = sHistoNames[h];
            TString sHistoLabel = sHistoLabels[h];

            TString sGraphName = Form("%s_%s_%s", sHistoName.Data(), sCorrName.Data(), sPadName.Data());
            // printf("%s\n",sGraphName.Data());

            gr = (TGraphErrors*) fileIn->Get(sGraphName.Data());
            if(!gr) {
                printf("Graph '%s' not found!\n",sGraphName.Data());
                // fileIn->ls();
                return;
            }

            SetGraphError(gr,0.0);

            SetMarkerPID(gr,kSpecies(h));
            gr->SetFillColorAlpha(gCol[h],dAlpha);
            // gr->Draw(sDrawOpt.Data());
            list->Add(gr);





            if(p < 1 && !sHistoLabel.IsNull()) {
                if(h < iNumPads / 2 + 1) {
                    leg1->AddEntry(gr,sHistoLabel.Data(),"pf");
                } else {
                    leg2->AddEntry(gr,sHistoLabel.Data(),"pf");
                }

            }

            // loading systematics

            TFile* fileInSyst = TFile::Open(Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/%s/syst_total_%s.root",sInputFileNameSyst.Data(),sHistoName.Data()),"READ");
            if(!fileInSyst) { printf("Input file '%s' (%s) not openned!\n",sInputFileNameSyst.Data(),sHistoName.Data()); return; }
            // fileInSyst->ls();
            TString sSystName = Form("%s_hFlow4_harm2_gap-10_cent%d",sHistoName.Data(),p);
            TH1D* histSyst = (TH1D*) fileInSyst->Get(sSystName.Data());
            if(!histSyst) { printf("Input histo %s syst not openned\n",sSystName.Data()); }
            grSyst = TransferSyst(gr,histSyst);
            if(!grSyst)  { printf("Syst graph not transfered histo %s \n",sSystName.Data()); }

            if(grSyst) {
                grSyst->SetFillColorAlpha(gCol[h],dAlpha);
                // grSyst->Draw("same 2");
                listSyst->Add(grSyst);
            }

        }



        for(Int_t h(0); h < iNumHistos-1; ++h) {
            grSyst = (TGraphErrors*) listSyst->At(h);
            grSyst->Draw("same 2");
        }

        for(Int_t h(0); h < iNumHistos-1; ++h) {
            gr = (TGraphErrors*) list->At(h);
            gr->Draw(sDrawOpt.Data());
        }
    }

    can->cd(1);
    leg1->Draw();
    leg2->Draw();
    can->SaveAs(Form("./plots/%s",sOutPlotName.Data()));
}

TGraphErrors* TransferSyst(TGraphErrors* histo, TH1D* syst)
{
    if(!histo) { printf("E::Transfer : Input not found!\n"); return nullptr; }
    if(!syst) { printf("E::Transfer : Input syst histo not found!\n"); return nullptr; }

    Double_t dSystRel = 0.0;

    Double_t dX = 0.0;
    Double_t dY = 0.0;
    Double_t dET = 0.0;

    TGraphErrors* err = (TGraphErrors*) histo->Clone();
    // SetEx(err,0.05);
    Int_t np = err->GetN();
    Int_t bins = syst->GetXaxis()->GetNbins();
    // printf("bins: %d | points %d\n",np,bins);

    for (Int_t i=0; i<np; i++) {
        err->GetPoint(i,dX,dY);
        dET = err->GetErrorY(i);
        dSystRel = syst->GetBinContent(i+1);
        err->SetPointError(i,0.08, dSystRel*dY);

        // printf("%f : value %f +- %f += %f (%f)\n",dX,dY,dET,dY*dSystRel,dSystRel);

    }

    return err;
}
