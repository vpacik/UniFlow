enum PosLegend {kLegTopLeft = 1, kLegTopRight, kLegBotLeft, kLegBotRight};
TLegend* MakeLegend(PosLegend pos);

void SetStyle(Bool_t graypalette = 0);

TGraphErrors* Transfer(TH1D* histo, Double_t dSyst);
TGraphErrors* Transfer(TH1D* histo, TH1D* syst);

void PlotVnAllWithSystHistosPreliminary()
{

    SetStyle();

    Double_t dAlpha = 0.35;
    TString sInDir = "./";
    TString sInDirSyst = "./syst/";
    TString sOutDir = "./out/";

    Double_t xMin = 0.0; Double_t xMax = 10.0;
    Double_t yMin = 0.0; Double_t yMax = 0.35;
    PosLegend pos = kLegBotLeft;
    TString yAxisTitle = "";
    Bool_t bPhi = 1;

    // TString sLabel = "v_{2}{2}"; TString sHistName = "hFlow2_harm2_gap-10";
    // TString sLabel = "v_{2}{2,|#Delta#eta| > 0.0}"; TString sHistName = "hFlow2_harm2_gap00";
    // TString sLabel = "v_{2}{2,|#Delta#eta| > 0.4}"; TString sHistName = "hFlow2_harm2_gap04";
    // TString sLabel = "v_{2}{2,|#Delta#eta| > 0.8}"; TString sHistName = "hFlow2_harm2_gap08";  yMin = 0.0; yMax = 0.4; pos = kLegTopLeft;
    TString sLabel = "v_{2}{4}"; TString sHistName = "hFlow4_harm2_gap-10"; yMin = 0.0; yMax = 0.35;  pos = kLegTopRight; yAxisTitle = "v_{2}{4}";
    // TString sLabel = "v_{2}{4,|#Delta#eta| > 0}"; TString sHistName = "hFlow4_harm2_gap00"; yMin = 0.0; yMax = 0.4;  pos = kLegTopLeft;

    TString sHistNameSyst = "hFlow4_harm2_gap-10";


    TString sCent[] = {

        "5-10% V0M",
        "10-20% V0M",
        "20-30% V0M",
        "30-40% V0M",
        "40-50% V0M",
        "50-60% V0M",
    };

    std::vector<TString> sSpecies = {};
    std::vector<TString> sSpesLabels = {};
    std::vector<Color_t> colors = {};
    std::vector<Int_t> markers = {};
    std::vector<Double_t> markerSizes = {};
    std::vector<Double_t> systs = {};


    sSpecies.push_back("Charged"); sSpesLabels.push_back(" h^{#pm}"); colors.push_back(kGray+2); markers.push_back(kOpenCircle); markerSizes.push_back(0.8); systs.push_back(0.03);
    sSpecies.push_back("Pion"); sSpesLabels.push_back(" #pi^{#pm}"); colors.push_back(kRed+1); markers.push_back(kFullCircle); markerSizes.push_back(0.8); systs.push_back(0.03);
    sSpecies.push_back("Kaon"); sSpesLabels.push_back(" K^{#pm}"); colors.push_back(kGreen+3); markers.push_back(kFullTriangleUp); markerSizes.push_back(1); systs.push_back(0.04);
    sSpecies.push_back("K0s"); sSpesLabels.push_back(" K^{0}_{S}"); colors.push_back(kGreen+1); markers.push_back(kFullTriangleDown); markerSizes.push_back(1); systs.push_back(0.07);
    sSpecies.push_back("Proton"); sSpesLabels.push_back( " p(#bar{p})"); colors.push_back(kBlue); markers.push_back(kFullSquare); markerSizes.push_back(0.7); systs.push_back(0.04);
    sSpecies.push_back("Lambda"); sSpesLabels.push_back(" #Lambda(#bar{#Lambda})"); colors.push_back(kOrange+1); markers.push_back(kFullDiamond); markerSizes.push_back(1.4); systs.push_back(0.05);
    if(bPhi) { sSpecies.push_back("Phi"); sSpesLabels.push_back(" #phi"); colors.push_back(kMagenta); markers.push_back(kFullCross); markerSizes.push_back(1.); systs.push_back(0.08); }


    TFile* fileIn = TFile::Open(Form("%s/Processed.root",sInDir.Data()),"READ");
    if(!fileIn){  printf("ERROR: fileIn not found!\n");  return; }
    TFile* fileInPhi = TFile::Open(Form("%s/Processed.root",sInDir.Data()),"READ");
    if(!fileInPhi){  printf("ERROR: fileInPhi not found!\n");  return; }

    gSystem->mkdir(sOutDir,1);

    Int_t iNumMult = sizeof(sCent) / sizeof(sCent[0]);
    for(Int_t iMult(0); iMult < iNumMult; ++iMult) {

        TCanvas* can = new TCanvas("can","can",800,600);
        can->cd();
        TH1* frame = (TH1*) gPad->DrawFrame(xMin, yMin, xMax, yMax);
        frame->SetTitle(Form("; #it{p}_{T} (GeV/#it{c}); %s",yAxisTitle.Data()));

        // TLegend* leg = MakeLegend(pos);
        // TLegend* leg = new TLegend(0.14,0.12,0.5,0.4);
        TLegend* leg = new TLegend(0.67,0.6,0.84,0.82);
        leg->SetTextSize(0.035);
        leg->SetBorderSize(0);
        leg->SetFillColorAlpha(0.0,0.0);

        TLegend* leg2 = new TLegend(0.77,0.6,0.93,0.82);
        leg2->SetTextSize(0.035);
        leg2->SetBorderSize(0);
        leg2->SetFillColorAlpha(0.0,0.0);

        TList* listSyst = new TList();
        TList* listHistos = new TList();

        Int_t iNumPlots = sSpecies.size();
        for(Int_t i(0); i < iNumPlots; ++i) {

            Bool_t bPhi = (sSpecies[i].Contains("Phi"));

            TString sName = Form("%s_%s_cent%d",sSpecies[i].Data(), sHistName.Data(), iMult);
            TString sNameSyst = Form("%s_%s_cent%d",sSpecies[i].Data(), sHistNameSyst.Data(),1);

            TH1D* hist = (TH1D*) fileIn->Get(sName.Data());

            if(bPhi) {
                // if(iMult == 0 || iMult == iNumMult-1) continue;
                // sName = Form("%s_%s_cent%d",sSpecies[i].Data(), sHistName.Data(), iMult-1);
                hist = (TH1D*) fileInPhi->Get(sName.Data());
            }

            if(!hist) { printf("ERROR: histo '%s' not found!\n",sName.Data()); fileIn->ls(); return; }

            TFile* fileSyst = TFile::Open(Form("%s/syst_total_%s.root",sInDirSyst.Data(),sSpecies[i].Data()),"READ");
            if(!fileSyst) {  printf("ERROR: syst file not found!\n");  return; }

            TH1D* histSyst = (TH1D*) fileSyst->Get(sNameSyst.Data());
            if(!histSyst) { printf("ERROR: syst histo '%s' not found!\n",sNameSyst.Data()); fileSyst->ls(); return; }

            // TGraphErrors* graphSyst = Transfer(hist,systs.at(i));
            TGraphErrors* graphSyst = Transfer(hist,histSyst);
            if(!graphSyst) { printf("ERROR: syst not found!\n"); return; }

            graphSyst->SetLineColor(colors[i]);
            graphSyst->SetFillColorAlpha(colors[i],dAlpha);
            // graphSyst->SetMarkerColor(colors[i]);
            // graphSyst->SetMarkerStyle(markers[i]);
            // graphSyst->SetMarkerSize(markerSizes[i]);
            // graphSyst->Draw("same 2");
            listSyst->Add(graphSyst);

            hist->SetLineColor(colors[i]);
            hist->SetFillColorAlpha(colors[i],dAlpha);
            hist->SetMarkerColor(colors[i]);
            hist->SetMarkerStyle(markers[i]);
            hist->SetMarkerSize(markerSizes[i]);
            // hist->Draw("hist ape1x0 same");
            listHistos->Add(hist);

            if(i < 4) {
                leg->AddEntry(hist,sSpesLabels[i].Data(),"pf");
            } else {
                leg2->AddEntry(hist,sSpesLabels[i].Data(),"pf");
            }
        }

        TObject* o = nullptr;
        leg2->AddEntry(o," ","");
        if(!bPhi) { leg2->AddEntry(o," ",""); } // nophi


        // leg->SetHeader(Form("%s (%s)",sLabel.Data(),sCent[iMult].Data()));
        leg->Draw();
        leg2->Draw();

        TLatex * text = new TLatex();
        TLatex * text2 = new TLatex();
        text->SetTextFont(42);
        text->SetTextSize(0.04);
        text2->SetTextFont(42);
        text2->SetTextSize(0.04);

        // frame_canPlot->SetTitle(Form("; #it{p}_{T} (GeV/#it{c}); v_{2}^{sub} {2, |#Delta#eta| > %s}",sGapVal.Data()));
        text2->DrawLatexNDC(0.186,0.84,Form("ALICE Preliminary"));
        text->DrawLatexNDC(0.186,0.79,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV"));
        text->DrawLatexNDC(0.186,0.74,Form("|#it{#eta}| < 0.8"));
        text->DrawLatexNDC(0.72,0.84,sCent[iMult].Data());

        // Drawing systematics first, then histos

        for(Int_t i(0); i < iNumPlots; ++i) {
            TGraphErrors* gr = (TGraphErrors*) listSyst->At(i);
            gr->Draw("same 2");
        }
        for(Int_t i(0); i < iNumPlots; ++i) {
            TH1D* h = (TH1D*) listHistos->At(i);
            h->Draw("hist ape1x0 same");
        }


        can->SaveAs(Form("%s/can_%s_all_cent%d.pdf",sOutDir.Data(),sLabel.Data(),iMult),"pdf");

        delete leg;
        delete frame;
        delete can;
    }

    return;
}

TLegend* MakeLegend(PosLegend pos)
{
    switch(pos) {
        case kLegTopLeft: return new TLegend(0.14,0.65,0.5,0.88); break;
        case kLegTopRight: return new TLegend(0.68,0.65,0.88,0.88); break;
        case kLegBotLeft: return new TLegend(0.14,0.12,0.5,0.4); break;
        case kLegBotRight: return new TLegend(0.68,0.12,0.88,0.4); break;
        default: return nullptr;
    }
}

TGraphErrors* Transfer(TH1D* histo, Double_t dSyst)
{
    if(!histo) { printf("E::Transfer : Input not found!\n"); return nullptr; }

    Double_t dSystRel = dSyst;

    Double_t dX = 0.0;
    Double_t dY = 0.0;

    TGraphErrors* err = new TGraphErrors(histo);
    // SetEx(err,0.05);
    Int_t np = err->GetN();
    for (Int_t i=0; i<np; i++) {
        err->GetPoint(i,dX,dY);
        err->SetPointError(i,0.09, dSystRel*dY);
    }

    return err;
}

TGraphErrors* Transfer(TH1D* histo, TH1D* syst)
{
    if(!histo) { printf("E::Transfer : Input not found!\n"); return nullptr; }
    if(!syst) { printf("E::Transfer : Input syst histo not found!\n"); return nullptr; }

    Double_t dSystRel = 0.0;

    Double_t dX = 0.0;
    Double_t dY = 0.0;

    TGraphErrors* err = new TGraphErrors(histo);
    // SetEx(err,0.05);
    Int_t np = err->GetN();
    for (Int_t i=0; i<np; i++) {
        err->GetPoint(i,dX,dY);
        dSystRel = syst->GetBinContent(i+1);
        err->SetPointError(i,0.09, dSystRel*dY);
    }

    return err;
}


void SetStyle(Bool_t graypalette) {
  // cout << "Setting style!" << std::end;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}
//_____________________________________________________________________________
