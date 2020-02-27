// Macro for plotting two panel figure with inv. mass distribution & corr-mass figure
// Author: Vojtech Pacik (2019)

void SetStyle(Bool_t graypalette=kFALSE); // official ALICE fig. template

void PlotFlowMass() {

    SetStyle();

    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/non_linear_modes/output/Phi/";
    TString sListName = "Phi_<<3>>(4,-2,-2)_2sub(0)_cent2_pt3;1";
    Bool_t bPlotSig = 1;
    Bool_t bPrel = kFALSE;
    Bool_t bDrawLegs = 1;


    TFile* fileIn = TFile::Open(Form("%s/fits.root",sPath.Data()),"READ");
    if(!fileIn) { printf("FileIn not open!\n"); return; }

    TList* listIn = (TList*) fileIn->Get(sListName.Data());
    if(!listIn) { printf("List '%s' not found!\n",sListName.Data()); fileIn->ls(); return; }

    TH1D* hMass = (TH1D*) listIn->FindObject("histMass");
    if(!hMass) { printf("hMass not found!\n"); listIn->ls(); return; }

    TH1D* hCorr = (TH1D*) listIn->FindObject("histCorr");
    if(!hCorr) { printf("hCorr not found!\n"); listIn->ls(); return; }

    hMass->SetMarkerStyle(kFullCircle);
    hMass->SetMarkerSize(0.65);
    hMass->SetMarkerColor(kBlue+1);
    hMass->SetLineColor(kBlue+1);
    hMass->SetLineWidth(2);

    hCorr->SetLineColor(kGreen+2);
    hCorr->SetMarkerColor(kGreen+2);
    hCorr->SetMarkerStyle(21);

    TH1D* hMassSig = nullptr;
    TF1* fitMassSig = nullptr;

    if(bPlotSig) {
        hMassSig = (TH1D*) listIn->FindObject("histMass_sig");
        if(!hMassSig) { printf("hMassSig not found!\n"); listIn->ls(); return; }
        hMassSig->SetMarkerStyle(kOpenCircle);
        hMassSig->SetMarkerSize(0.65);
        hMassSig->SetMarkerColor(kRed);
        hMassSig->SetLineColor(kRed);
        hMassSig->SetLineWidth(2);

        fitMassSig = (TF1*) listIn->FindObject("fitMass_sig");
        if(!fitMassSig) { printf("fitMassSig not found!\n"); listIn->ls(); return; }
        fitMassSig->SetLineStyle(2);
        fitMassSig->SetLineWidth(2);
        fitMassSig->SetLineColor(kBlack);

    }


    TF1* fitMass = (TF1*) listIn->FindObject("fitMass");
    if(!fitMass) { printf("fitMass not found!\n"); listIn->ls(); return; }

    TF1* fitCorr = (TF1*) listIn->FindObject("fitCorr");
    if(!fitCorr) { printf("fitCorr not found!\n"); listIn->ls(); return; }


    // Setting the canvas
    Double_t dXmin = fitMass->GetXmin();
    Double_t dXmax = fitMass->GetXmax();
    // Double_t dXmin = hMass->GetXaxis()->GetXmin();
    // Double_t dXmax = hMass->GetXaxis()->GetXmax();

    Double_t dYmin = -5000;
    Double_t dYmax = 220000;

    printf("X: min : %f | max %f\n",dXmin,dXmax);
    printf("Y: min : %f | max %f\n",dYmin,dYmax);


    Double_t dXmin2 = dXmin;
    Double_t dXmax2 = dXmax;
    Double_t dYmin2 = 20e-5;
    Double_t dYmax2 = 32e-5;
    // Double_t dYmin2 = hCorr->GetMinimum();
    // Double_t dYmax2 = hCorr->GetMaximum();

    TLatex* xALICE = new TLatex(0.73,0.8, bPrel ? "ALICE Preliminary" : "ALICE");
    xALICE->SetNDC();
    xALICE->SetTextFont(42);
    xALICE->SetTextSize(0.08);

    TLatex* xSyst = new TLatex(0.19,0.8, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    xSyst->SetNDC();
    xSyst->SetTextFont(42);
    xSyst->SetTextSize(0.06);

    TLatex* xKine = new TLatex(0.19,0.66, "3 < #it{p}_{T} < 4.5 GeV/#it{c}");
    xKine->SetNDC();
    xKine->SetTextFont(42);
    xKine->SetTextSize(0.06);

    TLatex* xKine2 = new TLatex(0.19,0.73, "10-20%, |#eta| < 0.8");
    xKine2->SetNDC();
    xKine2->SetTextFont(42);
    xKine2->SetTextSize(0.06);

    TLegend* leg1 = new TLegend(0.55,0.1,0.85,0.5);
    leg1->SetFillColorAlpha(0,0);
    // leg1->SetHeader("0.2 < #it{p}_{T} < 3.0 GeV/#it{c}^{2}, |#eta| < 0.8");
    leg1->AddEntry(hMass,"K^{+}+K^{-} yield","p");
    // leg1->AddEntry(fitMass,"fit","l");
    if(bPlotSig) {
        leg1->AddEntry(hMassSig,"#phi #rightarrow K^{+} + K^{-}","p");
        leg1->AddEntry(fitMassSig,"Breit-Wigner fit","l");
    }
    leg1->AddEntry(fitCorr,"#it{d}_{4,22}^{total} fit","l");

    TLegend* leg2 = new TLegend(0.22,0.2,0.4,0.5);
    leg2->SetFillColorAlpha(0,0);
    // leg2->SetHeader("0.2 < #it{p}_{T} < 3.0 GeV/#it{c}^{2}, |#eta| < 0.8");
    // leg2->AddEntry(histMass,"K0s cand","pel");
    // leg2->AddEntry(fitCorr,"fit","l");

    // Define the Canvas
    TCanvas* c = new TCanvas("c", "canvas", 600, 600);


    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.4, 1, 1.0);
    pad1->SetBottomMargin(0.00); // 0.01
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.4);
    pad2->SetTopMargin(0.00); //0.01
    pad2->SetBottomMargin(0.3);
    pad2->Draw();

    pad1->cd();
    TH1* frame1 = (TH1*) gPad->DrawFrame(dXmin, dYmin, dXmax, dYmax);
    hMass->Draw("same ep");
    // fitMass->Draw("same");
    if(bPlotSig) {
        hMassSig->Draw("same ep");
        fitMassSig->Draw("same");
    }
    xALICE->Draw();
    xSyst->Draw();
    xKine->Draw();
    xKine2->Draw();
    if(bDrawLegs) { leg1->Draw(); }


    pad2->cd();
    TH1* frame2 = (TH1*) gPad->DrawFrame(dXmin2, dYmin2, dXmax2, dYmax2);
    hCorr->Draw("same ep");
    fitCorr->Draw("same");
    // if(bDrawLegs) { leg2->Draw(); }

    // TGaxis *axis = new TGaxis( dXmin, dYmin, dXmin, dYmax, dYmin,dYmax,510,"");
    // axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    // axis->SetLabelSize(15);
    // axis->Draw();

    // TGaxis::SetExponentOffset(0.0,0.0);

    // Main frame (frame1) setting
    frame1->GetXaxis()->SetLabelSize(0.0); // // hiding X axis label (due to margin)

    frame1->GetYaxis()->SetTitle("Counts");
    frame1->GetYaxis()->SetTitleSize(20);
    frame1->GetYaxis()->SetTitleFont(43);
    frame1->GetYaxis()->SetTitleOffset(1.55);
    frame1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame1->GetYaxis()->SetLabelSize(15);
    //
    // // Ratio plot (frame2) settings
    frame2->SetTitle(""); // Remove the ratio title

    frame2->GetYaxis()->SetTitle("#it{d}^{total}_{4,22} (#times10^{-3})");
    // frame2->GetYaxis()->SetNoExponent();
    frame2->GetYaxis()->SetNdivisions(505);
    frame2->GetYaxis()->SetTitleSize(20);
    frame2->GetYaxis()->SetTitleFont(43);
    frame2->GetYaxis()->SetTitleOffset(1.55);
    frame2->GetYaxis()->SetLabelFont(43);
    frame2->GetYaxis()->SetLabelSize(15);

    frame2->GetXaxis()->SetTitle("#it{M}_{inv} (GeV/#it{c}^{2})");
    frame2->GetXaxis()->SetTitleSize(20);
    frame2->GetXaxis()->SetTitleFont(43);
    frame2->GetXaxis()->SetTitleOffset(2.4);
    frame2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    frame2->GetXaxis()->SetLabelSize(15);

    c->SaveAs(Form("%s/flowmass_new_%s.pdf",sPath.Data(),sListName.Data()),"pdf");

}
// =============================================================================
void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;

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
  // gStyle->SetLabelFont(43,"xyz");
  // gStyle->SetLabelSize(15,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  // gStyle->SetTitleSize(20,"xyz");
  // gStyle->SetTitleFont(43,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  // gStyle->SetTextFont(43);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}
