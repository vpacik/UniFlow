// Macro for plotting fits within fits.root file from ProcessUniFlow
// 2019 for thesis

#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"

void SetPad();
void SetHistOff(TH1* hist);
void SetFuncAtt(TF1* func, Color_t color, Style_t style, Int_t width);
void SetHistAtt(TH1* hist, Color_t color, Style_t markStyle, Double_t markSize);

void PlotFitsSingleMass(
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/6815/thesis/",
    TString sCorrName = "<<2>>(2,-2)",
    TString sSpecies = "K0s",
    Int_t iPt = 8,
    Int_t iCent = 5,
    Bool_t bPhi = 0)
{

    TString sFileName = "fits.root";
    TString gOutFormat = "pdf";

    TString sName = Form("%s_%s_cent%d_pt%d",sSpecies.Data(),sCorrName.Data(),iCent,iPt);

    TFile* fInput = TFile::Open(Form("%s/%s",sPath.Data(),sFileName.Data()),"READ");
    if(!fInput) { printf("ERROR: File not open!\n"); return; }

    TList* list = (TList*) fInput->Get(sName.Data());
    if(!list) { printf("ERROR: Input list '%s' not found!\n",sName.Data()); fInput->ls(); return; }

    // list->ls();

    TString sCorr = ((TNamed*) list->FindObject("corr"))->GetTitle();
    TString sSpec = ((TNamed*) list->FindObject("spec"))->GetTitle();
    TString sCent = ((TNamed*) list->FindObject("cent"))->GetTitle();
    TString sPt = ((TNamed*) list->FindObject("pt"))->GetTitle();
    TString sPtLow = ((TNamed*) list->FindObject("ptLow"))->GetTitle();
    TString sPtHigh = ((TNamed*) list->FindObject("ptHigh"))->GetTitle();

    // ========== Processing correlation =======================================

    TH1D* histMass = (TH1D*) list->FindObject("histMass");
    if(!histMass) { printf("ERROR: Object '%s' not found!\n", "histMass"); list->ls(); return kFALSE; }

    TF1* fitMass = (TF1*) list->FindObject("fitMass");
    // TF1* fitMass = (TF1*) list->FindObject("fitMass_tot");
    if(!fitMass) { printf("ERROR: Object '%s' not found!\n", "fitMass"); list->ls(); return kFALSE; }

    TF1* fitMassSig = (TF1*) list->FindObject("fitMass_partSig");
    // TF1* fitMassSig = (TF1*) list->FindObject("fitMass_sig");
    if(!fitMassSig) { printf("ERROR: Object '%s' not found!\n", "fitMassSig"); list->ls(); return kFALSE; }

    TF1* fitMassBg = (TF1*) list->FindObject("fitMass_partBg");
    // TF1* fitMassBg = (TF1*) list->FindObject("fitMass_bg");
    if(!fitMassBg) { printf("ERROR: Object '%s' not found!\n", "fitMassBg"); list->ls(); return kFALSE; }

    TH1D* histFracSig = (TH1D*) list->FindObject("histMass_fracSig");
    if(!histFracSig) { printf("ERROR: Object '%s' not found!\n", "histFracSig"); list->ls(); return kFALSE; }

    TH1D* histFracBg = (TH1D*) list->FindObject("histMass_fracBg");
    if(!histFracBg) { printf("ERROR: Object '%s' not found!\n", "histFracBg"); list->ls(); return kFALSE; }

    TF1* fitFracSig = (TF1*) list->FindObject("fitMass_fracSig");
    if(!fitFracSig) { printf("ERROR: Object '%s' not found!\n", "fitFracSig"); list->ls(); return kFALSE; }

    TF1* fitFracBg = (TF1*) list->FindObject("fitMass_fracBg");
    if(!fitFracBg) { printf("ERROR: Object '%s' not found!\n", "fitFracBg"); list->ls(); return kFALSE; }

    Double_t dXlow = fitMass->GetXmin();
    Double_t dXhigh = fitMass->GetXmax();

    // histMass->SetTitle(sTitle.Data());
    histMass->SetTitle("");
    histMass->GetYaxis()->SetTitle("Counts");

    histMass->GetXaxis()->SetRangeUser(dXlow,dXhigh);
    histMass->SetStats(0);

    Double_t dCorMax = histMass->GetMaximum();

    // histMass->SetMinimum(-0.05*dCorMax);
    // histMass->SetMinimum(0.0);
    histMass->SetMinimum(0.0);
    histMass->SetMaximum(1.2*dCorMax);

    SetHistAtt(histMass, kBlack, kFullCircle, 0.7);
    SetFuncAtt(fitMass, kRed, kSolid, 3);
    SetFuncAtt(fitMassSig, kGreen-2, 7, 3);
    SetFuncAtt(fitMassBg, kBlue, 2, 3);

    SetHistAtt(histFracSig, kBlack, kFullCircle, 0.7);
    SetHistAtt(histFracBg, kBlack, kFullCircle, 0.7);
    SetFuncAtt(fitFracSig, kGreen-2, 7, 3);
    SetFuncAtt(fitFracBg, kBlue, 2, 3);


    TLatex events;
    // events.SetNDC(1);
    events.SetTextFont(42);
    events.SetTextSize(0.04);
    // events.SeteventsAlign(30);

    TLatex species;
    species.SetTextFont(42);
    species.SetTextSize(0.06);

    gSystem->mkdir(Form("%s/fits/mass/single_nice/",sPath.Data()),1);

    TCanvas* canMass = new TCanvas("canMass","canMass",600,600);
    canMass->cd();
    SetPad();
    // gPad->SetLogy();
    histMass->DrawCopy();
    fitMass->Draw("same");
    fitMassSig->Draw("same");
    fitMassBg->Draw("same");

    // TString sLabel = Form("#splitline}{Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}",sPtLow.Data(),sPtHigh.Data(),sCent.Data());

    events.DrawLatexNDC(0.20,0.89,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (%s%%)", sCent.Data()));
    events.DrawLatexNDC(0.20,0.83,Form("%s < #it{p}_{T} < %s GeV/#it{c}",sPtLow.Data(),sPtHigh.Data()));
    species.DrawLatexNDC(0.84,0.87,sSpec);
    // events.DrawLatexNDC(0.5,0.88,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (%s%%)", sCent.Data()));
    // events.DrawLatexNDC(0.55,0.82,Form("%s < #it{p}_{T} < %s GeV/#it{c}",sPtLow.Data(),sPtHigh.Data()));
    // events.DrawLatexNDC(0.2,0.88,sSpec);

    canMass->SaveAs(Form("%s/fits/mass/single_nice/histMass_%s.%s",sPath.Data(),sName.Data(),gOutFormat.Data()),gOutFormat.Data());

     // MAKING RATIO PLOT WITH FRACTIONS

     TLine lineZero = TLine();
     lineZero.SetLineColor(kBlack);
     lineZero.SetLineWidth(1);


    TCanvas* canDiv = new TCanvas("canDiv","canDiv",600,600);
    canDiv->cd();
    TPad* padRatio = new TPad("padRatio","padRatio", 0, 0.0, 1, 0.3);
    padRatio->SetTopMargin(0.03);
    padRatio->SetBottomMargin(0.27);
    padRatio->SetRightMargin(0.05);
    padRatio->SetLeftMargin(0.15);
    padRatio->Draw();
    TPad* padMain = new TPad("padMain","padMain", 0, 0.3, 1, 1.0);
    padMain->SetBottomMargin(0.01);
    padMain->SetRightMargin(0.05);
    padMain->SetLeftMargin(0.15);
    padMain->SetTopMargin(0.07);
    padMain->Draw();

    histMass->SetMinimum(-0.05*dCorMax);
    histMass->SetMaximum(1.3*dCorMax);

    TH1D* histFracDummy = (TH1D*) histFracSig->Clone("dummy");
    histFracDummy->Reset();
    histFracDummy->SetTitle("");
    histFracDummy->SetMinimum(0.0);
    histFracDummy->SetMaximum(1.0);

    TAxis* axisFracX = histFracDummy->GetXaxis();
    TAxis* axisFracY = histFracDummy->GetYaxis();
    // axis_y_main->SetLabelSize(0.0);
    // axisFracY->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    // axisFracY->SetLabelSize(15);
    axisFracY->SetNdivisions(5,5,0,kTRUE);
    // axisFracX->SetNdivisions(510);
    axisFracY->SetTitle("Fraction     ");
    axisFracY->SetNdivisions(505);
    axisFracY->SetTitleSize(20);
    axisFracY->SetTitleFont(43);
    axisFracY->SetTitleOffset(1.6);
    axisFracY->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axisFracY->SetLabelSize(18);
    axisFracY->SetLabelOffset(0.01);

    // X axis ratio plot settings
    axisFracX->SetTitleSize(20);
    axisFracX->SetTitleFont(43);
    axisFracX->SetTitleOffset(3.);
    axisFracX->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axisFracX->SetLabelSize(18);
    axisFracX->SetTickSize(0.07);
    // axisFracX->SetTick(2);


    TAxis* axisMassX = histMass->GetXaxis();
    TAxis* axisMassY = histMass->GetYaxis();
    // axisMassY->SetNdivisions(5,5,0,kTRUE);
    axisMassY->SetTitle("Counts");
    axisMassY->SetTitleSize(20);
    axisMassY->SetTitleFont(43);
    axisMassY->SetTitleOffset(1.6);
    axisMassY->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axisMassY->SetLabelSize(17);
    axisMassY->SetLabelOffset(0.01);
    axisMassX->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    axisMassX->SetLabelSize(0.0);
    // axisMassX->SetLabelOffset(0.1);
    // axisMassX->SetTitleOffset(1.6);
    axisMassX->SetTickSize(0.03);
    axisMassY->SetMaxDigits(3);


    padMain->cd();
    histMass->DrawCopy();
    lineZero.DrawLine(dXlow,0.0,dXhigh,0.0);
    fitMass->Draw("same");
    fitMassSig->Draw("same");
    fitMassBg->Draw("same");
    events.DrawLatexNDC(0.20,0.87,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (%s%%)", sCent.Data()));
    events.DrawLatexNDC(0.20,0.81,Form("%s < #it{p}_{T} < %s GeV/#it{c}",sPtLow.Data(),sPtHigh.Data()));
    species.DrawLatexNDC(0.84,0.85,sSpec);


    padRatio->cd();
    histFracDummy->DrawCopy();
    gPad->SetTickx(1);
    // histFracSig->DrawCopy();
    // histFracBg->DrawCopy("same");
    fitFracSig->Draw("same");
    fitFracBg->Draw("same");


    canDiv->SaveAs(Form("%s/fits/mass/single_nice/histMassFrac_%s.%s",sPath.Data(),sName.Data(),gOutFormat.Data()),gOutFormat.Data());

    return;
}

void SetPad()
{
    // gPad->SetBottomMargin(0.2);
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.15);
}

// void SetHistOff(TH1* hist)
// {
//     hist->GetXaxis()->SetNdivisions(505,1);
//
//     hist->SetLabelFont(43,"XY");
//     hist->SetLabelSize(12,"XY");
//
//     hist->SetTitleFont(43,"XY");
//     hist->SetTitleSize(12,"XY");
//
//     hist->SetTitleOffset(4.3,"X");
//     // hist->SetTitleOffset(4.3,"Y");
//     // hist->SetTitleOffset(4.3,"t");
//     // hist->SetTitleFont(43,"t");
//     // hist->SetTitleSize(15,"t");
//
//     // hist->SetTitleOffset(2.2,"XY");
// }

void SetFuncAtt(TF1* func, Color_t color, Style_t style, Int_t width)
{
    func->SetNpx(1000);
    func->SetLineColor(color);
    func->SetLineStyle(style);
    func->SetLineWidth(width);
}

void SetHistAtt(TH1* hist, Color_t color, Style_t markStyle, Double_t markSize)
{
    hist->SetStats(0);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markStyle);
    hist->SetMarkerSize(markSize);

    hist->SetLineColor(color);
    // hist->SetLineStyle(lineStyle);
    // hist->SetLineWidth(lineWidth);
}
