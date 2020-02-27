// Macro for plotting fits within fits.root file from ProcessUniFlow
// 2019 for thesis

#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"

Double_t markSize = 1.1;
Int_t lineSize = 2;


void SetPad();
void SetHistOff(TH1* hist);
void SetFuncAtt(TF1* func, Color_t color, Style_t style, Int_t width = 2.0);
void SetHistAtt(TH1* hist, Color_t color, Style_t markStyle, Double_t markSize = 1.0);

void PlotFitsSingleCorr(
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/6815/thesis/",
    TString sCorrName = "<<2>>(2,-2)",
    TString sCorrTitle = "<<2>>_(2,-2)",
    TString sSpecies = "K0s",
    Int_t iPt = 8,
    Int_t iCent = 5)
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

    TH1D* histCorr = (TH1D*) list->FindObject("histCorr");
    if(!histCorr) { printf("ERROR: Object '%s' not found!\n", "histCorr"); list->ls(); return kFALSE; }

    TF1* fitCorr = (TF1*) list->FindObject("fitCorr");
    if(!fitCorr) { printf("ERROR: Object '%s' not found!\n", "fitCorr"); list->ls(); return kFALSE; }

    TF1* fitCorrSig = (TF1*) list->FindObject("fitCorSig");
    if(!fitCorrSig) { printf("ERROR: Object '%s' not found!\n", "fitCorrSig"); list->ls(); return kFALSE; }

    TF1* fitCorrBg = (TF1*) list->FindObject("fitCorBg");
    if(!fitCorrBg) { printf("ERROR: Object '%s' not found!\n", "fitCorrBg"); list->ls(); return kFALSE; }

    TF1* fitCorrBgCoef = (TF1*) list->FindObject("fitFlowBg");
    if(!fitCorrBgCoef) { printf("ERROR: Object '%s' not found!\n", "fitCorrBgCoef"); list->ls(); return kFALSE; }

    Double_t dXlow = fitCorr->GetXmin();
    Double_t dXhigh = fitCorr->GetXmax();

    // histCorr->SetTitle(sTitle.Data());
    histCorr->SetTitle("");
    histCorr->GetYaxis()->SetTitle(sCorrTitle.Data());

    histCorr->GetXaxis()->SetRangeUser(dXlow,dXhigh);
    histCorr->SetStats(0);

    Double_t dCorMax = histCorr->GetMaximum();

    histCorr->SetMinimum(-0.05*dCorMax);
    histCorr->SetMaximum(1.4*dCorMax);

    SetHistAtt(histCorr, kBlack, kFullCircle, 0.7);
    SetFuncAtt(fitCorr, kRed, kSolid, 3);
    SetFuncAtt(fitCorrSig, kGreen-2, 7, 3);
    SetFuncAtt(fitCorrBg, kBlue, 2, 3);
    SetFuncAtt(fitCorrBgCoef, kOrange-6, 5, 3);

    TLatex events;
    // events.SetNDC(1);
    events.SetTextFont(42);
    events.SetTextSize(0.04);
    // events.SeteventsAlign(30);

    TLatex species;
    species.SetTextFont(42);
    species.SetTextSize(0.06);

    gSystem->mkdir(Form("%s/fits/corr/single_nice/",sPath.Data()),1);

    TCanvas* canCorr = new TCanvas("canCorr","canCorr",600,600);
    canCorr->cd();
    SetPad();
    histCorr->DrawCopy();
    fitCorr->Draw("same");
    fitCorrSig->Draw("same");
    fitCorrBg->Draw("same");
    fitCorrBgCoef->Draw("same");

    // TString sLabel = Form("#splitline}{Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}",sPtLow.Data(),sPtHigh.Data(),sCent.Data());

    events.DrawLatexNDC(0.20,0.89,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (%s%%)", sCent.Data()));
    events.DrawLatexNDC(0.20,0.83,Form("%s < #it{p}_{T} < %s GeV/#it{c}",sPtLow.Data(),sPtHigh.Data()));
    species.DrawLatexNDC(0.84,0.87,sSpec);
    // events.DrawLatexNDC(0.5,0.88,Form("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (%s%%)", sCent.Data()));
    // events.DrawLatexNDC(0.55,0.82,Form("%s < #it{p}_{T} < %s GeV/#it{c}",sPtLow.Data(),sPtHigh.Data()));
    // events.DrawLatexNDC(0.2,0.88,sSpec);

    canCorr->SaveAs(Form("%s/fits/corr/single_nice/histCorr_%s.%s",sPath.Data(),sName.Data(),gOutFormat.Data()),gOutFormat.Data());

    // ========== Processing histFrac =======================================


    // TH1D* histMass_fracSig = (TH1D*) list->FindObject("histMass_fracSig");
    // if(!histMass_fracSig) { printf("ERROR: Object '%s' not found!\n", "histMass_fracSig"); list->ls(); return kFALSE; }
    //
    // TH1D* histMass_fracBg = (TH1D*) list->FindObject("histMass_fracBg");
    // if(!histMass_fracBg) { printf("ERROR: Object '%s' not found!\n", "histMass_fracBg"); list->ls(); return kFALSE; }
    //
    // TF1* fitMass_fracSig = (TF1*) list->FindObject("fitMass_fracSig");
    // if(!fitMass_fracSig) { printf("ERROR: Object '%s' not found!\n", "fitMass_fracSig"); list->ls(); return kFALSE; }
    //
    // TF1* fitMass_fracBg = (TF1*) list->FindObject("fitMass_fracBg");
    // if(!fitMass_fracBg) { printf("ERROR: Object '%s' not found!\n", "fitMass_fracBg"); list->ls(); return kFALSE; }
    //
    // // Double_t dXlow = fitCorr->GetXmin();
    // // Double_t dXhigh = fitCorr->GetXmax();
    //
    // histMass_fracSig->SetTitle(sTitle.Data());
    // histMass_fracSig->GetXaxis()->SetRangeUser(dXlow,dXhigh);
    // histMass_fracSig->SetStats(0);
    // histMass_fracSig->SetMinimum(-0.05);
    // histMass_fracSig->SetMaximum(1.05);
    //
    // SetHistAtt(histMass_fracSig, Color_t(colorSig+1), kFullCircle, markSize);
    // SetHistAtt(histMass_fracBg, Color_t(colorBg+2), kFullCircle, markSize);
    //
    // SetFuncAtt(fitMass_fracSig, colorSig, kSolid, lineSize);
    // SetFuncAtt(fitMass_fracBg, colorBg, kSolid, lineSize);
    //
    // gSystem->mkdir(Form("%s/fits/mass/single/",sPath.Data()),1);
    //


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
    func->SetLineColor(color);
    func->SetLineStyle(style);
    func->SetLineWidth(width);
}

void SetHistAtt(TH1* hist, Color_t color, Style_t markStyle, Double_t markSize)
{
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markStyle);
    hist->SetMarkerSize(markSize);

    hist->SetLineColor(color);
    // hist->SetLineStyle(lineStyle);
    // hist->SetLineWidth(lineWidth);
}
