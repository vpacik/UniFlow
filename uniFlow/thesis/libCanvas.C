// Small library for handling canvases for figures

TLatex* gLatexCan = new TLatex();
TLatex* gLatexPad = new TLatex();
TLine* gLineZero = new TLine();

TLegendEntry* GetHeader(TLegend* leg);
void SetLegendStyle(TLegend* leg);


TList* PrepareCanvas1x1(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax);
TList* PrepareCanvas2x3(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax);
TList* PrepareCanvas3x2(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax);
TList* PrepareCanvas3x3(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax);

Bool_t SetCanvas(TCanvas* can, TGaxis* gaX, TGaxis* gaY);
Bool_t SetCanvasRelative(TCanvas* can, TGaxis* gaX, TGaxis* gaY);

// =============================================================================
TLegendEntry* GetHeader(TLegend* leg)
{
    if(!leg) { printf("GetHeader :: TLegend does not exits\n"); return nullptr; }
    return ((TLegendEntry*) leg->GetListOfPrimitives()->First());
}

void SetLegendStyle(TLegend* leg)
{
    if(!leg) { printf("SetLegendStyle :: TLegend does not exits\n"); return; }
    leg->SetFillColorAlpha(0,0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(18);

    TLegendEntry* header = GetHeader(leg);
    if(!header) return;
    header->SetTextFont(43);
    header->SetTextSize(18);
}

TList* PrepareCanvas1x1(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax)
{
    if(!SetCanvas(can,gaX,gaY)) { printf("PrepareCanvas1x1 :: Canvas setting failed!\n"); return nullptr; }

    gaX->SetLabelFont(42);
    gaX->SetLabelSize(0.035);
    gaX->SetLabelOffset(0.01);

    gaY->SetLabelFont(42);
    gaY->SetLabelSize(0.035);
    gaY->SetLabelOffset(0.01);

    gaX->CenterTitle(0);
    gaX->SetTitleFont(42);
    gaX->SetTitleSize(0.04);
    gaX->SetTitleOffset(1.15);

    gaY->CenterTitle(0);
    gaY->SetTitleFont(42);
    gaY->SetTitleSize(0.04);
    gaY->SetTitleOffset(1.5);

    can->SetLeftMargin(0.13);
    can->SetRightMargin(0.03);
    can->SetTopMargin(0.06);
    can->SetBottomMargin(0.1);


    TList* list = new TList();

    const Int_t iNumPads = 1;

    for(Int_t i(0); i < iNumPads; i++) {
        can->cd(i+1);

        gPad->SetTicks(1,1);

        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);

        TAxis* ax = frame->GetXaxis();
        TAxis* ay = frame->GetYaxis();

        ax->SetLabelSize(0);
        ax->SetTitleSize(0);

        ay->SetLabelSize(0);
        ay->SetTitleSize(0);

        list->Add(frame);

        if(ymin < 0) { gLineZero->DrawLine(xmin, 0.0, xmax, 0.0); }
    }

    can->cd(1);
    gaY->Draw();
    gaX->Draw();

    return list;
}

TList* PrepareCanvas2x3(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax)
{
    if(!SetCanvas(can,gaX,gaY)) { printf("PrepareCanvas3x2 :: Canvas setting failed!\n"); return nullptr; }

    gaX->SetTitleOffset(3.5);
    gaY->SetTitleOffset(4.0);
    can->SetTopMargin(0.13);

    gaX->SetLabelSize(18);
    gaY->SetLabelSize(18);

    can->Divide(2,3,0.0,0.0);

    TList* list = new TList();

    const Int_t iNumPads = 6;

    for(Int_t i(0); i < iNumPads; i++) {
        can->cd(i+1);

        gPad->SetTicks(1,1);

        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);

        TAxis* ax = frame->GetXaxis();
        TAxis* ay = frame->GetYaxis();

        ax->SetLabelSize(0);
        ax->SetTitleSize(0);

        ay->SetLabelSize(0);
        ay->SetTitleSize(0);

        list->Add(frame);

        // Drawing lines at zero
        if(ymin < 0) { gLineZero->DrawLine(xmin, 0.0, xmax, 0.0); }
    }

    // Drawing axes
    can->cd(1);
    gaY->Draw();
    can->cd(3);
    gaY->Draw();
    can->cd(5);
    gaY->Draw();

    // can->cd(4);
    // gaX->Draw();
    can->cd(5);
    gaX->Draw();
    can->cd(6);
    gaX->Draw();

    return list;
}



TList* PrepareCanvas3x2(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax)
{
    if(!SetCanvas(can,gaX,gaY)) { printf("PrepareCanvas3x2 :: Canvas setting failed!\n"); return nullptr; }

    gaX->SetTitleOffset(2.2);
    gaY->SetTitleOffset(3.0);
    can->SetTopMargin(0.13);

    gaX->SetLabelSize(18);
    gaY->SetLabelSize(18);

    can->Divide(3,2,0.0,0.0);

    TList* list = new TList();

    const Int_t iNumPads = 6;

    for(Int_t i(0); i < iNumPads; i++) {
        can->cd(i+1);

        gPad->SetTicks(1,1);

        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);

        TAxis* ax = frame->GetXaxis();
        TAxis* ay = frame->GetYaxis();

        ax->SetLabelSize(0);
        ax->SetTitleSize(0);

        ay->SetLabelSize(0);
        ay->SetTitleSize(0);

        list->Add(frame);

        // Drawing lines at zero
        if(ymin < 0) { gLineZero->DrawLine(xmin, 0.0, xmax, 0.0); }
    }

    // Drawing axes
    can->cd(1);
    gaY->Draw();
    can->cd(4);
    gaY->Draw();
    // can->cd(7);
    // gaY->Draw();

    can->cd(4);
    gaX->Draw();
    can->cd(5);
    gaX->Draw();
    can->cd(6);
    gaX->Draw();

    return list;
}

TList* PrepareCanvas3x3(TCanvas* can, TGaxis* gaX, TGaxis* gaY, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax)
{
    if(!SetCanvas(can,gaX,gaY)) { printf("PrepareCanvas3x3 :: Canvas setting failed!\n"); return nullptr; }

    can->Divide(3,3,0.0,0.0);

    TList* list = new TList();

    const Int_t iNumPads = 9;


    for(Int_t i(0); i < iNumPads; i++) {
        can->cd(i+1);

        gPad->SetTicks(1,1);

        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);

        TAxis* ax = frame->GetXaxis();
        TAxis* ay = frame->GetYaxis();

        ax->SetLabelSize(0);
        ax->SetTitleSize(0);

        ay->SetLabelSize(0);
        ay->SetTitleSize(0);

        list->Add(frame);

        if(ymin < 0) { gLineZero->DrawLine(xmin, 0.0, xmax, 0.0); }
    }

    // Drawing axes

    can->cd(1);
    gaY->Draw();
    can->cd(4);
    gaY->Draw();
    can->cd(7);
    gaY->Draw();

    can->cd(7);
    gaX->Draw();
    can->cd(8);
    gaX->Draw();
    can->cd(9);
    gaX->Draw();

    return list;
}

Bool_t SetCanvas(TCanvas* can, TGaxis* gaX, TGaxis* gaY)
{
    if(!can) { printf("SetCanvas :: No canvas found!\n"); return kFALSE; }
    if(!gaX) { printf("SetCanvas :: gaX found!\n"); return kFALSE; }
    if(!gaY) { printf("SetCanvas :: gaY found!\n"); return kFALSE; }

    gaX->SetLabelFont(43);
    gaX->SetLabelSize(14);
    gaX->SetLabelOffset(0.01);

    gaX->CenterTitle();
    gaX->SetTitleFont(43);
    gaX->SetTitleSize(18);
    gaX->SetTitleOffset(2.8);

    gaY->SetLabelFont(43);
    gaY->SetLabelSize(14);
    gaY->SetLabelOffset(0.01);

    gaY->CenterTitle();
    gaY->SetTitleFont(43);
    gaY->SetTitleSize(18);
    gaY->SetTitleOffset(3.5);

    can->SetLeftMargin(0.2);
    can->SetRightMargin(0.1);
    can->SetTopMargin(0.15);
    can->SetBottomMargin(0.18);

    if(!gLatexCan) { gLatexCan = new TLatex(); }
    if(!gLatexPad) { gLatexPad = new TLatex(); }

    gLatexCan->SetTextFont(43);
    gLatexCan->SetTextSize(20);
    gLatexCan->SetTextAlign(11);

    gLatexPad->SetTextFont(43);
    gLatexPad->SetTextSize(18);
    gLatexPad->SetTextAlign(33);

    if(!gLineZero) { gLineZero = new TLine(); }
    gLineZero->SetLineColor(kGray+1);
    gLineZero->SetLineStyle(kDashed);

    return kTRUE;
}

Bool_t SetCanvasRelative(TCanvas* can, TGaxis* gaX, TGaxis* gaY)
{
    if(!can) { printf("SetCanvas :: No canvas found!\n"); return kFALSE; }
    if(!gaX) { printf("SetCanvas :: gaX found!\n"); return kFALSE; }
    if(!gaY) { printf("SetCanvas :: gaY found!\n"); return kFALSE; }

    gaX->SetLabelFont(42);
    gaX->SetLabelSize(0.07);
    gaX->SetLabelOffset(0.01);

    gaX->CenterTitle();
    gaX->SetTitleFont(42);
    gaX->SetTitleSize(0.075);
    gaX->SetTitleOffset(1.13);

    gaY->SetLabelFont(42);
    gaY->SetLabelSize(0.07);
    gaY->SetLabelOffset(0.01);

    gaY->CenterTitle();
    gaY->SetTitleFont(42);
    gaY->SetTitleSize(0.75);
    gaY->SetTitleOffset(1.13);

    can->SetLeftMargin(0.2);
    can->SetRightMargin(0.1);
    can->SetTopMargin(0.15);
    can->SetBottomMargin(0.18);

    return kTRUE;
}
