// Small library for transfering TH1* objects into TGraphErrors

// // Input file Handling
// TFile* gFile = nullptr;
//
// Bool_t InputFile(
//     TString sFile
// )
// {
//     TFile* f = TFile::Open(sFile.Data(), "READ");
//     if(!f) {
//         printf("InputFile :: File '%s' not found!\n",sFile.Data());
//         return kFALSE;
//     }
//
//     gFile = f;
//
//     return kTRUE;
// }

void WriteToFile(TFile* fileOut, TObject* obj);

void SetColor(TGraph* gr,Color_t col);
// void SetMarker(TGraph* gr,Int_t marker,Double_t size = 1.0);
void SetMarker(TGraph* gr,Int_t marker,Double_t size);
void SetMarker(TGraph* gr,Int_t marker,Double_t size,Color_t col);
void SetStyle(TGraph* gr,Int_t marker, Double_t size, Color_t col);
// void SetStyle(TGraph* gr,Int_t marker, Double_t size = -10.0, Color_t col = -10);

TH1* LoadHisto(TFile* file,TString sName);
TH1* LoadHisto(TString sFile,TString sName);

TGraphErrors* LoadHistoToGraph(
    TFile* file,
    TString sHisto,
    TString sNewName,
    TString sNewTitle,
    Int_t marker,
    Double_t size,
    Color_t col
    // TString sNewName = "",
    // TString sNewTitle = "",
    // Int_t marker = -10,
    // Double_t size = -10.0,
    // Color_t col = -10
);
TGraphErrors* LoadHistoToGraph(
    TString sFile,
    TString sHisto,
    TString sNewName,
    TString sNewTitle,
    Int_t marker,
    Double_t size,
    Color_t col
);
//     TString sNewName = "",
//     TString sNewTitle = "",
//     Int_t marker = -10,
//     Double_t size = -10.0,
//     Color_t col = -10
// );

void SetGraphError(TGraphErrors* gr, Double_t dErrX, Double_t dErrY);

// =============================================================================

void WriteToFile(
    TFile* fileOut,
    TObject* obj
)
{
    if(!fileOut) {
        printf("WriteToFile :: FileOut not valid!\n");
        return;
    }

    if(!obj) {
        printf("WriteToFile :: Object not valid!\n ");
    }

    fileOut->cd();
    obj->Write();
}

void SetColor(
    TGraph* gr,
    Color_t col
)
{
    if(col < 0) {
        return;
    }

    if(!gr) {
        printf("SetColor :: Graph not valid!\n");
        return;
    }

    gr->SetLineColor(col);
    gr->SetMarkerColor(col);
    gr->SetFillColor(col);
}


void SetMarker(
    TGraph* gr,
    Int_t marker,
    Double_t size = 1.0
)
{
    if(!gr) {
        printf("SetMarker :: Graph not valid!\n");
        return;
    }

    gr->SetMarkerStyle(marker);

    if(size > 0.0) {
        gr->SetMarkerSize(size);
    }
}

void SetMarker(
    TGraph* gr,
    Int_t marker,
    Double_t size,
    Color_t col
)
{
    if(!gr) {
        printf("SetMarker :: Graph not valid!\n");
        return;
    }

    SetMarker(gr,marker,size);
    gr->SetMarkerColor(col);
}

void SetStyle(
    TGraph* gr,
    Int_t marker,
    Double_t size = -10.0,
    Color_t col = -10
)
{
    if(!gr) {
        printf("SetStyle :: Graph not valid!\n");
        return;
    }

    SetColor(gr, col);
    SetMarker(gr,marker,size);
}

TH1* LoadHisto(
    TFile* file,
    TString sName
)
{
    TH1* histo = (TH1*) file->Get(sName.Data());
    if(!histo) {
        printf("LoadHisto :: Histo '%s' not found!\n",sName.Data());
        file->ls();
        return nullptr;
    }

    return histo;
}

TH1* LoadHisto(
    TString sFile,
    TString sName
)
{
    TFile* file = TFile::Open(sFile.Data(),"READ");
    if(!file) {
        printf("LoadHisto :: File '%s' not openned!\n",sFile.Data());
        return nullptr;
    }

    return LoadHisto(file,sName);
}

TGraphErrors* LoadHistoToGraph(
    TFile* file,
    TString sHisto,
    TString sNewName = "",
    TString sNewTitle = "",
    Int_t marker = -10,
    Double_t size = -10.0,
    Color_t col = -10
)
{
    TH1* h = LoadHisto(file,sHisto);
    if(!h) { return nullptr; }

    TGraphErrors* g = new TGraphErrors(h);

    if(!sNewName.IsNull()) {
        g->SetName(sNewName);
    } else {
        g->SetName(sHisto);
    }

    if(!sNewTitle.IsNull()) {
        g->SetTitle(sNewTitle);
    }

    SetStyle(g, marker, size, col);

    return g;
}


TGraphErrors* LoadHistoToGraph(
    TString sFile,
    TString sHisto,
    TString sNewName = "",
    TString sNewTitle = "",
    Int_t marker = -10,
    Double_t size = -10.0,
    Color_t col = -10
)
{
    TFile* file = TFile::Open(sFile.Data(),"READ");
    if(!file) {
        printf("LoadHistoToGraph :: File '%s' not openned!\n",sFile.Data());
        return nullptr;
    }

    return LoadHistoToGraph(file, sHisto, sNewName, sNewTitle, marker, size, col);

}

void SetGraphError(TGraphErrors* gr, Double_t dErrX, Double_t dErrY = -1.0)
{
    if(!gr) { printf("SetErrorX :: Graph does not exists!\n"); return; }
    if(dErrX < 0.0) { printf("SetErrorX :: Error negative (%f)!\n",dErrX); return; }

    Int_t pts = gr->GetN();
    for(Int_t p(0); p < pts; ++p) {
        if(dErrY < 0.0) {
            Double_t dErrY_old = gr->GetErrorY(p);

            gr->SetPointError(p, dErrX, dErrY_old);
        } else {
            gr->SetPointError(p, dErrX, dErrY);
        }
    }

    return;
}
