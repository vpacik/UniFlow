TGraphErrors* Scale(TGraphErrors* orig,Double_t dFactorX, Double_t dFactorY);
TGraphErrors* TransKEt(TGraphErrors* orig, Double_t dMass);


void ScaleVn()
{
    TFile* fileIn = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/results/PbPb/cums/v24_prel/graphs_PbPb_v22_v24.root","READ");
    if(!fileIn) { printf("File in not open!\n"); return; }

    TFile* fileOut = TFile::Open("/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/results/PbPb/cums/v24_prel/graphs_PbPb_v22_v24_scaled.root","RECREATE");
    if(!fileOut) { printf("File out not open!\n"); return; }

    fileIn->ls();

    const Int_t iHistos = 5;
    TString sHistoName[iHistos] = {"v2{2}", "v2{2,gap00}", "v2{2,gap08}","v2{4}", "v2{4,gap00}"};

    const Int_t iSpecies = 6;
    TString sSpecies[iSpecies] = { "Pion", "Kaon", "K0s", "Proton", "Lambda", "Phi" };
    Double_t dScale[iSpecies] = { 1./2, 1./2, 1./2, 1./3, 1./3, 1./2 };
    Double_t dMass[iSpecies] = { 0.13957000, 0.49367700, 0.49761400, 0.93827200, 1.1156800, 1.0194550 };

    const Int_t iCent = 6;
    TString sCent[iCent] = {"0510","1020","2030","3040","4050","5060"};

    for(Int_t s(0); s < iSpecies; ++s)
    {
        for(Int_t h(0); h < iHistos; ++h)
        {
            for(Int_t c(0); c < iCent; ++c)
            {
                TString sName = Form("%s_%s_cent%s", sSpecies[s].Data(), sHistoName[h].Data(), sCent[c].Data());

                TGraphErrors* gr = (TGraphErrors*) fileIn->Get(sName.Data());

                // ncq-like scaling
                TGraphErrors* scaled = Scale(gr,dScale[s],dScale[s]);
                scaled->SetName(Form("%s_NCQ",gr->GetName()));
                scaled->SetTitle(Form("%s (v_{n}/n_{q})",gr->GetTitle()));
                // scaled->SetLineColor(kBlue);
                // scaled->Draw("pel");

                // KEt-like scaling
                TGraphErrors* trans = Scale(TransKEt(gr,dMass[s]),dScale[s],dScale[s]);
                trans->SetName(Form("%s_KEt",gr->GetName()));
                trans->SetTitle(Form("%s (KE_{T}/n_{q})",gr->GetTitle()));
                // trans->SetLineColor(kGreen);
                // trans->Draw("pel");

                // writing to output file
                fileOut->cd();
                scaled->Write();
                trans->Write();
            }
        }
    }

    fileOut->ls();

    return;
}

TGraphErrors* Scale(TGraphErrors* orig, Double_t dFactorX, Double_t dFactorY)
{
    if(!orig) { printf("Scale :: Input graph not valid!\n"); return nullptr; }

    TGraphErrors* scaled = (TGraphErrors*) orig->Clone(Form("scaled"));

    Int_t points = scaled->GetN();
    // printf("points %d\n", points);

    for(Int_t p(0); p < points; ++p)
    {
        Double_t x,y;
        orig->GetPoint(p,x,y);

        Double_t ex = orig->GetErrorX(p);
        Double_t ey = orig->GetErrorY(p);

        // printf("%d : %f,%f :: %f,%f\n", p, x, y, ex, ey);

        Double_t dSX = x * dFactorX;
        Double_t dSY = y * dFactorY;

        Double_t dSEX = TMath::Sqrt(dFactorX*dFactorX * ex*ex);
        Double_t dSEY = TMath::Sqrt(dFactorY*dFactorY * ey*ey);

        scaled->SetPoint(p,dSX,dSY);
        scaled->SetPointError(p,dSEX,dSEY);
    }

    return scaled;
}

TGraphErrors* TransKEt(TGraphErrors* orig, Double_t dMass)
{

        TGraphErrors* trans = (TGraphErrors*) orig->Clone(Form("trans"));

        Int_t points = trans->GetN();
        // printf("points %d\n", points);

        for(Int_t p(0); p < points; ++p)
        {
            Double_t x,y;
            orig->GetPoint(p,x,y);

            Double_t ex = orig->GetErrorX(p);
            Double_t ey = orig->GetErrorY(p);

            // printf("%d : %f,%f :: %f,%f\n", p, x, y, ex, ey);

            Double_t dSX = TMath::Sqrt(x*x + dMass*dMass) - dMass;

            Double_t dSEX = TMath::Sqrt(x*x * ex*ex / (x*x + dMass*dMass));

            trans->SetPoint(p,dSX,y);
            trans->SetPointError(p,dSEX,ey);
        }

        return trans;

}
