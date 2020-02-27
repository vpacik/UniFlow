
enum kSpecies { kCharged = 0, kPion, kKaon, kProton, kK0s, kLambda, kPhi, kRefs, kNum};
Color_t gCol[kNum];
Int_t gMark[kNum];
Double_t gSize[kNum];
TString gSpecies[kNum];
TString gSpeciesLabels[kNum];


void SetMarkerPID(TGraph* gr, kSpecies sp);

void SetStylePID();
void SetStyleCumsPbPb();


void SetMarkerPID(TGraph* gr, kSpecies sp)
{
    if(!gr) { printf("SetMarkerPID :: Graph does not exists!\n"); return; }

    gr->SetMarkerStyle(gMark[sp]);
    gr->SetMarkerColor(gCol[sp]);
    gr->SetMarkerSize(gSize[sp]);
    gr->SetLineColor(gCol[sp]);

    return;
}

void SetStylePID()
{
    gSpecies[kCharged] = "Charged";
    gSpecies[kPion] = "Pion";
    gSpecies[kKaon] = "Kaon";
    gSpecies[kK0s] = "K0s";
    gSpecies[kProton] = "Proton";
    gSpecies[kLambda] = "Lambda";
    gSpecies[kPhi] = "Phi";
    gSpecies[kRefs] = "Refs";

    gSpeciesLabels[kCharged] = "h^{#pm}";
    gSpeciesLabels[kPion] = "#pi^{#pm}";
    gSpeciesLabels[kKaon] = "K^{#pm}";
    gSpeciesLabels[kK0s] = "K_{S}^{0}";
    gSpeciesLabels[kProton] = "p/#bar{p}";
    gSpeciesLabels[kLambda] = "#Lambda/#bar{#Lambda}";
    gSpeciesLabels[kPhi] = "#phi";
    gSpeciesLabels[kRefs] = "RFP";

    gMark[kCharged] = kOpenSquare;
    gMark[kPion] = kFullCircle;
    gMark[kKaon] = kFullTriangleUp;
    gMark[kK0s] = kFullTriangleDown;
    gMark[kProton] = kFullSquare;
    gMark[kLambda] = kFullDiamond;
    gMark[kPhi] = kFullCross;
    gMark[kRefs] = kFullCircle;

    gCol[kCharged] = kGray+2;
    gCol[kPion] = kRed+1;
    gCol[kKaon] = kGreen+3;
    gCol[kK0s] = kGreen+1;
    gCol[kProton] = kBlue;
    gCol[kLambda] = kOrange+1;
    gCol[kPhi] = kMagenta;
    gCol[kRefs] = kBlack;

    gSize[kCharged] = 0.8;
    gSize[kPion] = 0.8;
    gSize[kKaon] = 0.8;
    gSize[kK0s] = 0.8;
    gSize[kProton] = 0.7;
    gSize[kLambda] = 1.1;
    gSize[kPhi] = 1.0;
    gSize[kRefs] = 1.0;

    // gSpecies.push_back("Charged"); sSpesLabels.push_back(" h^{#pm}"); gColors.push_back(kGray+2); gMarkers.push_back(kOpenCircle); gMarkergSizes.push_back(0.8); systs.push_back(0.03);
    // gSpecies.push_back("Pion"); sSpesLabels.push_back(" #pi^{#pm}"); gColors.push_back(kRed+1); gMarkers.push_back(kFullCircle); gMarkergSizes.push_back(0.8); systs.push_back(0.03);
    // gSpecies.push_back("Kaon"); sSpesLabels.push_back(" K^{#pm}"); gColors.push_back(kGreen+3); gMarkers.push_back(kFullTriangleUp); gMarkergSizes.push_back(1); systs.push_back(0.04);
    // gSpecies.push_back("K0s"); sSpesLabels.push_back(" K^{0}_{S}"); gColors.push_back(kGreen+1); gMarkers.push_back(kFullTriangleDown); gMarkergSizes.push_back(1); systs.push_back(0.07);
    // gSpecies.push_back("Proton"); sSpesLabels.push_back( " p(#bar{p})"); gColors.push_back(kBlue); gMarkers.push_back(kFullSquare); gMarkergSizes.push_back(0.7); systs.push_back(0.04);
    // gSpecies.push_back("Lambda"); sSpesLabels.push_back(" #Lambda(#bar{#Lambda})"); gColors.push_back(kOrange+1); gMarkers.push_back(kFullDiamond); gMarkergSizes.push_back(1.4); systs.push_back(0.05);
    // if(bPhi) { gSpecies.push_back("Phi"); sSpesLabels.push_back(" #phi"); gColors.push_back(kMagenta); gMarkers.push_back(kFullCross); gMarkergSizes.push_back(1.); systs.push_back(0.08); }
    //

}




void SetStyleCumsPbPb()
{
    gSpecies[kCharged] = "Charged";
    gSpecies[kPion] = "Pion";
    gSpecies[kKaon] = "Kaon";
    gSpecies[kK0s] = "K0s";
    gSpecies[kProton] = "Proton";
    gSpecies[kLambda] = "Lambda";
    gSpecies[kPhi] = "Phi";
    gSpecies[kRefs] = "Refs";

    gSpeciesLabels[kCharged] = "h^{#pm}";
    gSpeciesLabels[kPion] = "#pi^{#pm}";
    gSpeciesLabels[kKaon] = "K^{#pm}";
    gSpeciesLabels[kK0s] = "K_{S}^{0}";
    gSpeciesLabels[kProton] = "p/#bar{p}";
    gSpeciesLabels[kLambda] = "#Lambda/#bar{#Lambda}";
    gSpeciesLabels[kPhi] = "#phi";
    gSpeciesLabels[kRefs] = "RFP";

    gMark[kCharged] = kOpenCircle;
    gMark[kPion] = kFullCircle;
    gMark[kKaon] = kFullTriangleUp;
    gMark[kK0s] = kFullTriangleDown;
    gMark[kProton] = kFullSquare;
    gMark[kLambda] = kFullDiamond;
    gMark[kPhi] = kFullCross;
    gMark[kRefs] = kFullCircle;

    gCol[kCharged] = kGray+2;
    gCol[kPion] = kRed+1;
    gCol[kKaon] = kGreen+3;
    gCol[kK0s] = kGreen+1;
    gCol[kProton] = kBlue;
    gCol[kLambda] = kOrange+1;
    gCol[kPhi] = kMagenta;
    gCol[kRefs] = kBlack;

    gSize[kCharged] = 0.8;
    gSize[kPion] = 0.8;
    gSize[kKaon] = 1.0;
    gSize[kK0s] = 1.0;
    gSize[kProton] = 0.7;
    gSize[kLambda] = 1.4;
    gSize[kPhi] = 1.0;
    gSize[kRefs] = 1.0;

    // gSpecies.push_back("Charged"); sSpesLabels.push_back(" h^{#pm}"); gColors.push_back(kGray+2); gMarkers.push_back(kOpenCircle); gMarkergSizes.push_back(0.8); systs.push_back(0.03);
    // gSpecies.push_back("Pion"); sSpesLabels.push_back(" #pi^{#pm}"); gColors.push_back(kRed+1); gMarkers.push_back(kFullCircle); gMarkergSizes.push_back(0.8); systs.push_back(0.03);
    // gSpecies.push_back("Kaon"); sSpesLabels.push_back(" K^{#pm}"); gColors.push_back(kGreen+3); gMarkers.push_back(kFullTriangleUp); gMarkergSizes.push_back(1); systs.push_back(0.04);
    // gSpecies.push_back("K0s"); sSpesLabels.push_back(" K^{0}_{S}"); gColors.push_back(kGreen+1); gMarkers.push_back(kFullTriangleDown); gMarkergSizes.push_back(1); systs.push_back(0.07);
    // gSpecies.push_back("Proton"); sSpesLabels.push_back( " p(#bar{p})"); gColors.push_back(kBlue); gMarkers.push_back(kFullSquare); gMarkergSizes.push_back(0.7); systs.push_back(0.04);
    // gSpecies.push_back("Lambda"); sSpesLabels.push_back(" #Lambda(#bar{#Lambda})"); gColors.push_back(kOrange+1); gMarkers.push_back(kFullDiamond); gMarkergSizes.push_back(1.4); systs.push_back(0.05);
    // if(bPhi) { gSpecies.push_back("Phi"); sSpesLabels.push_back(" #phi"); gColors.push_back(kMagenta); gMarkers.push_back(kFullCross); gMarkergSizes.push_back(1.); systs.push_back(0.08); }
    //

}
