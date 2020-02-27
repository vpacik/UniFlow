
enum kSpecies { kCharged = 0, kPion, kKaon, kProton, kK0s, kLambda, kPhi, kRefs, kNum};
Color_t col[kNum];
Int_t mark[kNum];
Double_t size[kNum];
TString sSpecies[kNum];
TString sSpecLab[kNum];

void SetStyleCumsPbPb();


void SetStyleCumsPbPb()
{
    // enum kSpecies { kCharged = 0, kPion, kKaon, kProton, kK0s, kLambda, kPhi, kRefs, kNum};
    // Color_t col[kNum];
    // Int_t mark[kNum];
    // Double_t size[kNum];
    // TString sSpecies[kNum];
    // TString sSpecLab[kNum];



    sSpecies[kCharged] = "Charged";
    sSpecies[kPion] = "Pion";
    sSpecies[kKaon] = "Kaon";
    sSpecies[kK0s] = "K0s";
    sSpecies[kProton] = "Proton";
    sSpecies[kLambda] = "Lambda";
    sSpecies[kPhi] = "Phi";
    sSpecies[kRefs] = "Refs";

    sSpecLab[kCharged] = "h^{#pm}";
    sSpecLab[kPion] = "#pi^{#pm}";
    sSpecLab[kKaon] = "K^{#pm}";
    sSpecLab[kK0s] = "K_{S}^{0}";
    sSpecLab[kProton] = "p/#bar{p}";
    sSpecLab[kLambda] = "#Lambda/#bar{#Lambda}";
    sSpecLab[kPhi] = "#phi";
    sSpecLab[kRefs] = "RFP";

    mark[kCharged] = kOpenCircle;
    mark[kPion] = kFullCircle;
    mark[kKaon] = kFullTriangleUp;
    mark[kK0s] = kFullTriangleDown;
    mark[kProton] = kFullSquare;
    mark[kLambda] = kFullDiamond;
    mark[kPhi] = kFullCross;
    mark[kRefs] = kFullCircle;

    col[kCharged] = kGray+2;
    col[kPion] = kRed+1;
    col[kKaon] = kGreen+3;
    col[kK0s] = kGreen+1;
    col[kProton] = kBlue;
    col[kLambda] = kOrange+1;
    col[kPhi] = kMagenta;
    col[kRefs] = kBlack;

    size[kCharged] = 0.8;
    size[kPion] = 0.8;
    size[kKaon] = 1.0;
    size[kK0s] = 1.0;
    size[kProton] = 0.7;
    size[kLambda] = 1.4;
    size[kPhi] = 1.0;
    size[kRefs] = 1.0;

    //
    // sSpecies.push_back("Charged"); sSpesLabels.push_back(" h^{#pm}"); colors.push_back(kGray+2); markers.push_back(kOpenCircle); markerSizes.push_back(0.8); systs.push_back(0.03);
    // sSpecies.push_back("Pion"); sSpesLabels.push_back(" #pi^{#pm}"); colors.push_back(kRed+1); markers.push_back(kFullCircle); markerSizes.push_back(0.8); systs.push_back(0.03);
    // sSpecies.push_back("Kaon"); sSpesLabels.push_back(" K^{#pm}"); colors.push_back(kGreen+3); markers.push_back(kFullTriangleUp); markerSizes.push_back(1); systs.push_back(0.04);
    // sSpecies.push_back("K0s"); sSpesLabels.push_back(" K^{0}_{S}"); colors.push_back(kGreen+1); markers.push_back(kFullTriangleDown); markerSizes.push_back(1); systs.push_back(0.07);
    // sSpecies.push_back("Proton"); sSpesLabels.push_back( " p(#bar{p})"); colors.push_back(kBlue); markers.push_back(kFullSquare); markerSizes.push_back(0.7); systs.push_back(0.04);
    // sSpecies.push_back("Lambda"); sSpesLabels.push_back(" #Lambda(#bar{#Lambda})"); colors.push_back(kOrange+1); markers.push_back(kFullDiamond); markerSizes.push_back(1.4); systs.push_back(0.05);
    // if(bPhi) { sSpecies.push_back("Phi"); sSpesLabels.push_back(" #phi"); colors.push_back(kMagenta); markers.push_back(kFullCross); markerSizes.push_back(1.); systs.push_back(0.08); }
    //

}
