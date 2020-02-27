#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStyle.C"

void Trans_FlowFluct()
{
    TString sInputFile = "./process/FlowFluctCor.root";
    TString sOutputFile = "./graphs_PbPb_FlowFluctCor.root";

    constexpr Int_t iCent = 6;
    TString sCent[iCent] = {"0010","1020","2030","3040","4050","5060"};
    TString sCentLab[iCent] = {"0-10","10-20","20-30","30-40","40-50","50-60"};

    SetStyleCumsPbPb();

    TFile* fileIn = TFile::Open(sInputFile.Data(), "READ");
    if(!fileIn) { printf("ToGraphs :: fileIn '%s' not open!\n",sInputFile.Data()); return nullptr; }

    TFile* fileOut = TFile::Open(sOutputFile.Data(), "RECREATE");
    if(!fileOut) { printf("ToGraphs :: fileOut '%s' not open!\n",sOutputFile.Data()); return nullptr; }

    // pt diff

    for(Int_t s(0); s < kRefs; ++s)
    {
        if(s == kSpecies(kPhi)) break;
        for(Int_t c(0); c < iCent; ++c)
        {
            WriteToFile(fileOut,
                LoadHistoToGraph(fileIn,
                    Form("%s_mean_cent%d",sSpecies[s].Data(),c),
                    Form("%s_mean_v2_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                    Form("%s: #langle v_{2} #rangle (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                    mark[s], size[s], col[s]));

            WriteToFile(fileOut,
                LoadHistoToGraph(fileIn,
                    Form("%s_std_cent%d",sSpecies[s].Data(),c),
                    Form("%s_std_v2_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                    Form("%s: #sigma_{2} (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                    mark[s], size[s], col[s]));

            WriteToFile(fileOut,
                LoadHistoToGraph(fileIn,
                    Form("%s_rel_cent%d",sSpecies[s].Data(),c),
                    Form("%s_Fv2_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                    Form("%s: F(v_{2}) (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                    mark[s], size[s], col[s]));


            // WriteToFile(fileOut, LoadHistoToGraph(fileIn, "Charged_hFlow2_harm2_gap00_c0", Form("Charged_v2{2,gap00}_c0510", "h^{#pm}: v_{2}{2,|#Delta#eta| > 0} (5-10% V0M)", mark[kCharged], size[kCharged], col[kCharged]));

        }

    }

    fileOut->ls();

    return;

    // Alternative method writing to file
    // {
    //     TGraphErrors* graph = LoadHistoToGraph(fileIn, "Charged_hFlow2_harm2_gap00_c0", "", "h^{#pm}: v_{2}{2,|#Delta#eta| > 0} (0-20% V0M)", kFullSquare, 1.5, kRed);
    //     WriteToFile(fileOut,graph);
    //
    //     TGraphErrors* graph2 = LoadHistoToGraph(fileIn, "Charged_hFlow2_harm2_gap00_c1", "", "h^{#pm}: v_{2}{2,|#Delta#eta| > 0} (20-40% V0M)", kFullSquare, 1.5, kBlack);
    //     WriteToFile(fileOut,graph2);
    //
    //     fileOut->cd();
    //     graph->Write();
    //     graph2->Write();
    // }
}
