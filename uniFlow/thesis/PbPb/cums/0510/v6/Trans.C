#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/macros/libGraphs.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/macros/libStyle.C"

void Trans()
{
    TString sInputFile = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/results/PbPb/cums/v6/process/Processed.root";
    TString sOutputFile = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/results/PbPb/cums/v6/graphs_PbPb_v6.root";

    constexpr Int_t iCent = 6;
    TString sCent[iCent] = {"0510","1020","2030","3040","4050","5060"};
    TString sCentLab[iCent] = {"05-10","10-20","20-30","30-40","40-50","50-60"};

    SetStyleCumsPbPb();

    TFile* fileIn = TFile::Open(sInputFile.Data(), "READ");
    if(!fileIn) { printf("ToGraphs :: fileIn '%s' not open!\n",sInputFile.Data()); return nullptr; }

    TFile* fileOut = TFile::Open(sOutputFile.Data(), "RECREATE");
    if(!fileOut) { printf("ToGraphs :: fileOut '%s' not open!\n",sOutputFile.Data()); return nullptr; }


    // refs

    WriteToFile(fileOut,
        LoadHistoToGraph(fileIn,
            Form("%s_hFlow2_harm3_gap-10",sSpecies[kRefs].Data()),
            Form("%s_v3{2}",sSpecies[kRefs].Data()),
            Form("%s: v_{3}{2}",sSpecLab[kRefs].Data()),
            mark[kRefs], size[kRefs], col[kRefs]));

    WriteToFile(fileOut,
        LoadHistoToGraph(fileIn,
            Form("%s_hFlow2_harm3_gap00",sSpecies[kRefs].Data()),
            Form("%s_v3{2,gap00}",sSpecies[kRefs].Data()),
            Form("%s: v_{3}{2,|#Delta#eta| > 0}",sSpecLab[kRefs].Data()),
            mark[kRefs], size[kRefs], col[kRefs]));

    WriteToFile(fileOut,
        LoadHistoToGraph(fileIn,
            Form("%s_hFlow2_harm3_gap08",sSpecies[kRefs].Data()),
            Form("%s_v3{2,gap08}",sSpecies[kRefs].Data()),
            Form("%s: v_{3}{2,|#Delta#eta| > 0.8}",sSpecLab[kRefs].Data()),
            mark[kRefs], size[kRefs], col[kRefs]));

    WriteToFile(fileOut,
        LoadHistoToGraph(fileIn,
            Form("%s_hFlow4_harm3_gap-10",sSpecies[kRefs].Data()),
            Form("%s_v3{4}",sSpecies[kRefs].Data()),
            Form("%s: v_{3}{4}",sSpecLab[kRefs].Data()),
            mark[kRefs], size[kRefs], col[kRefs]));

    WriteToFile(fileOut,
        LoadHistoToGraph(fileIn,
            Form("%s_hFlow4_harm3_gap00",sSpecies[kRefs].Data()),
            Form("%s_v3{4,gap00}",sSpecies[kRefs].Data()),
            Form("%s: v_{3}{4,|#Delta#eta| > 0}",sSpecLab[kRefs].Data()),
            mark[kRefs], size[kRefs], col[kRefs]));




    // pt diff

    // for(Int_t s(0); s < kRefs; ++s)
    for(Int_t s(0); s < kK0s; ++s)
    {
        for(Int_t c(0); c < iCent; ++c)
        {
            WriteToFile(fileOut,
                LoadHistoToGraph(fileIn,
                    Form("%s_hFlow2_harm3_gap-10_cent%d",sSpecies[s].Data(),c),
                    Form("%s_v3{2}_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                    Form("%s: v_{3}{2} (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                    mark[s], size[s], col[s]));

            WriteToFile(fileOut,
                LoadHistoToGraph(fileIn,
                    Form("%s_hFlow2_harm3_gap00_cent%d",sSpecies[s].Data(),c),
                    Form("%s_v3{2,gap00}_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                    Form("%s: v_{3}{2,|#Delta#eta| > 0} (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                    mark[s], size[s], col[s]));


            WriteToFile(fileOut,
                LoadHistoToGraph(fileIn,
                    Form("%s_hFlow2_harm3_gap08_cent%d",sSpecies[s].Data(),c),
                    Form("%s_v3{2,gap08}_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                    Form("%s: v_{3}{2,|#Delta#eta| > 0.8} (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                    mark[s], size[s], col[s]));

            WriteToFile(fileOut,
                LoadHistoToGraph(fileIn,
                    Form("%s_hFlow4_harm3_gap-10_cent%d",sSpecies[s].Data(),c),
                    Form("%s_v3{4}_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                    Form("%s: v_{3}{4} (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                    mark[s], size[s], col[s]));

                WriteToFile(fileOut,
                    LoadHistoToGraph(fileIn,
                        Form("%s_hFlow4_harm3_gap00_cent%d",sSpecies[s].Data(),c),
                        Form("%s_v3{4,gap00}_cent%s",sSpecies[s].Data(), sCent[c].Data()),
                        Form("%s: v_{3}{4,|#Delta#eta| > 0} (%s%% V0M)",sSpecLab[s].Data(),sCentLab[c].Data()),
                        mark[s], size[s], col[s]));


            // WriteToFile(fileOut, LoadHistoToGraph(fileIn, "Charged_hFlow2_harm3_gap00_c0", Form("Charged_v3{2,gap00}_c0510", "h^{#pm}: v_{3}{2,|#Delta#eta| > 0} (5-10% V0M)", mark[kCharged], size[kCharged], col[kCharged]));

        }

    }

    fileOut->ls();

    return;

    // Alternative method writing to file
    // {
    //     TGraphErrors* graph = LoadHistoToGraph(fileIn, "Charged_hFlow2_harm3_gap00_c0", "", "h^{#pm}: v_{3}{2,|#Delta#eta| > 0} (0-20% V0M)", kFullSquare, 1.5, kRed);
    //     WriteToFile(fileOut,graph);
    //
    //     TGraphErrors* graph2 = LoadHistoToGraph(fileIn, "Charged_hFlow2_harm3_gap00_c1", "", "h^{#pm}: v_{3}{2,|#Delta#eta| > 0} (20-40% V0M)", kFullSquare, 1.5, kBlack);
    //     WriteToFile(fileOut,graph2);
    //
    //     fileOut->cd();
    //     graph->Write();
    //     graph2->Write();
    // }
}
