#include "v2{2}(pT)_v2{4}(pT)_AMPT.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStylePID.C"

void Trans_cums_AMPT()
{
    TFile* fileOut = TFile::Open("./graphs_PbPb_v22_v24_AMPT.root","RECREATE");
    if(!fileOut) { return; }

    // const Int_t row_PT =10;
    // double PT[row_PT]= { 0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85};
    Double_t PT_error[row_PT] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Pions

    TGraphErrors* gr_pion_v22_0005 = new TGraphErrors(row_PT, PT, v22pt_AMPT_pion_005, PT_error, v22pt_AMPT_Error_pion_005);
    TGraphErrors* gr_pion_v22_0510 = new TGraphErrors(row_PT, PT, v22pt_AMPT_pion_510, PT_error, v22pt_AMPT_Error_pion_510);
    TGraphErrors* gr_pion_v22_1020 = new TGraphErrors(row_PT, PT, v22pt_AMPT_pion_1020, PT_error, v22pt_AMPT_Error_pion_1020);
    TGraphErrors* gr_pion_v22_2030 = new TGraphErrors(row_PT, PT, v22pt_AMPT_pion_2030, PT_error, v22pt_AMPT_Error_pion_2030);
    TGraphErrors* gr_pion_v22_3040 = new TGraphErrors(row_PT, PT, v22pt_AMPT_pion_3040, PT_error, v22pt_AMPT_Error_pion_3040);
    TGraphErrors* gr_pion_v22_4050 = new TGraphErrors(row_PT, PT, v22pt_AMPT_pion_4050, PT_error, v22pt_AMPT_Error_pion_4050);
    TGraphErrors* gr_pion_v22_5060 = new TGraphErrors(row_PT, PT, v22pt_AMPT_pion_5060, PT_error, v22pt_AMPT_Error_pion_5060);

    // TGraphErrors* gr_pion_v24_0005 = new TGraphErrors(row_PT, PT, v24pt_AMPT_pion_005, PT_error, v24pt_AMPT_Error_pion_005);
    TGraphErrors* gr_pion_v24_0510 = new TGraphErrors(row_PT, PT, v24pt_AMPT_pion_510, PT_error, v24pt_AMPT_Error_pion_510);
    TGraphErrors* gr_pion_v24_1020 = new TGraphErrors(row_PT, PT, v24pt_AMPT_pion_1020, PT_error, v24pt_AMPT_Error_pion_1020);
    TGraphErrors* gr_pion_v24_2030 = new TGraphErrors(row_PT, PT, v24pt_AMPT_pion_2030, PT_error, v24pt_AMPT_Error_pion_2030);
    TGraphErrors* gr_pion_v24_3040 = new TGraphErrors(row_PT, PT, v24pt_AMPT_pion_3040, PT_error, v24pt_AMPT_Error_pion_3040);
    TGraphErrors* gr_pion_v24_4050 = new TGraphErrors(row_PT, PT, v24pt_AMPT_pion_4050, PT_error, v24pt_AMPT_Error_pion_4050);
    TGraphErrors* gr_pion_v24_5060 = new TGraphErrors(row_PT, PT, v24pt_AMPT_pion_5060, PT_error, v24pt_AMPT_Error_pion_5060);

    fileOut->cd();
    gr_pion_v22_0005->Write("Pion_v22_cent0005");
    gr_pion_v22_0510->Write("Pion_v22_cent0510");
    gr_pion_v22_1020->Write("Pion_v22_cent1020");
    gr_pion_v22_2030->Write("Pion_v22_cent2030");
    gr_pion_v22_3040->Write("Pion_v22_cent3040");
    gr_pion_v22_4050->Write("Pion_v22_cent4050");
    gr_pion_v22_5060->Write("Pion_v22_cent5060");

    // gr_pion_v24_0005->Write("pion_v24_cent0005");
    gr_pion_v24_0510->Write("pion_v24_cent0510");
    gr_pion_v24_1020->Write("pion_v24_cent1020");
    gr_pion_v24_2030->Write("pion_v24_cent2030");
    gr_pion_v24_3040->Write("pion_v24_cent3040");
    gr_pion_v24_4050->Write("pion_v24_cent4050");
    gr_pion_v24_5060->Write("pion_v24_cent5060");


    // kaons

    TGraphErrors* gr_kaon_v22_0005 = new TGraphErrors(row_PT, PT, v22pt_AMPT_kaon_005, PT_error, v22pt_AMPT_Error_kaon_005);
    TGraphErrors* gr_kaon_v22_0510 = new TGraphErrors(row_PT, PT, v22pt_AMPT_kaon_510, PT_error, v22pt_AMPT_Error_kaon_510);
    TGraphErrors* gr_kaon_v22_1020 = new TGraphErrors(row_PT, PT, v22pt_AMPT_kaon_1020, PT_error, v22pt_AMPT_Error_kaon_1020);
    TGraphErrors* gr_kaon_v22_2030 = new TGraphErrors(row_PT, PT, v22pt_AMPT_kaon_2030, PT_error, v22pt_AMPT_Error_kaon_2030);
    TGraphErrors* gr_kaon_v22_3040 = new TGraphErrors(row_PT, PT, v22pt_AMPT_kaon_3040, PT_error, v22pt_AMPT_Error_kaon_3040);
    TGraphErrors* gr_kaon_v22_4050 = new TGraphErrors(row_PT, PT, v22pt_AMPT_kaon_4050, PT_error, v22pt_AMPT_Error_kaon_4050);
    TGraphErrors* gr_kaon_v22_5060 = new TGraphErrors(row_PT, PT, v22pt_AMPT_kaon_5060, PT_error, v22pt_AMPT_Error_kaon_5060);

    // TGraphErrors* gr_kaon_v24_0005 = new TGraphErrors(row_PT, PT, v24pt_AMPT_kaon_005, PT_error, v24pt_AMPT_Error_kaon_005);
    TGraphErrors* gr_kaon_v24_0510 = new TGraphErrors(row_PT, PT, v24pt_AMPT_kaon_510, PT_error, v24pt_AMPT_Error_kaon_510);
    TGraphErrors* gr_kaon_v24_1020 = new TGraphErrors(row_PT, PT, v24pt_AMPT_kaon_1020, PT_error, v24pt_AMPT_Error_kaon_1020);
    TGraphErrors* gr_kaon_v24_2030 = new TGraphErrors(row_PT, PT, v24pt_AMPT_kaon_2030, PT_error, v24pt_AMPT_Error_kaon_2030);
    TGraphErrors* gr_kaon_v24_3040 = new TGraphErrors(row_PT, PT, v24pt_AMPT_kaon_3040, PT_error, v24pt_AMPT_Error_kaon_3040);
    TGraphErrors* gr_kaon_v24_4050 = new TGraphErrors(row_PT, PT, v24pt_AMPT_kaon_4050, PT_error, v24pt_AMPT_Error_kaon_4050);
    TGraphErrors* gr_kaon_v24_5060 = new TGraphErrors(row_PT, PT, v24pt_AMPT_kaon_5060, PT_error, v24pt_AMPT_Error_kaon_5060);

    fileOut->cd();
    gr_kaon_v22_0005->Write("kaon_v22_cent0005");
    gr_kaon_v22_0510->Write("kaon_v22_cent0510");
    gr_kaon_v22_1020->Write("kaon_v22_cent1020");
    gr_kaon_v22_2030->Write("kaon_v22_cent2030");
    gr_kaon_v22_3040->Write("kaon_v22_cent3040");
    gr_kaon_v22_4050->Write("kaon_v22_cent4050");
    gr_kaon_v22_5060->Write("kaon_v22_cent5060");

    // gr_kaon_v24_0005->Write("kaon_v24_cent0005");
    gr_kaon_v24_0510->Write("kaon_v24_cent0510");
    gr_kaon_v24_1020->Write("kaon_v24_cent1020");
    gr_kaon_v24_2030->Write("kaon_v24_cent2030");
    gr_kaon_v24_3040->Write("kaon_v24_cent3040");
    gr_kaon_v24_4050->Write("kaon_v24_cent4050");
    gr_kaon_v24_5060->Write("kaon_v24_cent5060");


    // protons

    TGraphErrors* gr_proton_v22_0005 = new TGraphErrors(9, PT_AMPT_proton, v22pt_AMPT_proton_005, PT_error, v22pt_AMPT_Error_proton_005);
    TGraphErrors* gr_proton_v22_0510 = new TGraphErrors(9, PT_AMPT_proton, v22pt_AMPT_proton_510, PT_error, v22pt_AMPT_Error_proton_510);
    TGraphErrors* gr_proton_v22_1020 = new TGraphErrors(9, PT_AMPT_proton, v22pt_AMPT_proton_1020, PT_error, v22pt_AMPT_Error_proton_1020);
    TGraphErrors* gr_proton_v22_2030 = new TGraphErrors(9, PT_AMPT_proton, v22pt_AMPT_proton_2030, PT_error, v22pt_AMPT_Error_proton_2030);
    TGraphErrors* gr_proton_v22_3040 = new TGraphErrors(9, PT_AMPT_proton, v22pt_AMPT_proton_3040, PT_error, v22pt_AMPT_Error_proton_3040);
    TGraphErrors* gr_proton_v22_4050 = new TGraphErrors(9, PT_AMPT_proton, v22pt_AMPT_proton_4050, PT_error, v22pt_AMPT_Error_proton_4050);
    TGraphErrors* gr_proton_v22_5060 = new TGraphErrors(9, PT_AMPT_proton, v22pt_AMPT_proton_5060, PT_error, v22pt_AMPT_Error_proton_5060);

    // TGraphErrors* gr_proton_v24_0005 = new TGraphErrors(row_PT, PT, v24pt_AMPT_proton_005, PT_error, v24pt_AMPT_Error_proton_005);
    TGraphErrors* gr_proton_v24_0510 = new TGraphErrors(9, PT_AMPT_proton, v24pt_AMPT_proton_510, PT_error, v24pt_AMPT_Error_proton_510);
    TGraphErrors* gr_proton_v24_1020 = new TGraphErrors(9, PT_AMPT_proton, v24pt_AMPT_proton_1020, PT_error, v24pt_AMPT_Error_proton_1020);
    TGraphErrors* gr_proton_v24_2030 = new TGraphErrors(9, PT_AMPT_proton, v24pt_AMPT_proton_2030, PT_error, v24pt_AMPT_Error_proton_2030);
    TGraphErrors* gr_proton_v24_3040 = new TGraphErrors(9, PT_AMPT_proton, v24pt_AMPT_proton_3040, PT_error, v24pt_AMPT_Error_proton_3040);
    TGraphErrors* gr_proton_v24_4050 = new TGraphErrors(9, PT_AMPT_proton, v24pt_AMPT_proton_4050, PT_error, v24pt_AMPT_Error_proton_4050);
    TGraphErrors* gr_proton_v24_5060 = new TGraphErrors(9, PT_AMPT_proton, v24pt_AMPT_proton_5060, PT_error, v24pt_AMPT_Error_proton_5060);

    fileOut->cd();
    gr_proton_v22_0005->Write("proton_v22_cent0005");
    gr_proton_v22_0510->Write("proton_v22_cent0510");
    gr_proton_v22_1020->Write("proton_v22_cent1020");
    gr_proton_v22_2030->Write("proton_v22_cent2030");
    gr_proton_v22_3040->Write("proton_v22_cent3040");
    gr_proton_v22_4050->Write("proton_v22_cent4050");
    gr_proton_v22_5060->Write("proton_v22_cent5060");

    // gr_proton_v24_0005->Write("proton_v24_cent0005");
    gr_proton_v24_0510->Write("proton_v24_cent0510");
    gr_proton_v24_1020->Write("proton_v24_cent1020");
    gr_proton_v24_2030->Write("proton_v24_cent2030");
    gr_proton_v24_3040->Write("proton_v24_cent3040");
    gr_proton_v24_4050->Write("proton_v24_cent4050");
    gr_proton_v24_5060->Write("proton_v24_cent5060");


}
