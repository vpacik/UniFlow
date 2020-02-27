#include "Fv2_AMPT.C"
#include "/Users/vpacik/Codes/ALICE/Flow/uniFlow/thesis/libStylePID.C"

void Trans_Fv2_AMPT()
{
    TFile* fileOut = TFile::Open("./graphs_PbPb_Fv2_AMPT.root","RECREATE");
    if(!fileOut) { return; }

    // const Int_t row_PT =10;
    // double PT[row_PT]= { 0.15,0.45,0.75,1.05,1.35,1.65,1.95,2.25,2.55,2.85};
    Double_t PT_error[row_PT] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Pions

    // TGraphErrors* gr_pion_Fv2_0005 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_pion_005, PT_error, Fv2total_AMPT_pion_Error_005);
    // TGraphErrors* gr_pion_Fv2_0510 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_pion_510, PT_error, Fv2total_AMPT_pion_Error_510);
    TGraphErrors* gr_pion_Fv2_1020 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_pion_1020, PT_error, Fv2total_AMPT_pion_Error_1020);
    TGraphErrors* gr_pion_Fv2_2030 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_pion_2030, PT_error, Fv2total_AMPT_pion_Error_2030);
    TGraphErrors* gr_pion_Fv2_3040 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_pion_3040, PT_error, Fv2total_AMPT_pion_Error_3040);
    TGraphErrors* gr_pion_Fv2_4050 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_pion_4050, PT_error, Fv2total_AMPT_pion_Error_4050);
    TGraphErrors* gr_pion_Fv2_5060 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_pion_5060, PT_error, Fv2total_AMPT_pion_Error_5060);

    fileOut->cd();
    // gr_pion_Fv2_0005->Write("Pion_Fv2_cent0005");
    // gr_pion_Fv2_0510->Write("Pion_Fv2_cent0510");
    gr_pion_Fv2_1020->Write("Pion_Fv2_cent1020");
    gr_pion_Fv2_2030->Write("Pion_Fv2_cent2030");
    gr_pion_Fv2_3040->Write("Pion_Fv2_cent3040");
    gr_pion_Fv2_4050->Write("Pion_Fv2_cent4050");
    gr_pion_Fv2_5060->Write("Pion_Fv2_cent5060");



    // kaons

    // TGraphErrors* gr_kaon_Fv2_0005 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_kaon_005, PT_error, Fv2total_AMPT_kaon_Error_005);
    // TGraphErrors* gr_kaon_Fv2_0510 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_kaon_510, PT_error, Fv2total_AMPT_kaon_Error_510);
    TGraphErrors* gr_kaon_Fv2_1020 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_kaon_1020, PT_error, Fv2total_AMPT_kaon_Error_1020);
    TGraphErrors* gr_kaon_Fv2_2030 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_kaon_2030, PT_error, Fv2total_AMPT_kaon_Error_2030);
    TGraphErrors* gr_kaon_Fv2_3040 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_kaon_3040, PT_error, Fv2total_AMPT_kaon_Error_3040);
    TGraphErrors* gr_kaon_Fv2_4050 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_kaon_4050, PT_error, Fv2total_AMPT_kaon_Error_4050);
    TGraphErrors* gr_kaon_Fv2_5060 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_kaon_5060, PT_error, Fv2total_AMPT_kaon_Error_5060);

    fileOut->cd();
    // gr_kaon_Fv2_0005->Write("kaon_Fv2_cent0005");
    // gr_kaon_Fv2_0510->Write("kaon_Fv2_cent0510");
    gr_kaon_Fv2_1020->Write("kaon_Fv2_cent1020");
    gr_kaon_Fv2_2030->Write("kaon_Fv2_cent2030");
    gr_kaon_Fv2_3040->Write("kaon_Fv2_cent3040");
    gr_kaon_Fv2_4050->Write("kaon_Fv2_cent4050");
    gr_kaon_Fv2_5060->Write("kaon_Fv2_cent5060");

    // protons

    // TGraphErrors* gr_proton_Fv2_0005 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_proton_005, PT_error, Fv2total_AMPT_proton_Error_005);
    // TGraphErrors* gr_proton_Fv2_0510 = new TGraphErrors(row_PT, PT, Fv2total_AMPT_proton_510, PT_error, Fv2total_AMPT_proton_Error_510);
    TGraphErrors* gr_proton_Fv2_1020 = new TGraphErrors(9, PT_AMPT_proton, Fv2total_AMPT_proton_1020, PT_error, Fv2total_AMPT_proton_Error_1020);
    TGraphErrors* gr_proton_Fv2_2030 = new TGraphErrors(9, PT_AMPT_proton, Fv2total_AMPT_proton_2030, PT_error, Fv2total_AMPT_proton_Error_2030);
    TGraphErrors* gr_proton_Fv2_3040 = new TGraphErrors(9, PT_AMPT_proton, Fv2total_AMPT_proton_3040, PT_error, Fv2total_AMPT_proton_Error_3040);
    TGraphErrors* gr_proton_Fv2_4050 = new TGraphErrors(9, PT_AMPT_proton, Fv2total_AMPT_proton_4050, PT_error, Fv2total_AMPT_proton_Error_4050);
    TGraphErrors* gr_proton_Fv2_5060 = new TGraphErrors(9, PT_AMPT_proton, Fv2total_AMPT_proton_5060, PT_error, Fv2total_AMPT_proton_Error_5060);


    fileOut->cd();
    // gr_proton_Fv2_0005->Write("proton_Fv2_cent0005");
    // gr_proton_Fv2_0510->Write("proton_Fv2_cent0510");
    gr_proton_Fv2_1020->Write("proton_Fv2_cent1020");
    gr_proton_Fv2_2030->Write("proton_Fv2_cent2030");
    gr_proton_Fv2_3040->Write("proton_Fv2_cent3040");
    gr_proton_Fv2_4050->Write("proton_Fv2_cent4050");
    gr_proton_Fv2_5060->Write("proton_Fv2_cent5060");



}
