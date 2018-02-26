/* RunF
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunFlowSub(const char* sOutputFilePath = "")
{
	// ##### Parameters setting ######
	Double_t dEtaGap = 0.0;
	// Double_t dMultBinning[] = {20,40};
	// Double_t dMultBinning[] = {0,20,40,60,100};
	// Double_t dMultBinning[] = {0.,5.,10.,20.,40.,60.,100.};

	// N_RFPs multiplicity
	// Double_t dMultBinning[] = {0,150};
	// Double_t dMultBinning[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150};
	// Double_t dMultBinning[] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150};
	Double_t dMultBinning[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
	printf("bins %d\n",sizeof(dMultBinning)/sizeof(dMultBinning[0]));


	// Double_t dMultBinning[] = {1,2,3,4,5,6,7,8,9,10};

	// Double_t dMultBinning[] = {0,1,2,3,4,5,6}; // fixed mult. binning

	// Double_t dPtBinningK0s[] = {3.,5.};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningK0s[] = {0.2,1.,2.,3.,4.,5.,7.,10.,20.};
	Double_t dPtBinningK0s[] = {0.2,0.5,1.,1.5,2.,2.5,3.,3.5,4.,5.,7.,10.,20.};
	// Double_t dPtBinningK0s[] = {0.0,2.0,4.0,6.0,8.0,10.0}; // large bins
	// Double_t dPtBinningK0s[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // Run1
	const char* sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-NchRFP-gap00";

	// ##### END Parameters setting ######

	// gROOT->AddIncludePath("-I~/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->LoadMacro("~vpacik/NBI/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");


 	FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
	taskK0s->SetAlexFitting(kTRUE);
	taskK0s->SetNumSamples(1);
	taskK0s->SetHarmonics(2);
	taskK0s->SetEtaGap(dEtaGap);
	taskK0s->SetInvMassRebin(2);
	taskK0s->SetFlowMassRebin(2);
	taskK0s->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskK0s->SetMergePosNeg();
	// taskK0s->SetFittingRange(0.45,0.55);
	// taskK0s->SetFittingRejectNumSigmas(3);
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskLambda = new FlowTask(FlowTask::kLambda	);
	taskLambda->SetAlexFitting(kTRUE);
	taskLambda->SetNumSamples(1);
	taskLambda->SetHarmonics(2);
	taskLambda->SetEtaGap(dEtaGap);
	taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	taskLambda->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskLambda->SetMergePosNeg();
	// taskLambda->SetFittingRange(0.45,0.55);
	// taskLambda->SetFittingRejectNumSigmas(3);
	// taskLambda->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskPhi = new FlowTask(FlowTask::kPhi);
	taskPhi->SetAlexFitting(kTRUE);
	taskPhi->SetNumSamples(1);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(dEtaGap);
	taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskPhi->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetFittingRejectNumSigmas(3);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");


	FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
	taskRefs->SetNumSamples(1);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(dEtaGap);
	// taskCharged->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskRefs->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");


	FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
	taskCharged->SetNumSamples(1);
	taskCharged->SetHarmonics(2);
	taskCharged->SetEtaGap(dEtaGap);
	taskCharged->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskCharged->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");


	FlowTask* taskPion = new FlowTask(FlowTask::kPion);
	taskPion->SetNumSamples(1);
	taskPion->SetHarmonics(2);
	taskPion->SetEtaGap(dEtaGap);
	taskPion->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskPion->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskKch = new FlowTask(FlowTask::kKaon);
	taskKch->SetNumSamples(1);
	taskKch->SetHarmonics(2);
	taskKch->SetEtaGap(dEtaGap);
	taskKch->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskKch->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskProton = new FlowTask(FlowTask::kProton);
	taskProton->SetNumSamples(1);
	taskProton->SetHarmonics(2);
	taskProton->SetEtaGap(dEtaGap);
	taskProton->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskProton->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");


	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath);
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	// process->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process->SetOutputFilePath(Form("%s/output_vn_5/GF_eventweighted",sInputPath));
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->AddTask(taskRefs);
	process->AddTask(taskCharged);
	// process->AddTask(taskPion);
	// process->AddTask(taskKch);
	// process->AddTask(taskProton);
	// process->AddTask(taskK0s);
	// process->AddTask(taskLambda);
	// process->AddTask(taskPhi);
	process->Run();

	ProcessUniFlow* process_noweight = new ProcessUniFlow();
	process_noweight->SetInputFilePath(sInputPath);
	process_noweight->SetInputFileName("AnalysisResults.root");
	process_noweight->SetTaskName("UniFlow");
	process_noweight->SetGlobalProfNameLabel("multScaled_");
	// process_noweight->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_noweight->SetOutputFilePath(Form("%s/output_vn_5/GF_noneventweighted",sInputPath));
	process_noweight->SetOutputFileName("Processed.root");
	process_noweight->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_noweight->AddTask(taskRefs);
	process_noweight->AddTask(taskCharged);
	// process_noweight->AddTask(taskPion);
	// process_noweight->AddTask(taskKch);
	// process_noweight->AddTask(taskProton);
	// process_noweight->AddTask(taskK0s);
	// process_noweight->AddTask(taskLambda);
	// process_noweight->AddTask(taskPhi);
	process_noweight->Run();

	ProcessUniFlow* process_sub = new ProcessUniFlow();
	process_sub->SetInputFilePath(sInputPath);
	process_sub->SetInputFileName("AnalysisResults.root");
	process_sub->SetTaskName("UniFlow_sub");
	// process_sub->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub->SetOutputFilePath(Form("%s/output_vn_5/SP_nonscaled_noneventweighted",sInputPath));
	process_sub->SetOutputFileName("Processed.root");
	process_sub->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub->SetSaveMult("Mult.root");
	process_sub->AddTask(taskRefs);
	process_sub->AddTask(taskCharged);
	// process_sub->AddTask(taskPion);
	// process_sub->AddTask(taskKch);
	// process_sub->AddTask(taskProton);
	// process_sub->AddTask(taskK0s);
	// process_sub->AddTask(taskLambda);
	// process_sub->AddTask(taskPhi);
	process_sub->Run();

	ProcessUniFlow* process_sub_norm = new ProcessUniFlow();
	process_sub_norm->SetInputFilePath(sInputPath);
	process_sub_norm->SetInputFileName("AnalysisResults.root");
	process_sub_norm->SetTaskName("UniFlow_sub");
	process_sub_norm->SetGlobalProfNameLabel("multScaled_");
	// process_sub_norm->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub_norm->SetOutputFilePath(Form("%s/output_vn_5/SP_scaled_noneventweighted",sInputPath));
	process_sub_norm->SetOutputFileName("Processed.root");
	process_sub_norm->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub_norm->AddTask(taskRefs);
	process_sub_norm->AddTask(taskCharged);
	// process_sub_norm->AddTask(taskPion);
	// process_sub_norm->AddTask(taskKch);
	// process_sub_norm->AddTask(taskProton);
	// process_sub_norm->AddTask(taskK0s);
	// process_sub_norm->AddTask(taskLambda);
	// process_sub_norm->AddTask(taskPhi);
	process_sub_norm->Run();

	ProcessUniFlow* process_sub_norm_weighted = new ProcessUniFlow();
	process_sub_norm_weighted->SetInputFilePath(sInputPath);
	process_sub_norm_weighted->SetInputFileName("AnalysisResults.root");
	process_sub_norm_weighted->SetTaskName("UniFlow_sub");
	process_sub_norm_weighted->SetGlobalProfNameLabel("multScaled_weighted_");
	// process_sub_norm_weighted->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub_norm_weighted->SetOutputFilePath(Form("%s/output_vn_5/SP_scaled_eventweighted",sInputPath));
	process_sub_norm_weighted->SetOutputFileName("Processed.root");
	process_sub_norm_weighted->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub_norm_weighted->AddTask(taskRefs);
	process_sub_norm_weighted->AddTask(taskCharged);
	// process_sub_norm_weighted->AddTask(taskPion);
	// process_sub_norm_weighted->AddTask(taskKch);
	// process_sub_norm_weighted->AddTask(taskProton);
	// process_sub_norm_weighted->AddTask(taskK0s);
	// process_sub_norm_weighted->AddTask(taskLambda);
	// process_sub_norm_weighted->AddTask(taskPhi);
	process_sub_norm_weighted->Run();

	return;
}
