/* RunF
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

void RunFlowSubMultInt(const char* sOutputFilePath = "")
{
	// ##### Parameters setting ######
	Double_t dEtaGap = 0.0;
	// Double_t dMultBinning[] = {20,40};
	Double_t dMultBinning[] = {0,100};
	// Double_t dMultBinning[] = {0,20,40,60,100};
	// Double_t dMultBinning[] = {0.,5.,10.,20.,40.,60.,100.};
	// Double_t dMultBinning[] = {0,1,2,3,4,5,6}; // fixed mult. binning

	// Double_t dPtBinningK0s[] = {3.,5.};
	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	Double_t dPtBinningK0s[] = {0.2,1.,2.,3.,4.,5.,7.,10.,20.};
	// Double_t dPtBinningK0s[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // Run1
	const char* sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pp-run2-gap0";

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
	process->SetOutputFilePath(Form("%s/output_multint",sInputPath));
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->AddTask(taskRefs);
	process->AddTask(taskCharged);
	process->AddTask(taskPion);
	process->AddTask(taskKch);
	process->AddTask(taskProton);
	process->AddTask(taskK0s);
	process->AddTask(taskLambda);
	process->Run();

	ProcessUniFlow* process_noweight = new ProcessUniFlow();
	process_noweight->SetInputFilePath(sInputPath);
	process_noweight->SetInputFileName("AnalysisResults.root");
	process_noweight->SetTaskName("UniFlow");
	process_noweight->SetGlobalProfNameLabel("multScaled_");
	// process_noweight->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_noweight->SetOutputFilePath(Form("%s/output_multint_noweight",sInputPath));
	process_noweight->SetOutputFileName("Processed.root");
	process_noweight->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_noweight->AddTask(taskRefs);
	process_noweight->AddTask(taskCharged);
	process_noweight->AddTask(taskPion);
	process_noweight->AddTask(taskKch);
	process_noweight->AddTask(taskProton);
	process_noweight->AddTask(taskK0s);
	process_noweight->AddTask(taskLambda);
	process_noweight->Run();

	ProcessUniFlow* process_sub = new ProcessUniFlow();
	process_sub->SetInputFilePath(sInputPath);
	process_sub->SetInputFileName("AnalysisResults.root");
	process_sub->SetTaskName("UniFlow_sub");
	// process_sub->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub->SetOutputFilePath(Form("%s/output_multint_sub",sInputPath));
	process_sub->SetOutputFileName("Processed.root");
	process_sub->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub->AddTask(taskRefs);
	process_sub->AddTask(taskCharged);
	process_sub->AddTask(taskPion);
	process_sub->AddTask(taskKch);
	process_sub->AddTask(taskProton);
	process_sub->AddTask(taskK0s);
	process_sub->AddTask(taskLambda);
	process_sub->Run();

	ProcessUniFlow* process_sub_norm = new ProcessUniFlow();
	process_sub_norm->SetInputFilePath(sInputPath);
	process_sub_norm->SetInputFileName("AnalysisResults.root");
	process_sub_norm->SetTaskName("UniFlow_sub");
	process_sub_norm->SetGlobalProfNameLabel("multScaled_");
	// process_sub_norm->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub_norm->SetOutputFilePath(Form("%s/output_multint_sub_norm",sInputPath));
	process_sub_norm->SetOutputFileName("Processed.root");
	process_sub_norm->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub_norm->AddTask(taskRefs);
	process_sub_norm->AddTask(taskCharged);
	process_sub_norm->AddTask(taskPion);
	process_sub_norm->AddTask(taskKch);
	process_sub_norm->AddTask(taskProton);
	process_sub_norm->AddTask(taskK0s);
	process_sub_norm->AddTask(taskLambda);
	process_sub_norm->Run();

	ProcessUniFlow* process_sub_norm_weighted = new ProcessUniFlow();
	process_sub_norm_weighted->SetInputFilePath(sInputPath);
	process_sub_norm_weighted->SetInputFileName("AnalysisResults.root");
	process_sub_norm_weighted->SetTaskName("UniFlow_sub");
	process_sub_norm_weighted->SetGlobalProfNameLabel("multScaled_weighted_");
	// process_sub_norm_weighted->SetOutputFilePath("/Users/vpacik/NBI/Flow/uniFlow/results/KchK0s/K0s/plots");
	process_sub_norm_weighted->SetOutputFilePath(Form("%s/output_multint_sub_norm_weighted",sInputPath));
	process_sub_norm_weighted->SetOutputFileName("Processed.root");
	process_sub_norm_weighted->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process_sub_norm_weighted->AddTask(taskRefs);
	process_sub_norm_weighted->AddTask(taskCharged);
	process_sub_norm_weighted->AddTask(taskPion);
	process_sub_norm_weighted->AddTask(taskKch);
	process_sub_norm_weighted->AddTask(taskProton);
	process_sub_norm_weighted->AddTask(taskK0s);
	process_sub_norm_weighted->AddTask(taskLambda);
	process_sub_norm_weighted->Run();

	return;
}
