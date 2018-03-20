/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunProcess()
{
	Double_t dEtaGap = 0.0;
	Int_t iNumSamples = 10;
	TString sEtaGap = "gap00";
	TString sInputPath = "./test/";
	// TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pPb-16qt-nua/";
	// TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pp-16kl";
	// TString sOutputFilePath = sInputPath+"/output-2/"+sEtaGap+"/";
	TString sOutputFilePath = "./test/";
	// TString sOutputFilePath = sInputPath+"/output_forsubt-3/"+sEtaGap+"/";

	// Double_t dMultBinning[] = {0,100};
	Double_t dMultBinning[] = {0,10,20,40,60,100};
	printf("bins %d\n",sizeof(dMultBinning)/sizeof(dMultBinning[0]));

	// Double_t dPtBinningK0s[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,3.,3.5,4.,5.,6.,7.,8.,10.,20.};
	// Double_t dPtBinningK0s[] = {0.2,1.,2.,3.,4.,5.,7.,10.,20.};
	Double_t dPtBinningK0s[] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,7.0,10.0,15.0};
	// Double_t dPtBinningK0s[] = {0.3,0.5,0.75,1.,1.25,1.5,2.,2.5,3.}; // Run1
	Double_t dPtBinningPhi[] = {0.2,1.0,2.0,3.0,4.0,7.0,10.0,15.0};

	// ##### END Parameters setting ######

	// gROOT->AddIncludePath("-I~/NBI/Flow/uniFlow/processUniFlow/");
	gROOT->LoadMacro("~vpacik/NBI/Flow/uniFlow/processUniFlow/ProcessUniFlow.cpp++g");

	FlowTask* taskRefs = new FlowTask(FlowTask::kRefs);
	taskRefs->SetNumSamples(iNumSamples);
	// taskRefs->SetDesamplingUseRMS(kTRUE);
	taskRefs->SetHarmonics(2);
	taskRefs->SetEtaGap(dEtaGap);
	taskRefs->SetMergePosNeg();


	FlowTask* taskCharged = new FlowTask(FlowTask::kCharged);
	taskCharged->SetNumSamples(iNumSamples);
	taskCharged->SetHarmonics(2);
	taskCharged->SetEtaGap(dEtaGap);
	taskCharged->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskCharged->SetMergePosNeg();


	FlowTask* taskPion = new FlowTask(FlowTask::kPion);
	taskPion->SetNumSamples(iNumSamples);
	taskPion->SetHarmonics(2);
	taskPion->SetEtaGap(dEtaGap);
	taskPion->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskPion->SetMergePosNeg();

	FlowTask* taskKch = new FlowTask(FlowTask::kKaon);
	taskKch->SetNumSamples(iNumSamples);
	taskKch->SetHarmonics(2);
	taskKch->SetEtaGap(dEtaGap);
	taskKch->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskKch->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

	FlowTask* taskProton = new FlowTask(FlowTask::kProton);
	taskProton->SetNumSamples(iNumSamples);
	taskProton->SetHarmonics(2);
	taskProton->SetEtaGap(dEtaGap);
	taskProton->SetPtBins(dPtBinningK0s,sizeof(dPtBinningK0s)/sizeof(dPtBinningK0s[0]));
	taskProton->SetMergePosNeg();
	// taskK0s->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");

 	FlowTask* taskK0s = new FlowTask(FlowTask::kK0s);
	taskK0s->SetAlexFitting(kTRUE);
	taskK0s->SetNumSamples(iNumSamples);
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
	taskLambda->SetNumSamples(iNumSamples);
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
	taskPhi->SetNumSamples(iNumSamples);
	taskPhi->SetHarmonics(2);
	taskPhi->SetEtaGap(dEtaGap);
	taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	taskPhi->SetPtBins(dPtBinningPhi,sizeof(dPtBinningPhi)/sizeof(dPtBinningPhi[0]));
	taskPhi->SetMergePosNeg();
	// taskPhi->SetFittingRange(0.45,0.55);
	// taskPhi->SetFittingRejectNumSigmas(3);
	// taskPhi->SetAlternativeProfileName("fp3V0sCorrK0s_<2>_harm2_gap08_Neg");



	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName("AnalysisResults.root");
	process->SetTaskName("UniFlow");
	process->SetDebug();
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("Processed.root");
	process->SetMultiplicityBins(dMultBinning,sizeof(dMultBinning)/sizeof(dMultBinning[0]));
	process->SetSaveMult(kTRUE);
	process->SetFitCumulants(kTRUE);
	process->AddTask(taskRefs);
	process->AddTask(taskCharged);
	process->AddTask(taskPion);
	process->AddTask(taskKch);
	process->AddTask(taskProton);
	// process->AddTask(taskK0s);
	// process->AddTask(taskLambda);
	// process->AddTask(taskPhi);
	process->Run();

	return;
}
