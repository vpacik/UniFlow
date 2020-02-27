// NEED to run .X Load.C beforehand

/* RunProcess.C
 *
 * Steer macro for procesing flow results of AliAnalysisTaskUniFlow task.
 * See ProcessUniFlow.cpp for class implementation.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

void RunCumsV2()
{
	Int_t iNumSamples = 1;

	// // Pb-Pb
	TString sInputPath = "/mnt/Codes/Flow/uniFlow/thesis/results/PbPb/";
	TString sInputFile = "train_CF_PbPb_6815.root";

	Int_t iHarmonics = 2; TString sOutputFilePath = Form("/mnt/Codes/Flow/uniFlow/thesis/results/PbPb/cums/v2/process/");

	Bool_t bDoK0s = kFALSE;
	Bool_t bDoLambda = kFALSE;
	Bool_t bDoPhi = kFALSE;


	bDoK0s = kTRUE;
	bDoLambda = kTRUE;
	bDoPhi = kTRUE;


	// // ====== Starting points
	// v2{4}
	// std::vector<Double_t> dMultBinning = {50,150};
	std::vector<Double_t> dMultBinning = {0,10,20,30,40,50,60};
	// std::vector<Double_t> dMultBinningPhi = {5,10,20,30,40,50,60};

	std::vector<Double_t> dPtBinning[kUnknown] = {};
	dPtBinning[kCharged] =  {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
	dPtBinning[kPion] =  {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
	dPtBinning[kKaon] =  {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
	dPtBinning[kProton] =  {0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0};
	dPtBinning[kK0s] = {0.5,0.7,0.9,1.1,1.3,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0,6.0,10.0};
	dPtBinning[kLambda] = {0.6,1.0,1.25,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0,6.0,10.0};
	dPtBinning[kPhi] = {1.0,1.5,2.0,3.0,4.0,6.0};

	// ##### END Parameters setting ######

	ProcessUniFlow* process = new ProcessUniFlow();
	process->SetInputFilePath(sInputPath.Data());
	process->SetInputFileName(SetInputFile.Data());
	process->SetTaskName("UniFlowgap00");
	process->SetOutputFilePath(sOutputFilePath.Data());
	process->SetOutputFileName("ProcessedAll.root");
	process->SetMultiplicityBins(dMultBinning);
	process->SetSaveMult(0);
	process->SetSaveInterSteps(1);
	process->SetFitCumulants(kFALSE);
	process->SetDebug(0);

	process->AddTask(new FlowTask(kRefs, iHarmonics, 4, -1.0, iNumSamples));
	process->AddTask(new FlowTask(kRefs, iHarmonics, 4, 0.0, iNumSamples));
	process->AddTask(new FlowTask(kRefs, iHarmonics, 2, 0.8, iNumSamples));

	process->AddTask(new FlowTask(kCharged, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kCharged]));
	process->AddTask(new FlowTask(kCharged, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kCharged]));
	process->AddTask(new FlowTask(kCharged, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kCharged]));
	//
	process->AddTask(new FlowTask(kPion, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kPion]));
	process->AddTask(new FlowTask(kPion, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kPion]));
	process->AddTask(new FlowTask(kPion, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kPion]));

	process->AddTask(new FlowTask(kKaon, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kKaon]));
	process->AddTask(new FlowTask(kKaon, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kKaon]));
	process->AddTask(new FlowTask(kKaon, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kKaon]));

	process->AddTask(new FlowTask(kProton, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kProton]));
	process->AddTask(new FlowTask(kProton, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kProton]));
	process->AddTask(new FlowTask(kProton, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kProton]));
	//
 	FlowTask* taskK0s = new FlowTask(kK0s, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kK0s]);
	if(bDoK0s) process->AddTask(taskK0s);

	FlowTask* taskK0s3 = new FlowTask(kK0s, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kK0s]);
	if(bDoK0s) process->AddTask(taskK0s3);

	FlowTask* taskK0s2 = new FlowTask(kK0s, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kK0s]);
	if(bDoK0s) process->AddTask(taskK0s2);

	FlowTask* taskLambda = new FlowTask(kLambda, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kLambda]);
	// taskLambda->SetInvMassRebin(2);
	taskLambda->SetFlowMassRebin(2);
	if(bDoLambda) process->AddTask(taskLambda);

	FlowTask* taskLambda3 = new FlowTask(kLambda, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kLambda]);
	// taskLambda3->SetInvMassRebin(2);
	taskLambda3->SetFlowMassRebin(2);
	if(bDoLambda) process->AddTask(taskLambda3);

 	FlowTask* taskLambda2 = new FlowTask(kLambda, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kLambda]);
	// taskLambda2->SetInvMassRebin(2);
	taskLambda2->SetFlowMassRebin(2);
	if(bDoLambda) process->AddTask(taskLambda2);

	FlowTask* taskPhi = new FlowTask(kPhi, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kPhi]);
	// taskPhi->SetInvMassRebin(2);
	taskPhi->SetFlowMassRebin(2);
	if(bDoPhi) process->AddTask(taskPhi);

	FlowTask* taskPhi2 = new FlowTask(kPhi, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kPhi]);
	// taskPhi2->SetInvMassRebin(2);
	taskPhi2->SetFlowMassRebin(2);
	if(bDoPhi) process->AddTask(taskPhi2);

	FlowTask* taskPhi3 = new FlowTask(kPhi, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kPhi]);
	// taskPhi3->SetInvMassRebin(2);
	taskPhi3->SetFlowMassRebin(2);
	if(bDoPhi) process->AddTask(taskPhi3);


	// ProcessUniFlow* processPhi = new ProcessUniFlow();
	// processPhi->SetInputFilePath(sInputPath.Data());
	// processPhi->SetInputFileName("AnalysisResults.root");
	// processPhi->SetTaskName("UniFlowgap00");
	// processPhi->SetOutputFilePath(sOutputFilePath.Data());
	// processPhi->SetOutputFileName("ProcessedPhi.root");
	// processPhi->SetMultiplicityBins(dMultBinningPhi);
	// processPhi->SetSaveMult(0);
	// processPhi->SetFitCumulants(kFALSE);
	// processPhi->SetDebug(1);
	// processPhi->AddTask(new FlowTask(kRefs, iHarmonics, 4, -1.0, iNumSamples));
	// processPhi->AddTask(new FlowTask(kRefs, iHarmonics, 4, 0.0, iNumSamples));
	// processPhi->AddTask(new FlowTask(kRefs, iHarmonics, 2, 0.8, iNumSamples));
	//
	//
	//
	// FlowTask* taskPhi = new FlowTask(kPhi, iHarmonics, 4, -1.0, iNumSamples, dPtBinning[kPhi]);
	// // taskPhi->SetInvMassRebin(2);
	// // taskPhi->SetFlowMassRebin(2);
	// if(bDoPhi) processPhi->AddTask(taskPhi);
	//
	// FlowTask* taskPhi2 = new FlowTask(kPhi, iHarmonics, 4, 0.0, iNumSamples, dPtBinning[kPhi]);
	// // taskPhi2->SetInvMassRebin(2);
	// // taskPhi2->SetFlowMassRebin(2);
	// if(bDoPhi) processPhi->AddTask(taskPhi2);
	//
	// FlowTask* taskPhi3 = new FlowTask(kPhi, iHarmonics, 2, 0.8, iNumSamples, dPtBinning[kPhi]);
	// // taskPhi3->SetInvMassRebin(2);
	// // taskPhi3->SetFlowMassRebin(2);
	// if(bDoPhi) processPhi->AddTask(taskPhi3);

	process->Run();
	// processPhi->Run();

	return;
}
