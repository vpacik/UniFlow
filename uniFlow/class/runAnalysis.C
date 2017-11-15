void runAnalysis()
{
    Bool_t local = 1; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 0; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)

    TString sGridMode = "full";
    //TString sGridMode = "terminate";

    Bool_t bMergeViaJDL = kTRUE;
    //Bool_t bMergeViaJDL = kFALSE;

    TString sWorkDir = "uniflow_testD";
    TString sOutDir = "output";
    TString sPeriod = "LHC16q";

    // run switcher
    // Run2 8.16 TeV
    // RunList_LHC16r_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [12 runs ~16,6M ]
    //Int_t runNumber[] = {266318, 266317, 266316,   266208, 266197, 266196, 266187, 265754, 265744, 265607, 265596, 265594};
    // RunList_LHC16s_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [16 runs ~10,5M ]
    // Int_t runNumber[] = {267110, 267081, 267077, 267072, 267070, 267030, 266998, 266997, 266994, 266993, 266944, 266886, 266885, 266883, 266882, 266437};

    // Run2 5.02 TeV
    // RunList_LHC16t_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [4 runs]
    // Int_t runNumber[] = {267166, 267165, 267164, 267163};
    // RunList_LHC16q_pass1_CentralBarrelTracking_hadronPID_20170318_v1.txt [31 runs]
    // Int_t runNumber[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387};
    // Int_t runNumber[] = {265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};


    //test
    Int_t runNumber[] = {265385, 265384, 265383};


    // since we will compile a class, tell root where to look for headers
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");

    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    //gSystem->AddIncludePath("-I$ALICE_PHYSICS/../src/OADB/COMMON/MULTIPLICITY");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSISaliceBase");
    gSystem->Load("libCORRFW");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskUniFlow");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    // Physics selection as suggested by HMTF
    gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE,kTRUE);

    //Add this here: run before your task, but after definition of manager and input handler
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask* taskMultSelection = AddTaskMultSelection(kFALSE); // user mode:
    taskMultSelection->SetSelectedTriggerClass(AliVEvent::kINT7);

    // PID response needed for PID
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse* taskPIDResponse = AddTaskPIDResponse(kFALSE); // not MC

    /*
    // PID response QA by ALICE
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa* taskPIDqa = AddTaskPIDqa("PID_QA.root");
    */

    gROOT->LoadMacro("AliAnalysisTaskUniFlow.cxx++g"); // compile the class (locally)
    gROOT->LoadMacro("AddTaskUniFlow.C"); // load the addtask macro

    AliAnalysisTaskUniFlow* task1 = AddTaskUniFlow("UniFlow");
    // Analysis
    task1->SetRunMode(AliAnalysisTaskUniFlow::kFull);
    task1->SetNumEventsAnalyse(50);
    task1->SetAnalysisType(AliAnalysisTaskUniFlow::kAOD);
    // task1->SetSampling(kTRUE);
    // task1->SetFillQAhistos(kFALSE);
    task1->SetProcessCharged(kTRUE);
    task1->SetProcessPID(kTRUE);
    task1->SetProcessPhi(kTRUE);
    task1->SetProcessV0s(kTRUE);
    // Flow
    task1->SetFlowRFPsPtMin(0.2);
    task1->SetFlowRFPsPtMax(5.);
    task1->SetFlowDoFourCorrelations(kFALSE);
    task1->SetFlowFillWeights(kFALSE);
    // task1->SetUseWeigthsFile("alien:///alice/cern.ch/user/v/vpacik/weights_preliminary_16q.root",kFALSE);
    // Events selection
    // task1->SetUseAliEventCuts();
    task1->SetTrigger(0);
    task1->SetColisionSystem(AliAnalysisTaskUniFlow::kPPb);
    task1->SetMultEstimator("V0A");
    task1->SetPeriod(AliAnalysisTaskUniFlow::k16q);
    task1->SetPVtxZMax(10);
    // Charged selection
    task1->SetChargedEtaMax(0.8);
    // task1->SetChargedPtMin(0.2);
    // task1->SetChargedPtMax(5.);
    // task1->SetChargedDCAzMax(0.1);
    // task1->SetChargedDCAxyMax(0.2);
    task1->SetChargedNumTPCclsMin(70);
    task1->SetChargedTrackFilterBit(96);
    // PID selection
    task1->SetPIDUseAntiProtonOnly(kFALSE);
    task1->SetPIDNumSigmasPionMax(3);
    task1->SetPIDNumSigmasKaonMax(3);
    task1->SetPIDNumSigmasProtonMax(3);
    task1->SetUseBayesPID(kTRUE);
    task1->SetPIDBayesProbPionMin(0.95);
    task1->SetPIDBayesProbKaonMin(0.85);
    task1->SetPIDBayesProbProtonMin(0.85);
    // Phi selection
    // task1->SetPhiMotherEtaMax(0.8);
    // V0 selection cuts
    task1->SetV0sOnFly(kFALSE);
    task1->SetV0sTPCRefit(kTRUE);
    task1->SetV0sRejectKinks(kTRUE);
    task1->SetV0sUseCrossMassRejection(kTRUE);
    task1->SetV0sCrossMassCutK0s(0.005);
    task1->SetV0sCrossMassCutLambda(0.020);
    task1->SetV0sDCAPVMin(0.06);
    task1->SetV0sDCAPVMax(0.);
    // task1->SetV0sDCAPVzMax(1.);
    // task1->SetV0sDaughtersFilterBit(211);
    task1->SetV0sDCADaughtersMax(1.);
    task1->SetV0sDecayRadiusMin(0.5);
    task1->SetV0sDecayRadiusMax(0.);
    task1->SetV0sDaughterPtMin(0.1);
    task1->SetV0sDaughterEtaMax(0.8);
    task1->SetV0sMotherEtaMax(0.8);
    task1->SetV0sMotherRapMax(0.);
    task1->SetV0sK0sCPAMin(0.97);
    task1->SetV0sLambdaCPAMin(0.99);
    task1->SetV0sK0sNumTauMax(5);
    task1->SetV0sK0sArmenterosAlphaMin(0.2);
    task1->SetV0sLambdaNumTauMax(3.8);
    // task1->SetV0sK0sKaonNumTPCSigmaMax(3.);
    // task1->SetV0sLambdaPionNumTPCSigmaMax(3.);
    // task1->SetV0sLambdaProtonNumTPCSigmaMax(3.);

    /*
    AliAnalysisTaskUniFlow* task2 = AddTaskUniFlow("UniFlow_FB768_Nsigma");
    // tracks & event selection cuts
    task2->SetAODAnalysis(kTRUE);
    task2->SetPbPbAnalysis(kTRUE);
    task2->SetPPbAnalysis(kFALSE);
    task2->SetPPAnalysis(kFALSE);
    task2->SetPeriod10h(kFALSE);
    task2->SetCentFlag(0);
    task2->SetPVtxZMax(10);
    task2->SetTrackEtaMax(0.8);
    task2->SetTrackPtMin(0.2);
    task2->SetTrackPtMax(5.);
    task2->SetNumTPCclsMin(70);
    task2->SetTrackFilterBit(768);
    task2->SetPionNumSigmasMax(3);
    task2->SetKaonNumSigmasMax(3);
    task2->SetProtonNumSigmasMax(3);
    //task2->SetDoFlow(kTRUE);
    task2->SetDiffFlow(kTRUE);
    task2->SetPID(kTRUE);
    task2->SetDoV0s(kTRUE);
    task2->SetSampling(kTRUE);
    task2->SetDoFlowGenFramKatarina(kTRUE);
    task2->SetDoOldFlow(kFALSE);
    task2->SetUseBayesPID(kFALSE);
    task2->SetPIDBayesProbPionMin(0.95);
    task2->SetPIDBayesProbKaonMin(0.95);
    task2->SetPIDBayesProbProtonMin(0.95);
    // V0 selection cuts
    task2->SetV0sOnFly(kFALSE);
    task2->SetV0sTPCRefit(kTRUE);
    task2->SetV0sRejectKinks(kTRUE);
    task2->SetV0sDCAPVMin(0.1);
    task2->SetV0sDCAPVMax(0.);
    task2->SetV0sDCADaughtersMax(1.);
    task2->SetV0sDecayRadiusMin(5.);
    task2->SetV0sDecayRadiusMax(0.);
    task2->SetV0sDaughterPtMin(0.1);
    task2->SetV0sDaughterEtaMax(0.8);
    task2->SetV0sMotherEtaMax(0.8);
    task2->SetV0sMotherRapMax(0.);
    task2->SetV0sMotherPtMin(0.2);
    task2->SetV0sMotherPtMax(10.);
    task2->SetV0sK0sCPAMin(0.998);
    task2->SetV0sLambdaCPAMin(0.998);
    task2->SetV0sK0sNumTauMax(3.);
    task2->SetV0sK0sArmenterosAlphaMin(0.2);
    task2->SetV0sLambdaNumTauMax(3.);
    task2->SetV0sProtonNumSigmaMax(3.);
    task2->SetV0sProtonPIDPtMax(1.2);
    */
    if (!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    //mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 500);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        //chain->Add("~/NBI/Flow/flow/testData/2010/LHC10h/000138275/ESDs/pass2/AOD160/0803/AliAOD.root");
        // chain->Add("~/NBI/Flow/flow/testData/2016/LHC16l/000259888/pass1/AOD/015/AliAOD.root");
        chain->Add("~/NBI/Flow/data/2016/LHC16q/000265427/pass1_CENT_wSDD/AOD/001/AliAOD.root");
        mgr->StartAnalysis("local", chain); // start the analysis locally, reading the events from the TChain
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskUniFlow.cxx AliAnalysisTaskUniFlow.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskUniFlow.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20170820-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20160131-1");
        // select the input data
        alienHandler->SetGridDataDir(Form("/alice/data/2016/%s/",sPeriod.Data()));
        alienHandler->SetDataPattern("/pass1_CENT_wSDD/AOD/*/AliAOD.root");
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");
        // runnumber

        Int_t iNumRuns = sizeof(runNumber) / sizeof(runNumber[0]);

        for (Int_t i = 0; i < iNumRuns; i++)
        {
            //if (i == sizeof(runArray) / sizeof(runArray[1])) break;
            alienHandler->AddRunNumber(runNumber[i]);
        }

        alienHandler->SetMasterResubmitThreshold(90);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(5);
        alienHandler->SetExecutable("FlowPID.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(50000);
        alienHandler->SetJDLName("FlowPID.jdl");
        alienHandler->SetPrice(1);
        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);

        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(bMergeViaJDL);

        // define the output folders
        alienHandler->SetGridWorkingDir(sWorkDir.Data());
        alienHandler->SetGridOutputDir(sOutDir.Data());

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode(sGridMode.Data());
            mgr->StartAnalysis("grid");
        }
      }
}
