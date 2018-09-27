#ifdef __CLING__
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskUniFlow.h"
#endif

void runAnalysis()
{
    Bool_t local = 0; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 1; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)

    TString sGridMode = "full";
    //TString sGridMode = "terminate";

    Bool_t bMergeViaJDL = kTRUE;
    //Bool_t bMergeViaJDL = kFALSE;

    TString sWorkDir = "test";
    TString sOutDir = "output";

    // Pb-Pb Run2 5.02 TeV (Run2) : RunList_LHC15o_pass1_pidfix_CentralBarrelTracking_hadronPID_20161018_v0.txt [6 runs]
    TString sPeriod = "2015/LHC15o"; TString sPass = "pass1_pidfix"; Int_t runNumber[] = {
      245231, 245152, 245151, 245146, 245145
    };
    // Pb-Pb Run2 5.02 TeV (Run2) : RunList_LHC15o_pass1_CentralBarrelTracking_hadronPID_20161130_v6.txt [77 runs]
    // TString sPeriod = "2015/LHC15o"; TString sPass = "pass1"; Int_t runNumber[] = {
    //   246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246851, 246847,
    //   246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766,
    //   246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493,
    //   246488, 246487, 246434, 246431, 246424, 246276, 246275, 246272, 246271, 246225,
    //   246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151,
    //   246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042,
    //   246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923,
    //   245833, 245831, 245829, 245705, 245702, 245692, 245683
    // };

    #if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
    #else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    #endif

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskUniFlow");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE,kTRUE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"))));

    AliMultSelectionTask* taskMultSelection = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));
    taskMultSelection->SetSelectedTriggerClass(AliVEvent::kINT7);

    AliAnalysisTaskPIDResponse* taskPIDResponse = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"))));

    /*
    // PID response QA by ALICE
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa* taskPIDqa = AddTaskPIDqa("PID_QA.root");
    */

    #if !defined (__CINT__) || defined (__CLING__)
      gInterpreter->LoadMacro("AliAnalysisTaskUniFlow.cxx++g");
      AliAnalysisTaskUniFlow *task1 = reinterpret_cast<AliAnalysisTaskUniFlow*>(gInterpreter->ExecuteMacro("AddTaskUniFlow.C(\"UniFlow\")"));
    #else
      gROOT->LoadMacro("AliAnalysisTaskUniFlow.cxx++g");
      gROOT->LoadMacro("AddTaskUniFlow.C");
      AliAnalysisTaskUniFlow *task1 = AddTaskUniFlow("UniFlow");
    #endif

    // AliAnalysisTaskUniFlow* task1 = AddTaskUniFlow("UniFlow");
    // Analysis
    task1->AddTwo(2,-2);
    task1->AddTwoGap(2,-2,0.4);
    // task1->AddTwoGap(2,-2, 0.4);
    // task1->AddTwoGap(3,-3, 0.4);
    task1->AddFour(2,2,-2,-2);
    task1->AddFour(2,3,-2,-3,1,0);
    task1->AddFour(3,3,-3,-3,1,0);
    task1->AddFourGap(2,2,-2,-2,0.0);
    // task1->AddFourGap(2,3,-2,-3,0.0, FlowTask::kRFP);
    // task1->AddFourGap(3,3,-3,-3,0.0, FlowTask::kRFP);
    task1->AddFourGap(2,2,-2,-2,0.4);
    task1->AddFourGap(2,3,-2,-3,0.4,1,0);
    task1->AddFourGap(3,3,-3,-3,0.4,1,0);

    task1->AddThree(4,-2,-2, 0,1);
    task1->AddThree(5,-3,-2, 0,1);
    task1->AddThree(6,-3,-3, 0,1);
    // task1->AddThreeGap(4,-2,-2,0.0, FlowTask::kPOI);
    // task1->AddThreeGap(5,-3,-2,0.0, FlowTask::kPOI);
    // task1->AddThreeGap(6,-3,-3,0.0, FlowTask::kPOI);
    task1->AddThreeGap(4,-2,-2,0.4, 0,1);
    task1->AddThreeGap(5,-3,-2,0.4, 0,1);
    task1->AddThreeGap(6,-3,-3,0.4, 0,1);

    task1->SetAnalysisType(AliAnalysisTaskUniFlow::kAOD);
    task1->SetRunMode(AliAnalysisTaskUniFlow::kFull);
    task1->SetRejectAddPileUp(1);
    task1->SetNumEventsAnalyse(1);
    task1->SetMC(kFALSE);
    task1->SetSampling(0);
    task1->SetFillQAhistos(kTRUE);
    task1->SetProcessPID(kTRUE);
    task1->SetProcessPhi(kTRUE);
    task1->SetProcessV0s(kTRUE);
    // Flow
    // task1->SetFlowRFPsPtMin(0.2);
    // task1->SetFlowRFPsPtMax(5.0);
    task1->SetFlowFillWeights(kTRUE);
    // task1->SetUseWeigthsFile("alien:///alice/cern.ch/user/v/vpacik/weights-prel/weights_16l.root",kFALSE);
    // task1->SetUseWeigthsFile("./weights_16l.root",kTRUE);
    task1->SetUseWeights3D(kFALSE);
    // Events selection
    task1->SetTrigger(AliVEvent::kINT7);
    task1->SetCollisionSystem(AliAnalysisTaskUniFlow::kPbPb);
    task1->SetMultEstimator(AliAnalysisTaskUniFlow::kRFP);
    // task1->SetPVtxZMax(10);
    // Charged selection
    // task1->SetChargedEtaMax(0.8);
    // // task1->SetChargedPtMin(0.2);
    // // task1->SetChargedPtMax(5.);
    // // task1->SetChargedDCAzMax(0.1);
    // // task1->SetChargedDCAxyMax(0.2);
    // task1->SetChargedNumTPCclsMin(70);
    // task1->SetChargedTrackFilterBit(96);
    // // PID selection
    // task1->SetPIDUseAntiProtonOnly(kFALSE);
    // task1->SetPIDNumSigmasPionMax(3);
    // task1->SetPIDNumSigmasKaonMax(3);
    // task1->SetPIDNumSigmasProtonMax(3);
    // task1->SetUseBayesPID(kTRUE);
    // task1->SetPIDBayesProbPionMin(0.95);
    // task1->SetPIDBayesProbKaonMin(0.85);
    // task1->SetPIDBayesProbProtonMin(0.85);
    // Phi selection
    // task1->SetPhiMotherEtaMax(0.8);

    // // V0 selection cuts
    // task1->SetV0sOnFly(kFALSE);
    // task1->SetV0sTPCRefit(kTRUE);
    // task1->SetV0sRejectKinks(kTRUE);
    // task1->SetV0sUseCrossMassRejection(kTRUE);
    // task1->SetV0sCrossMassCutK0s(0.005);
    // task1->SetV0sCrossMassCutLambda(0.020);
    // task1->SetV0sDCAPVMin(0.06);
    // task1->SetV0sDCAPVMax(0.);
    // // task1->SetV0sDCAPVzMax(1.);
    // // task1->SetV0sDaughtersFilterBit(211);
    // task1->SetV0sDCADaughtersMax(1.);
    // task1->SetV0sDecayRadiusMin(0.5);
    // task1->SetV0sDecayRadiusMax(0.);
    // task1->SetV0sDaughterPtMin(0.1);
    // task1->SetV0sDaughterEtaMax(0.8);
    // task1->SetV0sMotherEtaMax(0.8);
    // task1->SetV0sMotherRapMax(0.);
    // task1->SetV0sK0sCPAMin(0.97);
    // task1->SetV0sLambdaCPAMin(0.99);
    // task1->SetV0sK0sNumTauMax(5);
    // task1->SetV0sK0sArmenterosAlphaMin(0.2);
    // task1->SetV0sLambdaNumTauMax(3.8);
    // // task1->SetV0sK0sKaonNumTPCSigmaMax(3.);
    // // task1->SetV0sLambdaPionNumTPCSigmaMax(3.);
    // // task1->SetV0sLambdaProtonNumTPCSigmaMax(3.);


    if (!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    //mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 500);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // chain->Add("~/Codes/ALICE/Flow/data/2016/LHC16l/000259888/pass1/AOD/001/AliAOD.root");
        chain->Add("~/Codes/ALICE/Flow/data/2015/LHC15o/000245232/pass1_pidfix/AOD194/0001/AliAOD.root");
        // chain->Add("~/NBI/Flow/data/2016/LHC16q/000265427/pass1_CENT_wSDD/AOD/001/AliAOD.root");
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
        alienHandler->SetAliPhysicsVersion("vAN-20180701-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20160131-1");
        // select the input data
        alienHandler->SetGridDataDir(Form("/alice/data/%s/",sPeriod.Data()));
        alienHandler->SetDataPattern(Form("/%s/AOD/*/AliAOD.root",sPass.Data()));
        // alienHandler->SetDataPattern("/pass1_CENT_wSDD/AOD/*/AliAOD.root");
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
        alienHandler->SetTTL(20000);
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
