#ifdef __CLING__
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskUniFlow.cxx"
#include "AliUniFlowCorrTask.cxx"
#endif

void runDev()
{
    Bool_t local = 0; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 0; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)

    TString sGridMode = "full";
    // TString sGridMode = "terminate";

    Bool_t bMergeViaJDL = kTRUE;
    // Bool_t bMergeViaJDL = kFALSE;

    TString sWorkDir = "dummy";
    TString sOutDir = "output";

    // Pb-Pb Run2 5.02 TeV (Run2) : RunList_LHC15o_pass1_CentralBarrelTracking_hadronPID_20161130_v6.txt [77 runs]
    TString sPeriod = "2015/LHC15o"; TString sPass = "pass1"; Int_t runNumber[] = {
      246994
      // 246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246851, 246847,
      // 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766,
      // 246765 ,246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493,
      // 246488, 246487, 246434, 246431, 246424, 246276, 246275, 246272, 246271, 246225,
      //,
      // 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151,
      // 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042,
      // 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923,
      // 245833, 245831, 245829, 245705, 245702, 245692, 245683
    };



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

    AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE,kFALSE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"))));

    AliMultSelectionTask* taskMultSelection = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));
    taskMultSelection->SetSelectedTriggerClass(AliVEvent::kINT7);

    AliAnalysisTaskPIDResponse* taskPIDResponse = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"))));

    /*
    // PID response QA by ALICE
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa* taskPIDqa = AddTaskPIDqa("PID_QA.root");
    */

    #if !defined (__CINT__) || defined (__CLING__)
      // printf("\n CLING \n\n");
      // gInterpreter->LoadMacro("AliAnalysisTaskUniFlow.cxx++g");
      // AliAnalysisTaskUniFlow *task1 = reinterpret_cast<AliAnalysisTaskUniFlow*>(gInterpreter->ExecuteMacro("AddTaskUniFlow.C(AliAnalysisTaskUniFlow::kPbPb,\"alien:///alice/cern.ch/user/v/vpacik/weights/lhc15o/6519/weights.root\")"));
      // AliAnalysisTaskUniFlow *task1 = reinterpret_cast<AliAnalysisTaskUniFlow*>(gInterpreter->ExecuteMacro("AddTaskUniFlow.C(AliAnalysisTaskUniFlow::kPbPb,\"weights.root\")"));
      AliAnalysisTaskUniFlow *task1 = reinterpret_cast<AliAnalysisTaskUniFlow*>(gInterpreter->ExecuteMacro("AddTaskUniFlow.C(AliAnalysisTaskUniFlow::kPbPb)"));
    #else
      // printf("\n CINT \n\n");
      gROOT->LoadMacro("AliAnalysisTaskUniFlow.cxx++g");
      gROOT->LoadMacro("AddTaskUniFlow.C");
      AliAnalysisTaskUniFlow *task1 = AddTaskUniFlow(AliAnalysisTaskUniFlow::kPbPb);
    #endif

    if(!task1) { printf("E-runPbPb: Task not initialised!\n"); return; }

    // AliAnalysisTaskUniFlow* task1 = AddTaskUniFlow("UniFlow");
    // Analysis
    task1->SetRunMode(AliAnalysisTaskUniFlow::kTest);
    task1->SetNumEventsAnalyse(10);
    task1->SetSampling(0);
    task1->SetFillQAhistos(0);
    task1->SetProcessPID(1);
    task1->SetProcessPhi(1);
    task1->SetProcessV0s(1);
    task1->SetCentrality(AliAnalysisTaskUniFlow::kRFP,0,3000,90);
    // task1->SetAddCentCut(AliAnalysisTaskUniFlow::kV0M,0,10);
    // task1->SetFlowPOIsPtBins({1.0,4.0}, AliAnalysisTaskUniFlow::kK0s);
    // task1->SetFlowPOIsPtBins({2.0,3.0}, AliAnalysisTaskUniFlow::kLambda);
    // task1->SetFlowPOIsPtBins({0.1,2.0}, AliAnalysisTaskUniFlow::kPhi);
    // task1->SetFlowPOIsPtBins({1.0,3.0,3.2,5.}, AliAnalysisTaskUniFlow::kCharged);
    // // weigths
    task1->SetFlowFillWeights(0);
    task1->SetFlowFillAfterWeights(1);
    task1->SetUseWeigthsRunByRun(0);
    task1->SetUseWeights3D(kFALSE);
    // correlations
    task1->AddCorr({2,-2}, {0.0});
    // task1->AddCorr({3,-3}, {0.0});

    // task1->AddCorr({4,-2,-2}, {}, 0,1);
    // task1->AddCorr({5,-3,-2}, {}, 0,1);
    // task1->AddCorr({6,-3,-3}, {}, 0,1);

    // task1->AddCorr({2,2,-2,-2}, {});
    // task1->AddCorr({2,3,-2,-3}, {}, 1,0);
    // task1->AddCorr({3,3,-3,-3}, {});
    //
    // task1->AddCorr({2,-2}, {0.0});
    // task1->AddCorr({3,-3}, {0.0});
    //
    // task1->AddCorr({4,-2,-2}, {0.0}, 0,1);
    // task1->AddCorr({5,-3,-2}, {0.0}, 0,1);
    // task1->AddCorr({6,-3,-3}, {0.0}, 0,1);
    //
    // task1->AddCorr({2,2,-2,-2}, {0.0});
    // task1->AddCorr({2,3,-2,-3}, {0.0}, 1,0);
    // task1->AddCorr({3,3,-3,-3}, {0.0});
    //
    // task1->AddCorr({2,-2}, {0.4});
    // task1->AddCorr({3,-3}, {0.4});
    //
    // task1->AddCorr({4,-2,-2}, {0.4}, 0,1);
    // task1->AddCorr({5,-3,-2}, {0.4}, 0,1);
    // task1->AddCorr({6,-3,-3}, {0.4}, 0,1);
    //
    // task1->AddCorr({2,2,-2,-2}, {0.4});
    // task1->AddCorr({2,3,-2,-3}, {0.4}, 1,0);
    // task1->AddCorr({3,3,-3,-3}, {0.4});
    //
    // task1->AddCorr({2,-2}, {0.8});
    // task1->AddCorr({3,-3}, {0.8});
    //
    // task1->AddCorr({4,-2,-2}, {0.8}, 0,1);
    // task1->AddCorr({5,-3,-2}, {0.8}, 0,1);
    // task1->AddCorr({6,-3,-3}, {0.8}, 0,1);
    //
    // task1->AddCorr({2,2,-2,-2}, {0.8});
    // task1->AddCorr({2,3,-2,-3}, {0.8}, 1,0);
    // task1->AddCorr({3,3,-3,-3}, {0.8});


    if (!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    //mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 500);

    TStopwatch watch = TStopwatch();
    watch.Stop();
    watch.Reset();

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        // chain->Add("~/Codes/Flow/data/2016/LHC16l/000259888/pass1/AOD/001/AliAOD.root");
        // chain->Add("~/Codes/Flow/data/2016/LHC16q/000265427/pass1_CENT_wSDD/AOD/001/AliAOD.root");
        chain->Add("~/Codes/Flow/data/2015/LHC15o/000246153/pass1/AOD194/0002/AliAOD.root");
        watch.Start();
        mgr->StartAnalysis("local", chain); // start the analysis locally, reading the events from the TChain
        watch.Stop();
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("PWGCFFLOWGF.par libPWGEMCALbase.so");
        // alienHandler->SetAnalysisSource("AliUniFlowCorrTask_cxx.so AliAnalysisTaskUniFlow_cxx.so");
        // alienHandler->SetAnalysisSource("AliUniFlowCorrTask.cxx  AliAnalysisTaskUniFlow.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20190819_ROOT6-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20160131-1");
        // select the input data
        alienHandler->SetGridDataDir(Form("/alice/data/%s",sPeriod.Data()));
        alienHandler->SetDataPattern(Form("%s/AOD194/*/AliAOD.root",sPass.Data()));
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
        alienHandler->SetSplitMaxInputFileNumber(100);
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
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode(sGridMode.Data());
        }
        watch.Start();
        mgr->StartAnalysis("grid");
        watch.Stop();
      }

      printf("============================\n");
      printf("  Running for: %f (CPU) | %f (real)\n", watch.CpuTime(), watch. RealTime());
      printf("============================\n");
}
