/* ProcessUniFlow class
 *
 * Class implemented for processing results of AliAnalysisTaskUniFlow task.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */
#include <vector>
#include "TROOT.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

void SetStyle(Bool_t graypalette=kFALSE);
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
void DrawLogo (Int_t logo=0, Double_t xmin =  0.28, Double_t ymin= 0.68) ;
void FakeHistosOnlyForExample(TH1*&hstat, TH1*&hsyst, TH1*&hsystCorr);
void LoadLibs();


class FlowTask
{
  public:
    enum    PartSpecies {kUnknown=0, kRefs, kCharged, kPion, kKaon, kProton, kK0s, kLambda, kPhi}; // list of all particle species of interest

                FlowTask(PartSpecies species = kUnknown, const char* name = "");
                ~FlowTask(); // default destructor

    void        PrintTask(); // listing values of internal properties

    TString     GetSpeciesName();
    TString     GetEtaGapString() { return TString(Form("%02.2g",10*fEtaGap)); }

    void        SetHarmonics(Int_t harm) { fHarmonics = harm; }
    void        SetEtaGap(Float_t eta) { fEtaGap = eta; }
    void        SetNumSamples(Short_t num) { fNumSamples = num; }
    void        SetInputTag(const char* name) { fInputTag = name; }
    void        SetPtBins(Double_t* array, const Short_t size); // setup the pt binning for this task, where size is number of elements in array
    void        SetShowMultDist(Bool_t show) { fShowMult = show; }
    void        SetRebinning(Bool_t rebin = kTRUE) { fRebinning = rebin; }
    void        SetMergePosNeg(Bool_t merge = kTRUE) {fMergePosNeg = merge; }
    void        SetDesamplingUseRMS(Bool_t use = kTRUE) { fDesampleUseRMS = use; }
    void        SuggestPtBinning(Bool_t bin = kTRUE, Double_t entries = 20000) { fSuggestPtBins = bin; fSuggestPtBinEntries = entries; } // suggest pt binning based on number of candidates

    // fitting
    void        SetFlowFitFixTerms(Bool_t fix = kTRUE) { fFlowFitFixTerms = fix; }
    void        SetFlowPhiSubtLS(Bool_t sub = kTRUE) { fFlowPhiSubtLS = sub; }
    void        SetFittingRange(Double_t low = 0., Double_t high = 0.) { if(low > 0. && high > 0. && low < high) { fFlowFitRangeLow = low; fFlowFitRangeHigh = high; } else printf("E-SetFittingRange: One of the edges unspecified: values not set."); }
    void        SetFittingRejectNumSigmas(UShort_t sigmas) { fFlowFitRejectNumSigmas = sigmas; }
    void        SetInvMassRebin(Short_t rebin = 2) { fRebinInvMass = rebin; }
    void        SetFlowMassRebin(Short_t rebin = 2) { fRebinFlowMass = rebin; }
    void        SetAlexFitting(Bool_t fit = kTRUE) { fFitByAlex = fit; }
  protected:
  private:

    TString     fTaskTag; // "unique" tag used primarily for storing output
    TString     fName; // task name
    PartSpecies fSpecies; // species involved
    TString     fInputTag; // alterinative tag appended to name of input histos & profiles
    Int_t       fHarmonics; // harmonics
    Double_t    fEtaGap; // eta gap
    Short_t     fNumSamples; // [10] number of samples
    static const Short_t fNumPtBinsMax = 100; // initialization (maximum) number of pt bins
    Double_t    fPtBinsEdges[fNumPtBinsMax]; // pt binning
    Short_t     fNumPtBins; // actual number of pT bins (not size of array) for rebinning
    Bool_t      fShowMult; // show multiplicity distribution
    Bool_t      fSampleMerging; // [kFALSE] flag for merging TProfiles (good for refs)
    Bool_t      fRebinning; // [kTRUE] flag for rebinning prior to desampling
    Bool_t      fDesampleUseRMS; // [kFALSE] flag for using RMS as uncertainty during desampling
    Bool_t      fMergePosNeg; // [kFALSE] flag for merging results corresponding to positive and negative POIs
    Bool_t      fSuggestPtBins; // suggest pt binning
    Double_t    fSuggestPtBinEntries; // suggest pt binning
    // Reconstructed fitting
    Bool_t      fFlowPhiSubtLS; // [kFALSE] flag for subtraction of like-sign background from the unlike-sign one
    Bool_t      fFlowFitFixTerms; // [kTRUE] flag for using sequential flow extraction (fixing all terms in master formula)
    Double_t    fFlowFitRangeLow; // lower edge for fitting during flow extraction
    Double_t    fFlowFitRangeHigh; // high edge for fitting during flow extraction
    UShort_t    fFlowFitRejectNumSigmas; // number of sigmas for rejection peak region
    Short_t     fRebinInvMass; // flag for rebinning inv-mass (and BG) histo
    Short_t     fRebinFlowMass; // flag for rebinning flow-mass profile
    Bool_t      fFitByAlex; // flag for fitting according to Alexandru procedure

    std::vector<TH1D*>* fVecHistInvMass; // container for sliced inv. mass projections
    std::vector<TH1D*>* fVecHistInvMassBG; // container for sliced inv. mass projections for BG candidates (phi)
    std::vector<TH1D*>* fVecHistFlowMass; // container for sliced flow-mass projections
    TCanvas*     fCanvas; // temporary canvas for mass plotting
};

//_____________________________________________________________________________
FlowTask::FlowTask(PartSpecies species, const char* name)
{
  fName = name;
  fSpecies = species;
  fHarmonics = 0;
  fEtaGap = 0;
  fInputTag = "";
  fNumSamples = 10;
  fNumPtBins = -1;
  fShowMult = kFALSE;
  fRebinning = kTRUE;
  fSampleMerging = kFALSE;
  fDesampleUseRMS = kFALSE;
  fMergePosNeg = kFALSE;
  fSuggestPtBins = kFALSE;
  fSuggestPtBinEntries = 20000;
  fFlowFitFixTerms = kTRUE;
  fFlowPhiSubtLS = kFALSE;
  fFitByAlex = kFALSE;
  fFlowFitRangeLow = 0;
  fFlowFitRangeHigh = 0;
  fFlowFitRejectNumSigmas = 0;
  fRebinFlowMass = 0;
  fRebinInvMass = 0;
  fVecHistInvMass = new std::vector<TH1D*>;
  fVecHistInvMassBG = new std::vector<TH1D*>;
  fVecHistFlowMass = new std::vector<TH1D*>;

  fTaskTag = this->GetSpeciesName();
  if(!fName.EqualTo("")) fTaskTag.Append(Form("_%s",fName.Data()));
}
//_____________________________________________________________________________
FlowTask::~FlowTask()
{
  if(fVecHistFlowMass) delete fVecHistFlowMass;
  if(fVecHistInvMass) delete fVecHistInvMass;
  if(fVecHistInvMassBG) delete fVecHistInvMassBG;
  if(fCanvas) delete fCanvas;
}
//_____________________________________________________________________________
void FlowTask::SetPtBins(Double_t* array, const Short_t size)
{
  if(size < 0 || size > fNumPtBinsMax) { Error("Wrong size of pt binning array.","SetPtBins"); return; }
  if(!array) { Error("Wrong array.","SetPtBins"); return; }

  fNumPtBins = size - 1;

  for(Short_t i(0); i < size; i++)
  {
    fPtBinsEdges[i] = array[i];
  }

  return;
}
//_____________________________________________________________________________
TString FlowTask::GetSpeciesName()
{
  TString name = TString();
  switch (fSpecies)
  {
    case kRefs : name.Append("Refs"); break;
    case kCharged : name.Append("Charged"); break;
    case kPion : name.Append("Pion"); break;
    case kKaon : name.Append("Kaon"); break;
    case kProton : name.Append("Proton"); break;
    case kPhi : name.Append("Phi"); break;
    case kK0s : name.Append("K0s"); break;
    case kLambda : name.Append("Lambda"); break;
    default: name.Append("Unknown");
  }
  return name;
}
//_____________________________________________________________________________
void FlowTask::PrintTask()
{
  printf("----- Printing task info ------\n");
  printf("   fTaskTag: \"%s\"\n",fTaskTag.Data());
  printf("   fName: \"%s\"\n",fName.Data());
  printf("   fInputTag: \"%s\"\n",fInputTag.Data());
  printf("   fSpecies: %s (%d)\n",GetSpeciesName().Data(),fSpecies);
  printf("   fHarmonics: %d\n",fHarmonics);
  printf("   fEtaGap: %g\n",fEtaGap);
  printf("   fShowMult: %s\n", fShowMult ? "true" : "false");
  printf("   fSuggestPtBins: %s\n", fSuggestPtBins ? "true" : "false");
  printf("   fMergePosNeg: %s\n", fMergePosNeg ? "true" : "false");
  printf("   fFlowFitRangeLow: %g\n",fFlowFitRangeLow);
  printf("   fFlowFitRangeHigh: %g\n",fFlowFitRangeHigh);
  printf("   fFlowFitRejectNumSigmas: %u\n",fFlowFitRejectNumSigmas);
  printf("   fNumPtBins: %d (limit %d)\n",fNumPtBins,fNumPtBinsMax);
  if(fNumPtBins > -1) { printf("   fPtBinsEdges: "); for(Short_t i(0); i < fNumPtBins+1; i++) printf("%g ",fPtBinsEdges[i]); printf("\n"); }
  printf("------------------------------\n");
  return;
}
//#############################################################################
//----------------- End of FlowTask class -------------------------------------
//#############################################################################

class ProcessUniFlow
{
  public:

                ProcessUniFlow(); // default constructor
                ~ProcessUniFlow(); // default destructor

    void        SetInputFilePath(const char* path) { fsInputFilePath = path; }
    void        SetInputFileName(const char* name) { fsInputFileName = name; }
    void        SetOutputFilePath(const char* path) { fsOutputFilePathRoot = path; }
    void        SetOutputFileName(const char* name) { fsOutputFileName = name; }
    void        SetOutputFileMode(const char* mode = "RECREATE") { fsOutputFileMode = mode; }
    void        SetTaskName(const char* name) { fsTaskName = name; }
    void        SetGlobalProfNameLabel(const char* label = "") { fsGlobalProfNameLabel = label; } // add global profile label for all tasks NOTE: for the purpose of Flow sub
    void        SetSaveMult(Bool_t bSave = kTRUE) { fbSaveMult = bSave; } // save reference multiplicity
    void        SetMultiplicityBins(Double_t* array, const Short_t size); // setup the global multiplicity binning, where size is number of elements in array
    void        SetFitCumulants(Bool_t cum = kTRUE) { fFlowFitCumulants = cum; } // use cn{2} vs m_inv instead of vn{2} vs. m_inv
    void        SetDebug(Bool_t debug = kTRUE) { fbDebug = debug; }
    void        AddTask(FlowTask* task = 0x0); // add task to internal lists of all tasks
    void        Run(); // running the task (main body of the class)
    void        Clear(); // clearing (removing tasks, etc.) after running
  protected:

  private:
    Bool_t      Initialize(); // initialization task
    Bool_t      LoadLists(); // loading flow lists from input file

    void        ProcessTask(FlowTask* task = 0x0); // process FlowTask according to it setting
    Bool_t      ProcessRefs(FlowTask* task = 0x0); // process reference flow task
    Bool_t      ProcessDirect(FlowTask* task = 0x0, Short_t iMultBin = 0); // process PID (pion,kaon,proton) flow task
    Bool_t      ProcessReconstructed(FlowTask* task = 0x0, Short_t iMultBin = 0); // process  V0s flow
    Bool_t      PrepareSlices(const Short_t multBin, FlowTask* task = 0x0, TProfile3D* p3Cor = 0x0, TH3D* h3Entries = 0x0, TH3D* h3EntriesBG = 0x0); // prepare
    Bool_t 	    ExtractFlowPhi(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowK0s(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass); // extract flow via flow-mass method for K0s candidates
		Bool_t 	    ExtractFlowLambda(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass); // extract flow via flow-mass method for Lambda candidates

    Bool_t 	    ExtractFlowK0sAlex(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass); // extract flow via flow-mass method for K0s candidates

    void        SuggestMultBinning(const Short_t numFractions);
    void        SuggestPtBinning(TH3D* histEntries = 0x0, TProfile3D* profFlowOrig = 0x0, FlowTask* task = 0x0, Short_t binMult = 0); //
    TProfile*   MergeTProfile(TList* list); // Merge TProfiles in a list and return it
    TH1D*       DesampleList(TList* list = 0x0, FlowTask* task = 0x0, Short_t iMultBin = 0); // Desample list of samples for estimating the uncertanity
    TH1D*       Desampling(TList* list, TH1D* hMerged, FlowTask* task, Short_t iMultBin = 0); // (new-refactored) desampling procedure
    TH1D*       TestRebin(TH1D* hOrig = 0x0, FlowTask* task = 0x0); // testing desample - manual rebin

    void        TestProjections(); // testing projection of reconstructed particles
    TProfile2D* Project3DProfile(const TProfile3D* prof3dorig = 0x0); // making projection out of TProfile3D
    TProfile2D* DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,bool originalRange, bool useUF, bool useOF) const;
    TH2D*       DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY, bool computeErrors, bool originalRange, bool useUF, bool useOF) const;


    // printing output methods
    void        Fatal(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Error(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Warning(TString sMsg, TString sMethod = ""); // printf the msg as warning
    void        Info(TString sMsg, TString sMethod = ""); // printf the msg as info
    void        Debug(TString sMsg, TString sMethod = ""); // printf the msg as info

    static const Short_t fiNumMultBinsGlobal = 200; // global initialization number of bins
    Double_t    fdMultBins[fiNumMultBinsGlobal]; // global multiplicity/centrality binning
    Short_t     fiNumMultBins; // number of multiplicity bins (not size of array)

    TString     fsInputFilePath; // path to the input folder with input file
    TString     fsInputFileName; // name of input file
    TString     fsOutputFilePath; // current (working) path to the output folder (process/task dependent)
    TString     fsOutputFilePathRoot; // root (global/top-level) path to the output folder
    TString     fsOutputFileName; // name of output file
    TString     fsOutputFileMode; // [RECREATE] mode of output file
    TString     fsTaskName; // name of task (inchluded in data structure names)
    TString     fsOutputFileFormat; // [pdf] format of output files (pictures)
    TString     fsGlobalProfNameLabel; // global profile label for all task
    Bool_t      fbSaveMult; // [kFALSE]
    Bool_t      fFlowFitCumulants; // [kFALSE]

    Bool_t      fbInit; // flag for initialization status
    Bool_t      fbDebug; // flag for debugging : if kTRUE Debug() messages are displayed
    TFile*      ffInputFile; //! input file container
    TFile*      ffOutputFile; //! output file container
    TFile*      ffDesampleFile; //! output file for results of desampling
    TList*      flFlowRefs; //! TList from input file with RFPs flow profiles
    TList*      flFlowCharged; //! TList from input file with Charged flow profiles
    TList*      flFlowPID; //! TList from input file with PID (pi,K,p) flow profiles
    TList*      flFlowPhi; //! TList from input file with Phi flow profiles
    TList*      flFlowK0s; //! TList from input file with K0s flow profiles
    TList*      flFlowLambda; //! TList from input file with Lambda flow profiles
    TList*      flQACharged; //! TList from input file with Charged QA plots / profiles
    TList*      flQAPID; //! TList from input file with PID (pi,K,p) QA plots / profiles
    TList*      flQAPhi; //! TList from input file with Phi QA plots / profiles
    TList*      flQAV0s; //! TList from input file with K0s QA plots / profiles
    std::vector<FlowTask*> fvTasks; // vector of task for individual species proccesing

};
//_____________________________________________________________________________
ProcessUniFlow::ProcessUniFlow() :
  fbDebug(kFALSE),
  fbInit(kFALSE),
  fbSaveMult(kFALSE),
  fFlowFitCumulants(kFALSE),
  ffInputFile(0x0),
  ffOutputFile(0x0),
  ffDesampleFile(0x0),
  flFlowRefs(0x0),
  flFlowCharged(0x0),
  flFlowPID(0x0),
  flFlowPhi(0x0),
  flFlowK0s(0x0),
  flFlowLambda(0x0),
  flQACharged(0x0),
  flQAPID(0x0),
  flQAPhi(0x0),
  flQAV0s(0x0)
{
  // default constructor
  fsInputFilePath = TString("");
  fsInputFileName = TString("AnalysisResults.root");
  fsOutputFilePathRoot = TString("");
  fsOutputFilePath = TString("");
  fsOutputFileName = TString("UniFlow.root");
  fsOutputFileMode = TString("RECREATE");
  fsTaskName = TString("UniFlow");
  fsOutputFileFormat = TString("pdf");
  fsGlobalProfNameLabel = TString("");
  fvTasks = std::vector<FlowTask*>();

  for(Short_t i(0); i < fiNumMultBinsGlobal; i++) fdMultBins[i] = 0;
}
//_____________________________________________________________________________
ProcessUniFlow::~ProcessUniFlow()
{
  // default destructor
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;

  if(flFlowRefs) delete flFlowRefs;
  if(flFlowCharged) delete flFlowCharged;
  if(flFlowPID) delete flFlowPID;
  if(flFlowPhi) delete flFlowPhi;
  if(flFlowK0s) delete flFlowK0s;
  if(flFlowLambda) delete flFlowLambda;

  // deleting the FlowTasks
  const Short_t iNumTasks = fvTasks.size();
  for(Short_t index(0); index < iNumTasks; index++)
  {
    delete fvTasks.at(index);
  }
}
//_____________________________________________________________________________
void ProcessUniFlow::Clear()
{
  Info("Cleaning ProcessUniFlow instance","Clear");
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;

  flFlowRefs = 0x0;
  flFlowCharged = 0x0;
  flFlowPID = 0x0;
  flFlowPhi = 0x0;
  flFlowK0s = 0x0;
  flFlowLambda = 0x0;

  const Short_t iNumTasks = fvTasks.size();
  for(Short_t index(0); index < iNumTasks; index++)
  {
    if(fvTasks.at(index)) delete fvTasks.at(index);
  }
  Info("Cleaning done!","Clear");
}
//_____________________________________________________________________________
void ProcessUniFlow::Run()
{
  gStyle->SetOptFit(1100);

  // main body of the class
  if(!Initialize()) { Fatal("Task not initialized","Run"); return; }

  const Short_t iNumTasks = fvTasks.size();

  Info("===== Running over tasks ======","Run");
  Info(Form("  Number of tasks: %d\n",iNumTasks),"Run");
  for(Short_t iTask(0); iTask < iNumTasks; iTask++)
  {
    FlowTask* currentTask = fvTasks.at(iTask);
    if(!currentTask) continue;
    ProcessTask(currentTask);
  }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::Initialize()
{
  // initialization of all necessery prerequisits
  Info("Initializating task","Initialize");

  // opening input file
  ffInputFile = new TFile(Form("%s/%s",fsInputFilePath.Data(),fsInputFileName.Data()),"READ");
  if(!ffInputFile || !ffInputFile->IsOpen())
  {
    Fatal(Form("Input file %s/%s not open",fsInputFilePath.Data(),fsInputFileName.Data()),"Initialize");
    return kFALSE;
  }

  // Setting current output path
  fsOutputFilePath = fsOutputFilePathRoot;

  // checking specified output folder & required sub-folders
  gSystem->mkdir(fsOutputFilePath.Data(),kTRUE);

  // opening output file
  ffOutputFile = TFile::Open(Form("%s/%s",fsOutputFilePath.Data(),fsOutputFileName.Data()),fsOutputFileMode.Data());
  if(!ffOutputFile || !ffOutputFile->IsOpen())
  {
    Fatal(Form("Output file %s/%s not open",fsOutputFilePath.Data(),fsOutputFileName.Data()),"Initialize");
    return kFALSE;
  }

  // creating output file for Desampling
  ffDesampleFile = TFile::Open(Form("%s/desampling.root",fsOutputFilePath.Data()),"RECREATE");
  if(!ffDesampleFile) { Fatal(Form("Output desampling file '%s/desampling.root' not open!","Initialize")); return kFALSE; }

  Info("Files loaded","Initialize");

  if(!LoadLists()) return kFALSE;
  Info("Flow lists loaded","Initialize");

  // initialization succesfull
  Info("Initialization succesfull","Initialize");
  return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::SetMultiplicityBins(Double_t* array, const Short_t size)
{
  if(size < 0 || size > fiNumMultBinsGlobal) { Fatal("Wrong number of bins selected","SetMultiplicityBins"); return; }
  if(!array) { Fatal("Multiplicity array not valid","SetMultiplicityBins"); return; }

  fiNumMultBins = size-1;

  for(Short_t i(0); i < size; i++)
  {
    fdMultBins[i] = array[i];
  }
  return;
}
//_____________________________________________________________________________
TProfile* ProcessUniFlow::MergeTProfile(TList* list)
{
  if(!list) { Error("Input list does not exists!","MergeTProfile"); return 0x0; }

  TProfile* merged = (TProfile*) list->At(0)->Clone();
  merged->Reset();
  Double_t mergeStatus = merged->Merge(list);
  if(mergeStatus == -1.0) { Error("Merging unsuccesfull","MergeTProfile"); return 0x0; }

  return merged;
}
//_____________________________________________________________________________
void ProcessUniFlow::SuggestMultBinning(const Short_t numFractions)
{
  if(numFractions < 1) { Error("Suggested number of fractions is too low","SuggestMultBinning"); return; }
  if(!fbInit) Initialize();

  // loading multiplicity dist.
  ffInputFile->cd(fsTaskName.Data());
  gDirectory->ls();

  TList* lQACharged = (TList*) gDirectory->Get(Form("QA_Charged_%s",fsTaskName.Data()));
  if(!lQACharged) return;

  lQACharged->ls();
  TH1D* hMult = (TH1D*) lQACharged->FindObject("fhQAChargedMult_After");

  TCanvas* canSuggestMult = new TCanvas("canSuggestMult","SuggestMultBinning",400,400);
  canSuggestMult->cd(1);
  gPad->SetLogy();
  hMult->Draw();

  Double_t dCumul = 0;
  Double_t dContent = 0;
  Double_t dWeight = 0;
  Double_t dSum = 0;
  Double_t dSumWeight = 0;

  for(Short_t bin(1); bin < hMult->GetNbinsX()+1; bin++)
  {
    dContent = hMult->GetBinContent(bin);
    dWeight = 1/(bin - 1);
    dSum += dContent*(dWeight); // weighed average
    dSumWeight += dWeight;
  }

  Double_t dEntries = hMult->GetEntries();
  // Double_t dSegment = dEntries / numFractions;
  Double_t dSegment = dSum / dSumWeight / numFractions;

  printf("entries %g | Sum %g | sumweight %g | segmend %g\n",dEntries,dSum,dSumWeight,dSegment);

  for(Short_t bin(1); bin < hMult->GetNbinsX()+1; bin++)
  {
    dContent = hMult->GetBinContent(bin);
    // dWeight = bin - 1;
    dWeight = 1/(bin-1);
    dCumul += dContent*(dWeight);
    // dSumWeight += dWeight;
    if(dCumul >= dSegment)
    {
      printf("bin %d : mult %g (cumul %g)\n",bin,hMult->GetBinLowEdge(bin),dCumul);
      dCumul = 0;
    }
  }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::LoadLists()
{
  // loading TLists into task
  if(!ffInputFile) { Fatal("Input file does not exists!","LoadLists"); return kFALSE; }
  ffInputFile->cd(fsTaskName.Data());

  flFlowRefs = (TList*) gDirectory->Get(Form("Flow_Refs_%s",fsTaskName.Data()));
  if(!flFlowRefs) { Fatal("flFlow_Refs list does not exists!","LoadLists"); ffInputFile->ls(); return kFALSE; }
  flFlowCharged = (TList*) gDirectory->Get(Form("Flow_Charged_%s",fsTaskName.Data()));
  if(!flFlowCharged) { Fatal("flFlow_Charged list does not exists!","LoadLists"); return kFALSE; }
  flFlowPID = (TList*) gDirectory->Get(Form("Flow_PID_%s",fsTaskName.Data()));
  if(!flFlowPID) { Fatal("flFlow_PID list does not exists!","LoadLists"); return kFALSE; }
  flFlowPhi = (TList*) gDirectory->Get(Form("Flow_Phi_%s",fsTaskName.Data()));
  if(!flFlowPhi) { Fatal("flFlow_Phi list does not exists!","LoadLists"); return kFALSE; }
  flFlowK0s = (TList*) gDirectory->Get(Form("Flow_K0s_%s",fsTaskName.Data()));
  if(!flFlowK0s) { Fatal("flFlow_K0s list does not exists!","LoadLists"); return kFALSE; }
  flFlowLambda = (TList*) gDirectory->Get(Form("Flow_Lambda_%s",fsTaskName.Data()));
  if(!flFlowLambda) { Fatal("flFlow_Lambda list does not exists!","LoadLists"); return kFALSE; }

  flQACharged = (TList*) gDirectory->Get(Form("QA_Charged_%s",fsTaskName.Data()));
  if(!flQACharged) { Fatal("flQACharged list does not exists!","LoadLists"); return kFALSE; }
  flQAPID = (TList*) gDirectory->Get(Form("QA_PID_%s",fsTaskName.Data()));
  if(!flQAPID) { Fatal("flQAPID list does not exists!","LoadLists"); return kFALSE; }
  flQAPhi = (TList*) gDirectory->Get(Form("QA_Phi_%s",fsTaskName.Data()));
  if(!flQAPhi) { Fatal("flQAPhi list does not exists!","LoadLists"); return kFALSE; }
  flQAV0s = (TList*) gDirectory->Get(Form("QA_V0s_%s",fsTaskName.Data()));
  if(!flQAV0s) { Fatal("flQAV0s list does not exists!","LoadLists"); return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::ProcessTask(FlowTask* task)
{
  fsOutputFilePath = Form("%s/%s",fsOutputFilePathRoot.Data(),task->fTaskTag.Data());
  gSystem->mkdir(fsOutputFilePath.Data(),kTRUE);

  Info(Form("Processing task: %s",task->fName.Data()),"ProcessTask");
  if(!task) { Error("Task not valid!","ProcessTask"); return; }

  task->PrintTask();

  Bool_t bProcessed = kFALSE;
  switch (task->fSpecies)
  {
    case FlowTask::kRefs:
      bProcessed = ProcessRefs(task);
      break;

    case FlowTask::kCharged:
    case FlowTask::kPion:
    case FlowTask::kKaon:
    case FlowTask::kProton:
      for(Short_t binMult(0); binMult < fiNumMultBins; binMult++) { bProcessed = ProcessDirect(task,binMult); }
      break;

    case FlowTask::kPhi:
    case FlowTask::kK0s:
    case FlowTask::kLambda:
      for(Short_t binMult(0); binMult < fiNumMultBins; binMult++) { bProcessed = ProcessReconstructed(task,binMult); }
      break;
    default: break;
  }

  if(!bProcessed) { Error(Form("Task '%s' not processed correctly!",task->fName.Data()),"ProcessTask"); return; }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessRefs(FlowTask* task)
{
  Info("Processing Refs task","ProcesRefs");
  if(!task) { Error("Task not valid!","ProcessRefs"); return kFALSE; }
  if(task->fSpecies != FlowTask::kRefs) { Error("Task species not kRefs!","ProcessRefs"); return kFALSE; }

  // preparing for desampling
  TProfile* prof = 0x0;
  TProfile* profRebin  = 0x0;
  TH1D* histProj = 0x0;
  TList* list = new TList();
  list->SetOwner(kTRUE);
  TList* listMerge = new TList();
  listMerge->SetOwner(kTRUE);

  // rebinning <multiplicity>
  if(fbSaveMult)
  {
    TProfile* profMult = (TProfile*) flQACharged->FindObject(Form("fpRefsMult"));
    if(!profMult) { Error("MeanMult profile not found!"); flFlowRefs->ls(); return kFALSE; }
    TProfile* profMult_rebin = (TProfile*) profMult->Rebin(fiNumMultBins,Form("%s_rebin",profMult->GetName()),fdMultBins);

    ffOutputFile->cd();
    profMult_rebin->Write(profMult_rebin->GetName());
  }


  for(Short_t i(0); i < task->fNumSamples; i++)
  {
    TProfile* prof = (TProfile*) flFlowRefs->FindObject(Form("fpRefs_%s<2>_harm%d_gap%02.2g_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,i));
    if(!prof) { Warning(Form("Profile sample %d does not exits. Skipping",i),"ProcesRefs"); return kFALSE; }

    if(task->fRebinning)
    {
      // rebinning the profiles
      profRebin = (TProfile*) prof->Rebin(fiNumMultBins,Form("%s_rebin",prof->GetName()),fdMultBins);
      histProj = profRebin->ProjectionX(Form("%s_proj",profRebin->GetName()));
      listMerge->Add(profRebin);
    }
    else
    {
      // no rebinning
      histProj = prof->ProjectionX(Form("%s_proj",prof->GetName()));
      listMerge->Add(prof);
    }
    list->Add(histProj);
  }

  // merging samples in listMerge to get merged histo
  // merging all samples together (NOTE: good for Refs)
  TProfile* merged = (TProfile*) listMerge->At(0)->Clone();
  merged->Reset();
  Double_t mergeStatus = merged->Merge(listMerge);
  merged->SetName(Form("hCum2_%s_harm%d_gap%s_merged",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));
  if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessRefs"); return kFALSE; }

  TH1D* hMerged = (TH1D*) merged->ProjectionX();
  hMerged->SetName(Form("hFlow2_%s_harm%d_gap%s_merged",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  // desampling (similarly to PID & Charged tracks)
  TH1D* hDesampled = Desampling(list, hMerged, task);
  if(!hDesampled) { Error("Desampling unsuccesfull","ProcessRefs"); return kFALSE; }
  hDesampled->SetName(Form("hCum2_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));
  hDesampled->SetTitle(Form("%s c_{%d}{2} | Gap %s ",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));


  // // comparing old and new method
  // TH1D* hDesampledList = DesampleList(list,task);
  // TCanvas* canTest = new TCanvas("canTest","canTest");
  // canTest->Divide(3,1);
  // canTest->cd(1);
  // hMerged->DrawCopy();
  // canTest->cd(2);
  // hDesampled->DrawCopy();
  // canTest->cd(3);
  // hDesampledList->DrawCopy();

  // estimating vn out of cn
  TH1D* hDesampledFlow = (TH1D*) hDesampled->Clone(Form("%s_flow",hDesampled->GetName()));
  hDesampledFlow->SetName(Form("hFlow2_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));
  hDesampledFlow->SetTitle(Form("%s v_{%d}{2} | Gap %s ",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  Double_t dContent = 0., dError = 0.;
  for(Short_t iBin(1); iBin < hDesampled->GetNbinsX()+1; iBin++)
  {
    dContent = hDesampled->GetBinContent(iBin);
    dError = hDesampled->GetBinError(iBin);

    if(dContent > 0. && dError >= 0.)
    {
      hDesampledFlow->SetBinContent(iBin, TMath::Sqrt(dContent));
      hDesampledFlow->SetBinError(iBin, TMath::Sqrt( TMath::Power(dError,2) / (4*dContent) ));
    }
    else
    {
      hDesampledFlow->SetBinContent(iBin, 99.9);
      hDesampledFlow->SetBinError(iBin, 99.9);
    }
  }

  // saving to output file
  ffOutputFile->cd();
  hDesampled->Write(Form("%s",hDesampled->GetName()));
  hDesampledFlow->Write(Form("%s",hDesampledFlow->GetName()));

  // TCanvas* canTest = new TCanvas("canTest","canTest");
  // canTest->Divide(2,1);
  // canTest->cd(1);
  // hDesampled->Draw();
  // canTest->cd(2);
  // hDesampledFlow->Draw();

  if(list) delete list;
  if(listMerge) delete listMerge;
  // if(hDesampled) delete hDesampled;
  // if(hDesampledFlow) delete hDesampledFlow;
  // if(profRebin) delete profRebin;
  // if(histProj) delete histProj;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessDirect(FlowTask* task, Short_t iMultBin)
{
  Info("Processing direct task","ProcesDirect");
  if(!task) { Error("Task not valid!","ProcessDirect"); return kFALSE; }

  TList* listInput = 0x0;
  switch (task->fSpecies)
  {
    case FlowTask::kCharged:
      listInput = flFlowCharged;
      break;

    case FlowTask::kPion:
    case FlowTask::kKaon:
    case FlowTask::kProton:
      listInput = flFlowPID;
      break;

    default:
      Error("Task species not direct!","ProcessDirect");
      return kFALSE;
  }

  // preparing vn' samples
  TList* listFlow = new TList();
  TList* listCum = new TList();
  TList* listMerge = new TList();
  TProfile2D* prof2 = 0x0;
  TProfile2D* prof2pos = 0x0;
  TProfile2D* prof2neg = 0x0;
  TProfile* prof = 0x0;
  TProfile* profRebin = 0x0;
  TProfile* profRef = 0x0;
  TProfile* profRefRebin = 0x0;
  TH1D* hFlow = 0x0;

  Short_t binMultLow = 0;
  Short_t binMultHigh = 0;
  Double_t dRef = 0;

  for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
  {
    if(task->fMergePosNeg)
    {
      // loading Pos & Neg if fMergePosNeg is ON
      prof2pos = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_Pos_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample));
      prof2neg = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_Neg_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample));

      if(!prof2pos || !prof2neg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessDirect"); return kFALSE; }
      listMerge->Clear();
      listMerge->Add(prof2pos);
      listMerge->Add(prof2neg);

      prof2 = (TProfile2D*) listMerge->At(0)->Clone();
      prof2->Reset();
      Double_t mergeStatus = prof2->Merge(listMerge);
      if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessDirect"); return kFALSE; }
    }
    else
    {
      // loading single profile
      if(task->fInputTag.EqualTo(""))
      { prof2 = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_Pos_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample)); }
      else
      { prof2 = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%02.2g_%s_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data(),iSample)); }
    }
    if(!prof2) { Error(Form("Profile sample %d does not exists.",iSample),"ProcessDirect"); return kFALSE; }

    // preparing Refs
    profRef = (TProfile*) flFlowRefs->FindObject(Form("fpRefs_%s<2>_harm%d_gap%02.2g_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,iSample));
    if(!profRef) { Error(Form("Profile sample %d does not exists",iSample),"ProcesPID"); return kFALSE; }

    profRefRebin = (TProfile*) profRef->Rebin(fiNumMultBins,Form("%s_rebin",profRef->GetName()),fdMultBins);
    dRef = profRefRebin->GetBinContent(iMultBin+1);
    profRefRebin->Draw();
    // printf("dRef %g\n",dRef);
    // NOTE: complains about Sumw2

    // rebinning according in mult bin
    binMultLow = prof2->GetXaxis()->FindFixBin(fdMultBins[iMultBin]);
    binMultHigh = prof2->GetXaxis()->FindFixBin(fdMultBins[iMultBin+1]) - 1;
    prof = prof2->ProfileY(Form("%s_cent%d",prof2->GetName(),iMultBin),binMultLow,binMultHigh);

    if(task->fNumPtBins > 0)
    {
      // rebinning according to pt bins
      profRebin = (TProfile*) prof->Rebin(task->fNumPtBins,Form("%s_rebin",prof->GetName()),task->fPtBinsEdges);
    }
    else
    {
      profRebin = (TProfile*) prof->Clone(Form("%s_rebin",prof->GetName()));
    }



    // making TH1D projection (to avoid handling of bin entries)
    TH1D* hCum = (TH1D*) profRebin->ProjectionX();
    hCum->SetName(Form("%s_Cum",prof->GetName()));
    listCum->Add(hCum);

    hFlow = (TH1D*) profRebin->ProjectionX();
    // printf("Pre %g\n",hFlow->GetBinContent(20));
    // printf("Scale %g\n",1/TMath::Sqrt(dRef));
    hFlow->SetName(Form("%s_Flow",prof->GetName()));
    hFlow->Scale(1/TMath::Sqrt(dRef));
    // printf("Post %g\n",hFlow->GetBinContent(20));
    listFlow->Add(hFlow);
  }
  delete listMerge;

  Debug(Form("Number of samples in list pre merging %d",listFlow->GetEntries()),"ProcessDirect");

  TH1D* hDesampled = DesampleList(listFlow,task,iMultBin);
  if(!hDesampled) { Error("Desampling unsuccesfull","ProcessDirect"); return kFALSE; }

  hDesampled->SetName(Form("hFlow2_%s_harm%d_gap%s_cent%d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
  hDesampled->SetTitle(Form("%s v_{%d}{2} | Gap %s | Cent %d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));

  // desampling cumulants
  TH1D* hDesampled_Cum = DesampleList(listCum,task,iMultBin);
  if(!hDesampled_Cum) { Error("Desampling unsuccesfull","ProcessDirect"); return kFALSE; }

  hDesampled_Cum->SetName(Form("hCum2_%s_harm%d_gap%s_cent%d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
  hDesampled_Cum->SetTitle(Form("%s v_{%d}{2} | Gap %s | Cent %d",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));


  // saving to output file
  ffOutputFile->cd();
  hDesampled->Write();
  hDesampled_Cum->Write();

  delete listFlow;
  delete prof;
  delete profRebin;
  delete profRefRebin;
  delete hFlow;
  delete hDesampled;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessReconstructed(FlowTask* task,Short_t iMultBin)
{
  Info("Processing task","ProcessReconstructed");
  if(!task) { Error("Task not valid!","ProcessReconstructed"); return kFALSE; }
  // if(task->fNumPtBins < 1) { Error("Num of pt bins too low!","ProcessReconstructed"); return kFALSE; }


  TList* listMerge = 0x0;
  // preparing particle dependent variables for switch
  //  -- input histos / profiles with entries and correlations
  TH3D* histEntries = 0x0;
  TH3D* histEntriesPos = 0x0;
  TH3D* histEntriesNeg = 0x0;
  TH3D* histBG = 0x0; // entries for BG (phi)
  TH3D* histBGPos = 0x0; // entries for BG (phi)
  TH3D* histBGNeg = 0x0; // entries for BG (phi)
  TProfile3D* profFlow = 0x0;
  TProfile3D* profFlowPos = 0x0;
  TProfile3D* profFlowNeg = 0x0;
  //  -- naming variables
  TString sSpeciesName; // in objects name
  TString sSpeciesLabel; // LaTeX for titles

  // checking particles species and assigning particle dependent variables
  switch (task->fSpecies)
  {
    case FlowTask::kPhi :
      sSpeciesName = TString("Phi");
      sSpeciesLabel = TString("#phi");

      if(task->fMergePosNeg)
      {
        // loading Pos & Neg if fMergePosNeg is ON
        // merging profiles
        profFlowPos = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        profFlowNeg = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_%s<2>_harm%d_gap%02.2g_Neg",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));

        if(!profFlowPos || !profFlowNeg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessDirect"); return kFALSE; }
        listMerge = new TList();
        listMerge->Add(profFlowPos);
        listMerge->Add(profFlowNeg);

        profFlow = (TProfile3D*) listMerge->At(0)->Clone();
        profFlow->Reset();
        Double_t mergeStatus = profFlow->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        // merging histos
        histEntriesPos = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        histEntriesNeg = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_Neg",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        if(!histEntriesPos || !histEntriesNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histEntriesPos);
        listMerge->Add(histEntriesNeg);

        histEntries = (TH3D*) listMerge->At(0)->Clone();
        histEntries->Reset();
        mergeStatus = histEntries->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        histBGPos = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_Pos",10*task->fEtaGap));
        histBGNeg = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_Neg",10*task->fEtaGap));
        if(!histBGPos || !histBGNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histBGPos);
        listMerge->Add(histBGNeg);

        histBG = (TH3D*) listMerge->At(0)->Clone();
        histBG->Reset();
        mergeStatus = histBG->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;
      }
      else
      {
        // loading single profile
        if(task->fInputTag.EqualTo(""))
        {
          histEntries = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
          histBG = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_Pos",10*task->fEtaGap));
          profFlow = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_<2>_%sharm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        }
        else
        {
          histEntries = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesSignal_%sgap%02.2g_%s",fsGlobalProfNameLabel.Data(),10*task->fEtaGap,task->fInputTag.Data()));
          histBG = (TH3D*) flFlowPhi->FindObject(Form("fh3PhiEntriesBG_gap%02.2g_%s",10*task->fEtaGap,task->fInputTag.Data()));
          profFlow = (TProfile3D*) flFlowPhi->FindObject(Form("fp3PhiCorr_%s<2>_harm%d_gap%02.2g_%s",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data()));
        }
      }
    break;

    case FlowTask::kK0s :
      sSpeciesName = TString("K0s");
      sSpeciesLabel = TString("K^{0}_{S}");

      if(task->fMergePosNeg)
      {
        // loading Pos & Neg if fMergePosNeg is ON
        // merging TProfiles
        profFlowPos = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        profFlowNeg = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_Neg",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        if(!profFlowPos || !profFlowNeg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(profFlowPos);
        listMerge->Add(profFlowNeg);

        profFlow = (TProfile3D*) listMerge->At(0)->Clone();
        profFlow->Reset();
        Double_t mergeStatus = profFlow->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging profiles unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        // merging histos
        histEntriesPos = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        histEntriesNeg = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_Neg",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        if(!histEntriesPos || !histEntriesNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histEntriesPos);
        listMerge->Add(histEntriesNeg);

        histEntries = (TH3D*) listMerge->At(0)->Clone();
        histEntries->Reset();
        mergeStatus = histEntries->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;
      }
      else
      {
        // loading single profile
        if(task->fInputTag.EqualTo(""))
        {
          histEntries = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
          profFlow = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        }
        else
        {
          histEntries = (TH3D*) flFlowK0s->FindObject(Form("fh3V0sEntriesK0s_%sgap%02.2g_%s",fsGlobalProfNameLabel.Data(),10*task->fEtaGap,task->fInputTag.Data()));
          profFlow = (TProfile3D*) flFlowK0s->FindObject(Form("fp3V0sCorrK0s_%s<2>_harm%d_gap%02.2g_%s",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data()));
        }
      }
    break;

    case FlowTask::kLambda :
      sSpeciesName = TString("Lambda");
      sSpeciesLabel = TString("#Lambda/#bar{#Lambda}");

      if(task->fMergePosNeg)
      {
        // loading Pos & Neg if fMergePosNeg is ON
        // merging profiles
        profFlowPos = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        profFlowNeg = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_Neg",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));

        if(!profFlowPos || !profFlowNeg) { Error("Pos OR Neg profile not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }
        listMerge = new TList();
        listMerge->Add(profFlowPos);
        listMerge->Add(profFlowNeg);

        profFlow = (TProfile3D*) listMerge->At(0)->Clone();
        profFlow->Reset();
        Double_t mergeStatus = profFlow->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;

        // merging histos
        histEntriesPos = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        histEntriesNeg = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_Neg",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
        if(!histEntriesPos || !histEntriesNeg) { Error("Pos OR Neg histo not found for Pos&Neg merging.","ProcessReconstructed"); return kFALSE; }

        listMerge = new TList();
        listMerge->Add(histEntriesPos);
        listMerge->Add(histEntriesNeg);

        histEntries = (TH3D*) listMerge->At(0)->Clone();
        histEntries->Reset();
        mergeStatus = histEntries->Merge(listMerge);
        if(mergeStatus == -1) { Error("Merging histos unsuccesfull","ProcessReconstructed"); return kFALSE; }
        delete listMerge;
      }
      else
      {
        // loading single profile
        if(task->fInputTag.EqualTo(""))
        {
          histEntries = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_Pos",fsGlobalProfNameLabel.Data(),10*task->fEtaGap));
          profFlow = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_Pos",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap));
        }
        else
        {
          histEntries = (TH3D*) flFlowLambda->FindObject(Form("fh3V0sEntriesLambda_%sgap%02.2g_%s",fsGlobalProfNameLabel.Data(),10*task->fEtaGap,task->fInputTag.Data()));
          profFlow = (TProfile3D*) flFlowLambda->FindObject(Form("fp3V0sCorrLambda_%s<2>_harm%d_gap%02.2g_%s",fsGlobalProfNameLabel.Data(),task->fHarmonics,10*task->fEtaGap,task->fInputTag.Data()));
        }
      }
    break;

    default:
      Error("Task species not V0s nor Phi!","ProcessReconstructed");
      return kFALSE;
  }

  if(!histEntries) { Error("Entries histos not found!","ProcessReconstructed"); return kFALSE; }
  if(!profFlow) { Error("Cumulant histos not found!","ProcessReconstructed"); return kFALSE; }

  // check if suggest pt binning flag is on if of Pt binning is not specified
  if(task->fSuggestPtBins || task->fNumPtBins < 1)
  {
    SuggestPtBinning(histEntries,profFlow,task,iMultBin);
  }

  if(task->fNumPtBins < 1) { Error("Num of pt bins too low!","ProcessReconstructed"); return kFALSE; }

  // task->PrintTask();

  if(!PrepareSlices(iMultBin,task,profFlow,histEntries,histBG)) return kFALSE;

  TH1D* hFlow = 0x0;
  if(!fFlowFitCumulants)
  {
    hFlow = new TH1D(Form("hFlow2_%s_harm%d_gap%02.2g_cent%d",sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin),Form("%s: v_{%d}{2 | Gap %g} (%g - %g); #it{p}_{T} (GeV/#it{c})",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1]), task->fNumPtBins,task->fPtBinsEdges);
  }
  else
  {
    hFlow = new TH1D(Form("hCum2_%s_harm%d_gap%02.2g_cent%d",sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin),Form("%s: d_{%d}{2 | Gap %g} (%g - %g); #it{p}_{T} (GeV/#it{c})",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1]), task->fNumPtBins,task->fPtBinsEdges);
  }

  TH1D* hInvMass = 0x0;
  TH1D* hInvMassBG = 0x0;
  TH1D* hFlowMass = 0x0;
  Double_t dFlow = 0, dFlowError = 0; // containers for flow extraction results
  TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1600,1200); // canvas for fitting results

  TCanvas* canFlowAll = new TCanvas("canFlowAll","canFlowAll",1600,1200);
  TCanvas* canInvMassAll = new TCanvas("canInvMassAll","canInvMassAll",1600,1200);
  canFlowAll->Divide(3,ceil(task->fNumPtBins/3.));
  canInvMassAll->Divide(3,ceil(task->fNumPtBins/3.));

  TLatex* latex = new TLatex();
  latex->SetTextSize(0.12);
  latex->SetNDC();

  TLatex* latex2 = new TLatex();
  // latex2->SetTextSize(0.12);
  latex2->SetNDC();

  for(Short_t binPt(0); binPt < task->fNumPtBins; binPt++)
  {
    hInvMass = task->fVecHistInvMass->at(binPt);
    hFlowMass = task->fVecHistFlowMass->at(binPt);

    hInvMass->SetTitle(Form("%s: InvMass dist (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
    hFlowMass->SetTitle(Form("%s: FlowMass (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));

    hInvMass->SetMarkerStyle(kFullCircle);
    hFlowMass->SetMarkerStyle(kFullCircle);

    // extracting flow
    switch (task->fSpecies)
    {
      case FlowTask::kPhi :
        hInvMassBG = task->fVecHistInvMassBG->at(binPt);
        if( !ExtractFlowPhi(task,hInvMass,hInvMassBG,hFlowMass,dFlow,dFlowError,canFitInvMass) ) { Warning("Flow extraction unsuccesfull","ProcessReconstructed"); return kFALSE; }
      break;

      case FlowTask::kK0s :
        if(task->fFitByAlex)
        {
          if( !ExtractFlowK0sAlex(task,hInvMass,hFlowMass,dFlow,dFlowError,canFitInvMass) ) { Warning("Flow extraction unsuccesfull (Alex's procedure)","ProcessReconstructed"); return kFALSE; }
        }
        else
        {
          if( !ExtractFlowK0s(task,hInvMass,hFlowMass,dFlow,dFlowError,canFitInvMass) ) { Warning("Flow extraction unsuccesfull","ProcessReconstructed"); return kFALSE; }
        }
      break;

      case FlowTask::kLambda :
        if( !ExtractFlowLambda(task,hInvMass,hFlowMass,dFlow,dFlowError,canFitInvMass) ) { Warning("Flow extraction unsuccesfull","ProcessReconstructed"); return kFALSE; }
      break;

      default :
        Error("Uknown species","ProcessReconstructed");
        return kFALSE;
    }

    gSystem->mkdir(Form("%s/fits/",fsOutputFilePath.Data()));

    canFitInvMass->cd(1);
    // if(task->fSpecies == FlowTask::kPhi) canFitInvMass->cd(2);
    latex2->DrawLatex(0.18,0.58,Form("pt %g-%g cent %g-%g%%",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));

    canFitInvMass->SaveAs(Form("%s/fits/Fit_%s_n%d2_gap%02.2g_cent%d_pt%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,binPt,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

    canFlowAll->cd(binPt+1);
    hFlowMass->DrawCopy();
    latex->DrawLatex(0.18,0.73,Form("pt%g-%g cent%g-%g",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));

    canInvMassAll->cd(binPt+1);
    hInvMass->DrawCopy();
    latex->DrawLatex(0.18,0.73,Form("pt%g-%g cent%g-%g",task->fPtBinsEdges[binPt],task->fPtBinsEdges[binPt+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));


    if(TMath::Abs(dFlow) > 1 )
    {
      hFlow->SetBinContent(binPt+1,0);
      hFlow->SetBinError(binPt+1,0);
    }
    else
    {
      hFlow->SetBinContent(binPt+1,dFlow);
      hFlow->SetBinError(binPt+1,dFlowError);
    }

  } // endfor {binPt}

  // task->fCanvas->Draw();
  canFlowAll->SaveAs(Form("%s/FlowMassFits_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
  canInvMassAll->SaveAs(Form("%s/InvMassFits_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());


  TCanvas* cFlow = new TCanvas("cFlow","cFlow");
  cFlow->cd();
  hFlow->SetStats(0);
  hFlow->Draw();
  cFlow->SaveAs(Form("%s/Flow_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

  ffOutputFile->cd();
  hFlow->Write();

  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::DesampleList(TList* list, FlowTask* task, Short_t iMultBin)
{
  if(!task) { Error("FlowTask does not exists","DesampleList"); return 0x0; }
  if(!list) { Error("List does not exists","DesampleList"); return 0x0; }
  if(list->GetEntries() < 1) { Error("List is empty","DesampleList"); return 0x0; }
  if(list->GetEntries() != task->fNumSamples) { Warning("Number of list entries is different from task number of samples","DesampleList"); }

  Debug(Form("Number of samples in list pre-desampling: %d",list->GetEntries()),"DesampleList");

  TH1D* hTempSample = (TH1D*) list->At(0);
  TH1D* hDesampled = (TH1D*) hTempSample->Clone(Form("%s_Desampled",hTempSample->GetName()));
  hDesampled->Reset();

  Double_t content = 0;
  Double_t error = 0;

  std::vector<Double_t> vecContents;
  std::vector<Double_t> vecErrors;
  Double_t dMean = 0;
  Double_t dRMS = 0;

  Double_t dSum = 0;
  Double_t dW = 0;
  Double_t dAverage = 0;
  Double_t dAve_err = 0;

  for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++)
  {
    dSum = 0;
    dW = 0;
    dAverage = 9.;
    dAve_err = 9.;

    vecContents.clear();
    vecErrors.clear();

    for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
    {
      hTempSample = (TH1D*) list->At(iSample);
      if(!hTempSample) { Warning(Form("Sample %d not found! Skipping!",iSample),"DesampleList"); continue; }

      content = hTempSample->GetBinContent(bin);
      error = hTempSample->GetBinError(bin);

      vecContents.push_back(content);
      vecErrors.push_back(error);

      if(error <= 0.) continue;

      dSum += content / TMath::Power(error,2);
      dW += 1 / TMath::Power(error,2);
      Debug(Form("Sample: %d | bin %d | %g +- %g",iSample,bin,content,error),"DesampleList");
    }

    Debug(Form(" --- bin %d | Sum %g +- %g",bin,dSum,dW),"DesampleList");

    if(dSum == 0. && dW == 0.) continue; // skipping empty bins

    if(dW > 0.)
    {
      dAverage = dSum / dW;
      dAve_err = TMath::Sqrt(1/dW);
    }

    Debug(Form("W average | bin %d | %g +- %g",bin,dAverage,dAve_err),"DesampleList");

    // working on weighted mean and RMS
    if(vecContents.size() != vecErrors.size()) { Warning("Something went wrong: Vector with contents and error have not the same size","DesampleList"); }

    dSum = 0;
    dW = 0;

    for(ULong_t iSample(0); iSample < vecContents.size(); iSample++)
    {
      content = vecContents.at(iSample);
      error = vecErrors.at(iSample);

      // weighted mean
      // dW += TMath::Power(error, -2);
      // dSum += content * TMath::Power(error, -2);

      // non-weighted mean
      dW++;
      dSum += content;
    }

    dMean = dSum/dW;

    dSum = 0;
    for(ULong_t iSample(0); iSample < vecContents.size(); iSample++)
    {
      content = vecContents.at(iSample);
      // error = vecErrors.at(iSample);
      error = 1.;

      dSum += TMath::Power( (dMean - content) / error ,2);
    }
    dRMS = TMath::Sqrt(dSum / (dW-1));
    Debug(Form("Mean | bin %d | %g +- %g",bin,dMean,dRMS),"DesampleList");

    //ratio->SetBinContent(bin, dAverage / merged->GetBinContent(bin));
    //ratioErr->SetBinContent(bin, dAve_err / merged->GetBinError(bin));

    if(task->fDesampleUseRMS)
    {
      hDesampled->SetBinContent(bin,dMean);
      hDesampled->SetBinError(bin,dRMS);
      // hDesampled->SetBinError(bin,dRMS/TMath::Sqrt(task->fNumSamples));
    }
    else
    {
      hDesampled->SetBinContent(bin,dAverage);
      hDesampled->SetBinError(bin,dAve_err);
    }
    Debug(Form("Desampled: bin %d | %g +- %g\n",bin,hDesampled->GetBinContent(bin),hDesampled->GetBinError(bin)),"DesampleList");
  }

  // at this point, hDesampled is ready

  // getting copy which does not affect histo which is returned
  TH1D* hDesampledClone = (TH1D*) hDesampled->Clone(Form("%sClone",hDesampled->GetName()));

  TList* listOutput = new TList(); // list for collecting all QA histos

  // doing QA plots with spread, etc.
  TCanvas* canDesample = new TCanvas(Form("canDesample_%s",task->fName.Data()),Form("canDesample_%s",task->fName.Data()),1200,400);
  canDesample->Divide(3,1);

  TH1D* hTempRatio = 0x0;
  TH1D* hTempError = 0x0;

  TLine* lineUnity = new TLine();
  lineUnity->SetLineColor(kRed);
  lineUnity->SetLineWidth(3);

  canDesample->cd(1);
  hDesampledClone->SetStats(kFALSE);
  hDesampledClone->SetFillColor(kBlue);
  hDesampledClone->SetStats(kFALSE);
  hDesampledClone->SetMarkerStyle(20);
  hDesampledClone->SetMarkerSize(0.5);
  hDesampledClone->SetMarkerColor(kRed);
  hDesampledClone->DrawCopy("E2");

  for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
  {
    hTempSample = (TH1D*) list->At(iSample);
    if(!hTempSample) { Warning(Form("Sample %d not found during plotting QA! Skipping!",iSample),"DesampleList"); continue; }

    canDesample->cd(1);
    hTempSample->SetStats(kFALSE);
    hTempSample->SetLineColor(30+2*iSample);
    hTempSample->SetMarkerColor(30+2*iSample);
    hTempSample->SetMarkerStyle(24);
    hTempSample->SetMarkerSize(0.5);
    hTempSample->DrawCopy("hist p same");

    hTempRatio = (TH1D*) hTempSample->Clone(Form("%s_ratio",hTempSample->GetName()));
    hTempRatio->Divide(hDesampled);
    hTempRatio->SetYTitle("Value: final / sample");
    hTempRatio->SetTitleOffset(1.2,"Y");

    canDesample->cd(2);
    hTempRatio->SetMinimum(0.6);
    hTempRatio->SetMaximum(1.4);
    hTempRatio->Draw("hist p same");

    hTempError = (TH1D*) hTempSample->Clone(Form("%s_error",hTempSample->GetName()));
    for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hTempSample->GetBinError(bin)); }

    canDesample->cd(3);
    hTempError->SetMinimum(0.);
    hTempError->SetMaximum(1.5*hTempError->GetMaximum());
    hTempError->SetYTitle("Uncertainty");
    hTempError->SetTitleOffset(1.2,"Y");

    hTempError->Draw("hist p same");

    listOutput->Add(hTempSample);
    listOutput->Add(hTempRatio);
    listOutput->Add(hTempError);
  }

  canDesample->cd(1);
  hDesampledClone->DrawCopy("hist p same");

  canDesample->cd(2);
  lineUnity->DrawLine(hTempRatio->GetXaxis()->GetXmin(),1,hTempRatio->GetXaxis()->GetXmax(),1);

  hTempError = (TH1D*) hDesampledClone->Clone(Form("%s_error",hDesampled->GetName()));
  for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hDesampledClone->GetBinError(bin)); }
  listOutput->Add(hTempError);

  canDesample->cd(3);
  hTempError->Draw("hist p same");

  // saving QA plots
  canDesample->SaveAs(Form("%s/Desampling_%s_harm%d_gap%g_cent%d_%s.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,task->fName.Data(),fsOutputFileFormat.Data()));

  Info("Saving desampling QA into output file","DesampleList");
  ffDesampleFile->cd();
  listOutput->Add(canDesample);
  listOutput->Write(Form("Desampling_%s_cent%d_%s",task->GetSpeciesName().Data(),iMultBin,task->fName.Data()),TObject::kSingleKey);

  // deleting created stuff
  delete listOutput;
  // delete canDesample;
  delete lineUnity;
  // if(hTempSample) delete hTempSample;
  // delete hTempRatio;
  // delete hTempError;
  // delete hDesampledClone;

  return hDesampled;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::Desampling(TList* list, TH1D* hMerged, FlowTask* task, Short_t iMultBin)
// Calculate the uncertainties in each bin from sqrt of variance of individual saples and uses the central values from merged histogram
{
  if(!task) { Error("FlowTask does not exists","Desampling"); return 0x0; }
  if(!hMerged) { Error("Merged histo does not exists","Desampling"); return 0x0; }
  if(!list) { Error("List does not exists","Desampling"); return 0x0; }
  if(list->GetEntries() < 1) { Error("List is empty","Desampling"); return 0x0; }
  if(list->GetEntries() != task->fNumSamples) { Warning("Number of list entries is different from task number of samples","DesampleList"); }

  Debug(Form("Number of samples in list pre-desampling: %d",list->GetEntries()),"DesampleList");

  TH1D* hDesampled = (TH1D*) hMerged->Clone(Form("%s_Desampled",hMerged->GetName()));
  hDesampled->Reset(); // just to be sure

  // filling std::vectors with contents and errors
  for(Int_t iBin(1); iBin < hDesampled->GetNbinsX()+1; ++iBin)
  {
    Double_t dMergedMean = hMerged->GetBinContent(iBin);

    Double_t dVariance = 0.0;
    for(Int_t iSample(0); iSample < list->GetEntries(); ++iSample)
    {
      TProfile* profSample = (TProfile*) list->At(iSample);
      if(!profSample) { Error(Form("Sample %d not found!",iSample),"Desampling"); return 0x0; }

      Double_t dCon = profSample->GetBinContent(iBin);
      Double_t dErr = profSample->GetBinError(iBin);

      dVariance += TMath::Power(dCon - dMergedMean,2.0);
    }
    // getting variance
    dVariance = dVariance / list->GetEntries();
    // standard error on the mean
    // dVariance = dVariance / list->GetEntries();

    hDesampled->SetBinContent(iBin,dMergedMean);
    hDesampled->SetBinError(iBin,TMath::Sqrt(dVariance));
  }

  return hDesampled;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PrepareSlices(const Short_t multBin, FlowTask* task, TProfile3D* p3Cor, TH3D* h3Entries, TH3D* h3EntriesBG)
{

  if(!task) { Error("Input task not found!","PrepareSlices"); return kFALSE; }
  if(!h3Entries) { Error("Input hist with entries not found!","PrepareSlices"); return kFALSE; }
  if(!p3Cor) { Error("Input profile with correlations not found!","PrepareSlices"); return kFALSE; }
  if(multBin < 0 || multBin > fiNumMultBins) { Error("Wrong multiplicity bin index (not in range)!","PrepareSlices"); return kFALSE; }

  // cleaning the vectros with flow-mass and inv. mass plots
  if(task->fVecHistInvMass->size() > 0) task->fVecHistInvMass->clear();
  if(task->fVecHistInvMassBG->size() > 0) task->fVecHistInvMassBG->clear();
  if(task->fVecHistFlowMass->size() > 0) task->fVecHistFlowMass->clear();

  const Short_t binMultLow = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin]);
  const Short_t binMultHigh = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin+1]) - 1;
  // printf("Mult: %g(%d) -  %g(%d)\n",fdMultBins[multBin],binMultLow,fdMultBins[multBin+1],binMultHigh);

  Double_t dRefFlow = 1.;
  Double_t dRefFlowErr = 1.;

  if(!fFlowFitCumulants)
  {
    // loading reference flow, if not found, it will be prepared
    TH1D* hRefFlow = 0x0;
    hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
    if(!hRefFlow)
    {
      Warning("Relevant Reference flow not found within output file.","PrepareSlices");
      ffOutputFile->ls();
      // return kFALSE;

      Info("Creating relevant reference flow task.","PrepareSlices");
      FlowTask* taskRef = new FlowTask(FlowTask::kRefs,"Ref");
      taskRef->SetHarmonics(task->fHarmonics);
      taskRef->SetEtaGap(task->fEtaGap);
      taskRef->SetNumSamples(task->fNumSamples);
      taskRef->SetInputTag(task->fInputTag);
      if(ProcessRefs(taskRef))
      {
        hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
        if(!hRefFlow) {  Error("Automated Refs task completed, but RefFlow not found!","PrepareSlices"); return kFALSE; }
      }
      else { Error("Something went wrong when running automatic refs flow task:","PrepareSlices"); taskRef->PrintTask(); return kFALSE; }
    }
    dRefFlow = hRefFlow->GetBinContent(multBin+1);
    dRefFlowErr = hRefFlow->GetBinError(multBin+1);
    Debug(Form("Ref (bin %d): %g +- %g\n",multBin,dRefFlow,dRefFlowErr),"PrepareSlices");
  }

  // loop over pt
  Short_t binPtLow = 0;
  Short_t binPtHigh = 0;
  TH1D* hInvMass_temp = 0x0;
  TH1D* hInvMassBG_temp = 0x0;
  TProfile3D* prof3Flow_temp = 0x0;
  TProfile2D* prof2FlowMass_temp = 0x0;
  TProfile* profFlowMass_temp = 0x0;
  TH1D* hFlowMass_temp = 0x0;

  prof3Flow_temp = (TProfile3D*) p3Cor->Clone(Form("prof3Flow_temp_cent%d",multBin));
  prof3Flow_temp->GetXaxis()->SetRange(binMultLow,binMultHigh);
  prof2FlowMass_temp = Project3DProfile(prof3Flow_temp);

  Short_t iNumPtBins = task->fNumPtBins;

  TCanvas* canFlowMass = new TCanvas("canFlowMass","FlowMass",1400,600);
  TCanvas* canInvMass = new TCanvas("canInvMass","InvMass",1400,600);
  TCanvas* canInvMassBG = new TCanvas("canInvMassBG","InvMassBG",1400,600);
  canFlowMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canInvMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canInvMassBG->Divide(5,std::ceil(iNumPtBins/5)+1);

  Double_t dContent = 0, dError = 0;

  for(Short_t binPt(0); binPt < iNumPtBins; binPt++)
  {
    // estimating pt edges
    binPtLow = h3Entries->GetYaxis()->FindFixBin(task->fPtBinsEdges[binPt]);
    binPtHigh = h3Entries->GetYaxis()->FindFixBin(task->fPtBinsEdges[binPt+1]) - 1; // for rebin both bins are included (so that one needs to lower)
    printf("   Pt: %g(%d) -  %g(%d)\n",task->fPtBinsEdges[binPt],binPtLow,task->fPtBinsEdges[binPt+1],binPtHigh);

    // rebinning entries based on mult & pt binning
    hInvMass_temp = (TH1D*) h3Entries->ProjectionZ(Form("hInvMass_cent%d_pt%d",multBin,binPt),binMultLow,binMultHigh,binPtLow,binPtHigh,"e");
    if(h3EntriesBG) hInvMassBG_temp = (TH1D*) h3EntriesBG->ProjectionZ(Form("hInvMassBG_cent%d_pt%d",multBin,binPt),binMultLow,binMultHigh,binPtLow,binPtHigh,"e");

    // checking if rebinning inv mass hist
    if(task->fRebinInvMass > 1) { hInvMass_temp->Rebin(task->fRebinInvMass); if(h3EntriesBG){ hInvMassBG_temp->Rebin(task->fRebinInvMass); } }


    // projection of flow-mass profile
    profFlowMass_temp = (TProfile*) prof2FlowMass_temp->ProfileX(Form("profFlowMass_cent%d_pt%d",multBin,binPt),binPtLow,binPtHigh);

    // checking for rebinning the flow-mass profile
    if(task->fRebinFlowMass > 1) { profFlowMass_temp->Rebin(task->fRebinFlowMass); }

    hFlowMass_temp = (TH1D*) profFlowMass_temp->ProjectionX(Form("hFlowMass_cent%d_pt%d",multBin,binPt));
    // NOTE: this is the ONLY (for some freaking reason) way how to get proper TH1 wth <<2>> out of TProfile3D

    // scaling flow-mass with reference flow
    if(!fFlowFitCumulants)
    {
      for(Short_t bin(1); bin < hFlowMass_temp->GetNbinsX()+1; bin++)
      {
        dContent = hFlowMass_temp->GetBinContent(bin);
        dError = hFlowMass_temp->GetBinError(bin);

        if(dContent == 0. || dError == 0.) continue;

        hFlowMass_temp->SetBinContent(bin,dContent/dRefFlow);
        hFlowMass_temp->SetBinError(bin, TMath::Sqrt( TMath::Power(dError/dRefFlow,2) + TMath::Power(dRefFlowErr*dContent/(dRefFlow*dRefFlow),2) - 2*(dError*dRefFlowErr*dContent*TMath::Power(dRefFlow,-3))) );
        // printf("%g | %g \n", TMath::Sqrt( TMath::Power(dError/dRefFlow,2) + TMath::Power(dRefFlowErr*dContent/(dRefFlow*dRefFlow),2) - 2*(dError*dRefFlowErr*dContent*TMath::Power(dRefFlow,-3))), TMath::Sqrt( TMath::Power(dError/dRefFlow,2) + TMath::Power(dRefFlowErr*dContent/(dRefFlow*dRefFlow),2)) );
      }
    }

    // hInvMass_temp->SetTitle(Form("%s: Inv. Mass (|#Delta#eta| > %g)",sSpeciesLabel.Data(),0.05));
    // hFlowMass_temp->SetTitle("TEST");

    // ready to fitting
    task->fVecHistFlowMass->push_back(hFlowMass_temp);
    task->fVecHistInvMass->push_back(hInvMass_temp);
    if(h3EntriesBG) task->fVecHistInvMassBG->push_back(hInvMassBG_temp);

    canFlowMass->cd(binPt+1);
    hFlowMass_temp->Draw();

    canInvMass->cd(binPt+1);
    hInvMass_temp->Draw();

    if(h3EntriesBG)
    {
      canInvMassBG->cd(binPt+1);
      hInvMassBG_temp->Draw();
    }

  } // endfor {binPt}: over Pt bins

  printf(" # of slices: InvMass: %lu | InvMassBG %lu | FlowMass %lu\n",task->fVecHistInvMass->size(),task->fVecHistInvMassBG->size(),task->fVecHistFlowMass->size());

  gSystem->mkdir(Form("%s/slices/",fsOutputFilePath.Data()));
  canFlowMass->SaveAs(Form("%s/slices/Slices_FlowMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  canInvMass->SaveAs(Form("%s/slices/Slices_InvMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  if(h3EntriesBG) canInvMassBG->SaveAs(Form("%s/slices/Slices_InvMassBG_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));

  if(task->fVecHistInvMass->size() < 1 || task->fVecHistFlowMass->size() < 1 || task->fVecHistFlowMass->size() != task->fVecHistInvMass->size()) { Error("Output vector empty. Something went wrong","PrepareSlices"); return kFALSE; }
  if(h3EntriesBG && (task->fVecHistInvMassBG->size() < 1 || task->fVecHistInvMassBG->size() != task->fVecHistInvMass->size()) ) { Error("Output vector empty. Something went wrong with BG histograms","PrepareSlices"); return kFALSE; }
  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::TestRebin(TH1D* hOrig, FlowTask* task)
{
  if(!hOrig) { Error("Original histogram not found!","TestRebin"); return 0x0; }
  if(!task) { Error("Task not found!","TestRebin"); return 0x0; }
  TH1D* hRebin = (TH1D*) hOrig->Rebin(fiNumMultBins,Form("%s_testRebin",hOrig->GetName()),fdMultBins);

  // rebinning
  Short_t numBins = fiNumMultBins;
  Double_t* multBins = fdMultBins;
  Short_t binIndex = 0;


  const Short_t iNumBins = hOrig->GetNbinsX();

  Double_t dSum = 0;
  Double_t dSumWeights = 0;
  Double_t dContent = 0;
  Double_t dWeight = 0;
  for(Short_t mult = 1; mult < iNumBins+1; mult++)
  {
    if(hOrig->GetBinLowEdge(mult) < multBins[binIndex]) continue;

    dContent = hOrig->GetBinContent(mult);
    dWeight =  TMath::Power(hOrig->GetBinError(mult),-2);

    dSumWeights += dWeight;
    dSum += dContent*dWeight;

    if( hOrig->GetBinLowEdge(mult+1) == multBins[binIndex+1] )
    {
      hRebin->SetBinContent(binIndex+1, dSum / dSumWeights);
      hRebin->SetBinError(binIndex+1, TMath::Sqrt(1/dSumWeights));

      dSumWeights = 0;
      dSum = 0;
      binIndex++;
    }
  }

  return hRebin;
}
//_____________________________________________________________________________
void ProcessUniFlow::AddTask(FlowTask* task)
{
  if(!task) return;

  if(task->fSpecies == FlowTask::kRefs) fvTasks.insert(fvTasks.begin(), task);
  else fvTasks.push_back(task);

  return;
}
//_____________________________________________________________________________
void ProcessUniFlow::SuggestPtBinning(TH3D* histEntries, TProfile3D* profFlowOrig, FlowTask* task, Short_t binMult)
{
  if(!histEntries) return;
  if(!profFlowOrig) return;
  if(!task) return;

  Short_t binMultLow = histEntries->GetXaxis()->FindFixBin(fdMultBins[binMult]);
  Short_t binMultHigh = histEntries->GetXaxis()->FindFixBin(fdMultBins[binMult+1]) - 1; // for rebin both bins are included (so that one needs to lower)

  TH1D* hPtProj = (TH1D*) histEntries->ProjectionY("hPtProj",binMultLow,binMultHigh);

  const Double_t dMinEntries = task->fSuggestPtBinEntries;
  std::vector<Double_t> vecBins;
  std::vector<Double_t> vecContents;

  printf("vector pre: %lu\n", vecBins.size());
  Double_t dCount = 0;
  const Short_t iNBins = hPtProj->GetNbinsX();
  Short_t iBin = 1;

  // if some of the pt bins are set, the suggestion started with first bin after the last edge
  if(task->fNumPtBins > 0 && task->fNumPtBins < iNBins)
  {
    iBin = hPtProj->FindFixBin(task->fPtBinsEdges[task->fNumPtBins]) + 1;
    printf("LastBin %d (%g) | iBin %d (%g)\n",task->fNumPtBins,task->fPtBinsEdges[task->fNumPtBins], iBin, hPtProj->GetBinLowEdge(iBin));
  }
  else { vecBins.push_back(hPtProj->GetXaxis()->GetXmin()); } // pushing first edge

  for(Short_t bin(iBin); bin < iNBins+1; bin++)
  {
    dCount += hPtProj->GetBinContent(bin);
    if(dCount > dMinEntries)
    {
      vecBins.push_back(hPtProj->GetXaxis()->GetBinUpEdge(bin));
      vecContents.push_back(dCount);
      dCount = 0;
    }
  }
  vecContents.push_back(dCount); // pushing last (remaining) bin content
  vecBins.push_back(hPtProj->GetXaxis()->GetXmax()); // pushing last edge (low edge of underflowbin)
  const Short_t iNumBinEdges = vecBins.size();
  printf("vector post: %hd\n", iNumBinEdges);

  // filling internal pt binning with suggestion result
  for(Short_t index(0); index < iNumBinEdges; index++)
  {
    task->fPtBinsEdges[index + task->fNumPtBins] = vecBins.at(index); // NOTE +1 ???
  }
  task->fNumPtBins += iNumBinEdges-1;


  // TODO hBinsContent reimplmente
  // TH1D* hBinsContent = new TH1D("hBinsContent","BinContent", task->fNumPtBins, task->fPtBinsEdges);
  // for(Short_t iBin(0); iBin < task->-1; iBin++)
  // {
  //   hBinsContent->SetBinContent(iBin+1,vecContents.at(iBin));
  // }

  TLine* line = new TLine();
  line->SetLineColor(kRed);

  TCanvas* cPtBins = new TCanvas("cPtBins","cPtBins",600,600);
  cPtBins->Divide(1,2);
  cPtBins->cd(1);
  hPtProj->Draw();
  cPtBins->cd(2);
  // hBinsContent->SetMinimum(0);
  // hBinsContent->Draw();
  // line->DrawLine(task->fPtBinsEdges[0],dMinEntries, task->fPtBinsEdges[task->fNumPtBins], dMinEntries);
  cPtBins->cd(1);

  Double_t dPt = 0;
  printf("Suggested binning: ");
  for(Short_t index(0); index < task->fNumPtBins+1; index++)
  {
    dPt = task->fPtBinsEdges[index];
    printf("%hd (%g) | ",index,dPt);
    line->DrawLine(dPt,0,dPt,1.05*hPtProj->GetMaximum());
  }
  printf("\n");

  cPtBins->SaveAs(Form("%s/suggestBins/SuggestPtBins_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,binMult,fsOutputFileFormat.Data()));
}
//_____________________________________________________________________________
void ProcessUniFlow::TestProjections()
{
  Info("Testing profile projections");

  ffInputFile->cd("UniFlow");
  // ffInputFile->ls();

  // TList* lFlow = (TList*) gDirectory->Get("UniFlow");
  // lFlow->ls();

  //
  // TList* lRef = (TList*) lFlow->FindObject("fFlowRefs");
  // TProfile* Ref = (TProfile*) lRef->FindObject("fpRefs_<2>_harm2_gap-10_sample0");
  // if(!Ref) { Error("NotFound"); return; }
  // // Ref->Draw();
  // // NOTE: reference flow works
  //
  // // projections charged flow
  // TList* lCharged = (TList*) lFlow->FindObject("fFlowCharged");
  // // lCharged->ls();
  // TProfile2D* p2Charged = (TProfile2D*) lCharged->FindObject("fp2Charged_<2>_harm2_gap-10_sample0");
  // if(!p2Charged) { Error("NotFound"); return; }
  // TProfile* p2ChargedProjY = p2Charged->ProfileY("p2ChargedProjY",10,10);

  // TCanvas* cCharged = new TCanvas("cCharged","cCharged",1000,1000);
  // cCharged->Divide(1,2);
  // cCharged->cd(1);
  // p2Charged->Draw("colz");
  // cCharged->cd(2);
  // p2ChargedProjY->Draw();

  // projections V0s
  // TList* lV0s = (TList*) gDirectory->Get("Flow_V0s_UniFlow");
  TList* lV0s = (TList*) gDirectory->Get("Flow_Phi_UniFlow");
  if(!lV0s) { Error("NotFound"); return; }
  lV0s->ls();

  // entries
  TH3D* h3K0sEntries = (TH3D*) lV0s->FindObject("fh3PhiEntriesSignal_gap-10")->Clone("h3K0sEntries");
  // TH3D* h3K0sEntries = (TH3D*) lV0s->FindObject("fh3V0sEntriesK0s_gap-10")->Clone("h3K0sEntries");
  // if(!h3K0sEntries) { Error("NotFound"); return; }
  // TH1D* h3K0sEntriesProjX = h3K0sEntries->ProjectionX("h3K0sEntriesProjX"); // whole projeciton
  // TH1D* h3K0sEntriesProjY = h3K0sEntries->ProjectionY("h3K0sEntriesProjY"); // whole projeciton
  TH1D* h3K0sEntriesProjZ = h3K0sEntries->ProjectionZ("h3K0sEntriesProjZ"); // whole projeciton
  // TH1D* h3K0sEntriesProjZsub = h3K0sEntries->ProjectionZ("h3K0sEntriesProjZsub",10,12,10,12); // projection in certain cent & pt range
  //
  // TCanvas* cK0s = new TCanvas("cK0s","cK0s",1000,1000);
  // cK0s->Divide(2,3);
  // cK0s->cd(1);
  // h3K0sEntries->Draw();
  // cK0s->cd(2);
  // h3K0sEntriesProjZ->Draw();
  // cK0s->cd(3);
  // h3K0sEntriesProjX->Draw();
  // cK0s->cd(4);
  // h3K0sEntriesProjY->Draw();
  // cK0s->cd(5);
  // h3K0sEntriesProjZsub->Draw();
  // cK0s->cd(6);
  // NOTE seems to work properly

  // correlations
  // TProfile3D* p3K0sCor = (TProfile3D*) lV0s->FindObject("fp3V0sCorrK0s_<2>_harm2_gap-10")->Clone("p3K0sCor");
  TProfile3D* p3K0sCor = (TProfile3D*) lV0s->FindObject("fp3PhiCorr_<2>_harm2_gap00")->Clone("p3K0sCor");
  if(!p3K0sCor) { Error("NotFound"); return; }

  p3K0sCor->GetXaxis()->SetRange(20,50);
  // p3K0sCor->GetZaxis()->SetRange(1,10);

  TProfile2D* profROOT = p3K0sCor->Project3DProfile("yz");
  TProfile2D* prof = Project3DProfile(p3K0sCor);

  TProfile* profROOTy = profROOT->ProfileX("profROOTy");
  TProfile* profy = prof->ProfileX("profy",40,70);
  // profy->Rebin(30);

  // TCanvas* cK0sCor = new TCanvas("cK0sCor","cK0sCor",1000,1000);
  // cK0sCor->Divide(2,3);
  // cK0sCor->cd(1);
  // p3K0sCor->Draw();
  // cK0sCor->cd(2);
  // h3K0sEntriesProjZ->Draw();
  // cK0sCor->cd(3);
  // profROOT->Draw("colz");
  // cK0sCor->cd(4);
  // prof->Draw("colz");
  // cK0sCor->cd(5);
  // profROOTy->Draw();
  // cK0sCor->cd(6);
  // profy->Draw();




  return;
}
//_____________________________________________________________________________
TProfile2D* ProcessUniFlow::Project3DProfile(const TProfile3D* prof3dorig)
{
  if(!prof3dorig) return 0x0;
  TProfile3D* prof3d = (TProfile3D*) prof3dorig->Clone();

  Int_t iBinFirst = prof3d->GetXaxis()->GetFirst();
  Int_t iBinLast = prof3d->GetXaxis()->GetLast();
  Int_t iNumBins = prof3d->GetNbinsX();
  Int_t iNumBinsAxis = prof3d->GetXaxis()->GetNbins();
  // printf("Bins:  %d - %d (%d | %d) \n", iBinFirst,iBinLast,iNumBins,iNumBinsAxis);

  // // making 3d hist from 3d profile
  // TH3D* hist3d = prof3d->ProjectionXYZ();   //NOTE do not care about range !!!
  // TH3D* hist3d_entry = prof3d->ProjectionXYZ("hist3d_entry","B");   //NOTE do not care about range !!!
  // TH3D* hist3d_weight = prof3d->ProjectionXYZ("hist3d_weight","W");   //NOTE do not care about range !!!

  TProfile2D* prof2d_test = DoProjectProfile2D(prof3d,"prof2d_test","",prof3d->GetYaxis(),prof3d->GetZaxis(),1,0,0);
  // prof2d_test->Draw("colz");

  return prof2d_test;

  // resulting profile
  // TProfile* result = new TProfile("result","result",100,0,10);
  // for(Int_t i(0); i < 10; i++) result->Fill(i,1);

  //
  // TCanvas* canTest = new TCanvas("canTest");
  // canTest->Divide(2,2);
  // canTest->cd(1);
  // prof3d->Draw("box");
  // canTest->cd(2);
  // hist3d->Draw("box");
  // canTest->cd(3);
  // hist3d_entry->Draw("box");
  // canTest->cd(4);
  // hist3d_weight->Draw("box");


  // return result;
}
//_____________________________________________________________________________
TProfile2D * ProcessUniFlow::DoProjectProfile2D(TProfile3D* h3, const char* name, const char * title, TAxis* projX, TAxis* projY,
                                           bool originalRange, bool useUF, bool useOF) const
{
// internal method to project to a 2D Profile
 // called from TH3::Project3DProfile but re-implemented in case of the TPRofile3D since what is done is different

 // projX, projY: axes of the orifinal histogram to which the projection is done (e.g. xy)

 // Get the ranges where we will work.
 Int_t ixmin = projX->GetFirst();
 Int_t ixmax = projX->GetLast();
 Int_t iymin = projY->GetFirst();
 Int_t iymax = projY->GetLast();
 if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
 if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
 Int_t nx = ixmax-ixmin+1;
 Int_t ny = iymax-iymin+1;

 // Create the projected profiles
 TProfile2D *p2 = 0;
 // Create always a new TProfile2D (not as in the case of TH3 projection)

 const TArrayD *xbins = projX->GetXbins();
 const TArrayD *ybins = projY->GetXbins();
 // assume all axis have variable bins or have fixed bins
 if ( originalRange ) {
    if (xbins->fN == 0 && ybins->fN == 0) {
       p2 = new TProfile2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                           ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
    } else {
       p2 = new TProfile2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
    }
 } else {
    if (xbins->fN == 0 && ybins->fN == 0) {
       p2 = new TProfile2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                           ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
    } else {
       p2 = new TProfile2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
    }
 }

 // new profile p2 is set according to axis ranges (keeping originals or not)

 // weights
 bool useWeights = (h3->fBinSumw2.fN != 0); //array elements
 if (useWeights) p2->Sumw2();

 // make projection in a 3D first // from 3D profile -> TH3
 TH3D * h3dW = h3->ProjectionXYZ("h3temp-W","W"); // getbincontent*getBinentries
 TH3D * h3dN = h3->ProjectionXYZ("h3temp-N","B"); // bin content is original profile = GetEntriesBin

 // fix ???
 h3dW->GetXaxis()->SetRange(h3->GetXaxis()->GetFirst(),h3->GetXaxis()->GetLast());
 h3dW->GetYaxis()->SetRange(h3->GetYaxis()->GetFirst(),h3->GetYaxis()->GetLast());
 h3dW->GetZaxis()->SetRange(h3->GetZaxis()->GetFirst(),h3->GetZaxis()->GetLast());
 h3dN->GetXaxis()->SetRange(h3->GetXaxis()->GetFirst(),h3->GetXaxis()->GetLast());
 h3dN->GetYaxis()->SetRange(h3->GetYaxis()->GetFirst(),h3->GetYaxis()->GetLast());
 h3dN->GetZaxis()->SetRange(h3->GetZaxis()->GetFirst(),h3->GetZaxis()->GetLast());




 h3dW->SetDirectory(0); h3dN->SetDirectory(0); // istograms does not bellow to any directorz ???

 // note that h3dW is always a weighted histogram - so we need to compute error in the projection
 TAxis * projX_hW = h3dW->GetXaxis();
 TAxis * projX_hN = h3dN->GetXaxis();
 if (projX == h3->GetYaxis() ) {  projX_hW =  h3dW->GetYaxis();  projX_hN =  h3dN->GetYaxis(); }
 if (projX == h3->GetZaxis() ) {  projX_hW =  h3dW->GetZaxis();  projX_hN =  h3dN->GetZaxis(); }
 TAxis * projY_hW = h3dW->GetYaxis();
 TAxis * projY_hN = h3dN->GetYaxis();
 if (projY == h3->GetXaxis() ) {  projY_hW =  h3dW->GetXaxis();  projY_hN =  h3dN->GetXaxis(); }
 if (projY == h3->GetZaxis() ) {  projY_hW =  h3dW->GetZaxis();  projY_hN =  h3dN->GetZaxis(); }
 // checking the axes

 // TH3 -> TH2
 TH2D * h2W = DoProject2D(h3dW,"htemp-W","",projX_hW, projY_hW, true, originalRange, useUF, useOF);
 TH2D * h2N = DoProject2D(h3dN,"htemp-N","",projX_hN, projY_hN, useWeights, originalRange, useUF, useOF);
 h2W->SetDirectory(0); h2N->SetDirectory(0);


 // fill the bin content
 R__ASSERT( h2W->fN == p2->fN );
 R__ASSERT( h2N->fN == p2->fN );
 R__ASSERT( h2W->GetSumw2()->fN != 0); // h2W should always be a weighted histogram since h3dW is weighted


 // filling the new tprofile2D
 for (int i = 0; i < p2->fN ; ++i) {
    //std::cout << " proj bin " << i << "  " <<  h2W->fArray[i] << "  " << h2N->fArray[i] << std::endl;
    p2->fArray[i] = h2W->fArray[i];   // array of profile is sum of all values
    p2->GetSumw2()->fArray[i]  = h2W->GetSumw2()->fArray[i];   // array of content square of profile is weight square of the W projected histogram
    p2->SetBinEntries(i, h2N->fArray[i] );
    if (useWeights) p2->GetBinSumw2()->fArray[i] = h2N->GetSumw2()->fArray[i];    // sum of weight squares are stored to compute errors in h1N histogram
 }
 // delete the created histograms
 delete h3dW;
 delete h3dN;
 delete h2W;
 delete h2N;

 // Also we need to set the entries since they have not been correctly calculated during the projection
 // we can only set them to the effective entries
 p2->SetEntries( p2->GetEffectiveEntries() );

 return p2;
}
//_____________________________________________________________________________
TH2D* ProcessUniFlow::DoProject2D(TH3D* h3, const char * name, const char * title, TAxis* projX, TAxis* projY,
                    bool computeErrors, bool originalRange,
                    bool useUF, bool useOF) const
{
  // internal method performing the projection to a 2D histogram
     // called from TH3::Project3D

     TH2D *h2 = 0;

     // Get range to use as well as bin limits
     Int_t ixmin = projX->GetFirst();
     Int_t ixmax = projX->GetLast();
     Int_t iymin = projY->GetFirst();
     Int_t iymax = projY->GetLast();
     if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
     if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
     Int_t nx = ixmax-ixmin+1;
     Int_t ny = iymax-iymin+1;

      const TArrayD *xbins = projX->GetXbins();
      const TArrayD *ybins = projY->GetXbins();
      if ( originalRange )
      {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
                          ,projX->GetNbins(),&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1]
                          ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
         } else {
            h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
         }
      } else {
         if (xbins->fN == 0 && ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else if (ybins->fN == 0) {
            h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
                          ,nx,&xbins->fArray[ixmin-1]);
         } else if (xbins->fN == 0) {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1]
                          ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
         } else {
            h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
         }
      }

    //  // Copy the axis attributes and the axis labels if needed.
    //  THashList* labels1 = 0;
    //  THashList* labels2 = 0;
    //  // "xy"
    //  h2->GetXaxis()->ImportAttributes(projY);
    //  h2->GetYaxis()->ImportAttributes(projX);
    //  labels1 = projY->GetLabels();
    //  labels2 = projX->GetLabels();
    //  if (labels1) {
    //     TIter iL(labels1);
    //     TObjString* lb;
    //     Int_t i = 1;
    //     while ((lb=(TObjString*)iL())) {
    //        h2->GetXaxis()->SetBinLabel(i,lb->String().Data());
    //        i++;
    //     }
    //  }
    //  if (labels2) {
    //     TIter iL(labels2);
    //     TObjString* lb;
    //     Int_t i = 1;
    //     while ((lb=(TObjString*)iL())) {
    //        h2->GetYaxis()->SetBinLabel(i,lb->String().Data());
    //        i++;
    //     }
    //  }
    //  h2->SetLineColor(this->GetLineColor());
    //  h2->SetFillColor(this->GetFillColor());
    //  h2->SetMarkerColor(this->GetMarkerColor());
    //  h2->SetMarkerStyle(this->GetMarkerStyle());

     // Activate errors
     if ( computeErrors) h2->Sumw2();

     // Set references to the axis, so that the bucle has no branches.
     TAxis* out = 0;
     if ( projX != h3->GetXaxis() && projY != h3->GetXaxis() ) {
        out = h3->GetXaxis();
     } else if ( projX != h3->GetYaxis() && projY != h3->GetYaxis() ) {
        out = h3->GetYaxis();
     } else {
        out = h3->GetZaxis();
     }

     Int_t *refX = 0, *refY = 0, *refZ = 0;
     Int_t ixbin, iybin, outbin;
     if ( projX == h3->GetXaxis() && projY == h3->GetYaxis() ) { refX = &ixbin;  refY = &iybin;  refZ = &outbin; }
     if ( projX == h3->GetYaxis() && projY == h3->GetXaxis() ) { refX = &iybin;  refY = &ixbin;  refZ = &outbin; }
     if ( projX == h3->GetXaxis() && projY == h3->GetZaxis() ) { refX = &ixbin;  refY = &outbin; refZ = &iybin;  }
     if ( projX == h3->GetZaxis() && projY == h3->GetXaxis() ) { refX = &iybin;  refY = &outbin; refZ = &ixbin;  }
     if ( projX == h3->GetYaxis() && projY == h3->GetZaxis() ) { refX = &outbin; refY = &ixbin;  refZ = &iybin;  }
     if ( projX == h3->GetZaxis() && projY == h3->GetYaxis() ) { refX = &outbin; refY = &iybin;  refZ = &ixbin;  }
     R__ASSERT (refX != 0 && refY != 0 && refZ != 0);

     // Fill the projected histogram excluding underflow/overflows if considered in the option
     // if specified in the option (by default they considered)
     Double_t totcont  = 0;

     Int_t outmin = out->GetFirst();
     Int_t outmax = out->GetLast();
     // GetFirst(), GetLast() can return (0,0) when the range bit is set artifically (see TAxis::SetRange)
     if (outmin == 0 && outmax == 0) { outmin = 1; outmax = out->GetNbins(); }
     // correct for underflow/overflows
     if (useUF && !out->TestBit(TAxis::kAxisRange) )  outmin -= 1;
     if (useOF && !out->TestBit(TAxis::kAxisRange) )  outmax += 1;

     for (ixbin=0;ixbin<=1+projX->GetNbins();ixbin++){
        if ( projX->TestBit(TAxis::kAxisRange) && ( ixbin < ixmin || ixbin > ixmax )) continue;
        Int_t ix = h2->GetYaxis()->FindBin( projX->GetBinCenter(ixbin) );

        for (iybin=0;iybin<=1+projY->GetNbins();iybin++){
           if ( projY->TestBit(TAxis::kAxisRange) && ( iybin < iymin || iybin > iymax )) continue;
           Int_t iy = h2->GetXaxis()->FindBin( projY->GetBinCenter(iybin) );

           Double_t cont = 0;
           Double_t err2 = 0;

           // loop on the bins to be integrated (outbin should be called inbin)
           for (outbin = outmin; outbin <= outmax; outbin++){

              Int_t bin = h3->GetBin(*refX,*refY,*refZ);

              // sum the bin contents and errors if needed
              cont += h3->GetBinContent(bin);
              if (computeErrors) {
                 Double_t exyz = h3->GetBinError(bin);
                 err2 += exyz*exyz;
              }

           }

           // remember axis are inverted
           h2->SetBinContent(iy , ix, cont);
           if (computeErrors) h2->SetBinError(iy, ix, TMath::Sqrt(err2) );
           // sum all content
           totcont += cont;

        }
     }

     // since we use fill we need to reset and recalculate the statistics (see comment in DoProject1D )
     // or keep original statistics if consistent sumw2
     bool resetStats = true;
     double eps = 1.E-12;
     if (h3->IsA() == TH3F::Class() ) eps = 1.E-6;
     if (h3->fTsumw != 0 && TMath::Abs( h3->fTsumw - totcont) <  TMath::Abs(h3->fTsumw) * eps) resetStats = false;

     bool resetEntries = resetStats;
     // entries are calculated using underflow/overflow. If excluded entries must be reset
     resetEntries |= !useUF || !useOF;

     if (!resetStats) {
        Double_t stats[TH1::kNstat];
        Double_t oldst[TH1::kNstat]; // old statistics
        for (Int_t i = 0; i < TH1::kNstat; ++i) { oldst[i] = 0; }
        h3->GetStats(oldst);
        std::copy(oldst,oldst+TH1::kNstat,stats);
        // not that projX refer to Y axis and projX refer to the X axis of projected histogram
        // nothing to do for projection in Y vs X
        if ( projY == h3->GetXaxis() && projX == h3->GetZaxis() ) {  // case XZ
           stats[4] = oldst[7];
           stats[5] = oldst[8];
           stats[6] = oldst[9];
        }
        if ( projY == h3->GetYaxis() ) {
           stats[2] = oldst[4];
           stats[3] = oldst[5];
           if ( projX == h3->GetXaxis() )  { // case YX
              stats[4] = oldst[2];
              stats[5] = oldst[3];
           }
           if ( projX == h3->GetZaxis() )  { // case YZ
              stats[4] = oldst[7];
              stats[5] = oldst[8];
              stats[6] = oldst[10];
           }
        }
        else if  ( projY == h3->GetZaxis() ) {
           stats[2] = oldst[7];
           stats[3] = oldst[8];
           if ( projX == h3->GetXaxis() )  { // case ZX
              stats[4] = oldst[2];
              stats[5] = oldst[3];
              stats[6] = oldst[9];
           }
           if ( projX == h3->GetYaxis() )  { // case ZY
              stats[4] = oldst[4];
              stats[5] = oldst[5];
              stats[6] = oldst[10];
           }
        }
        // set the new statistics
        h2->PutStats(stats);
     }
     else {
        // recalculate the statistics
        h2->ResetStats();
     }

     if (resetEntries) {
        // use the effective entries for the entries
        // since this  is the only way to estimate them
        Double_t entries =  h2->GetEffectiveEntries();
        if (!computeErrors) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
        h2->SetEntries( entries );
     }
     else {
        h2->SetEntries( h3->fEntries );
     }


     return h2;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowPhi(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass)
{
  //gStyle->SetOptFit(1110);

  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowPhi"); return kFALSE; }
  if(!hInvMass) { Error("InvMass (signal) histogram does not exists!","ExtractFlowPhi"); return kFALSE; }
  if(!hInvMassBG) { Error("InvMass (BG) histogram does not exists!","ExtractFlowPhi"); return kFALSE; }
	if(!hFlowMass) { Error("Flow-mass histogram does not exists!","ExtractFlowPhi"); return kFALSE;	}
	if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowPhi"); return kFALSE; }

  TString sFitOpt = TString("R");
  if(!fbDebug) sFitOpt.Append("Q");

  hFlowMass->SetMinimum(-0.1);
  // if(hFlowMass->GetMaximum() < 0.5) hFlowMass->SetMaximum(0.5);
  hFlowMass->SetMaximum(0.5);

	// Reseting the canvas (removing drawn things)
	canFitInvMass->Clear();

  TLatex* latex = new TLatex();
  latex->SetNDC();
  // latex->SetLineColor(kRed);

	// Fitting Phi
	const TString sOutputFormat = fsOutputFileFormat;
  // Overall fitting region
  Double_t dFittingLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dFittingHigh = hInvMass->GetXaxis()->GetXmax();

  // InvMass rejection region for sidebands fitting (based on mean ± iNumSigmas*sigma)
  Double_t dMassLow = 0.;
  Double_t dMassHigh = 0.;

	UShort_t iNumSigmas = 3; if(task->fFlowFitRejectNumSigmas > 0) iNumSigmas = task->fFlowFitRejectNumSigmas;
	Double_t dMeanShot = 0;
	Double_t dSigmaShot = 0;

  // subtraction of BG
  // normalisation BG to Signal in region of interest
  const Double_t dNormMassLow = 1.03;
  const Double_t dNormMassHigh = 1.06;

  Short_t iBinNormLow = hInvMassBG->FindFixBin(dNormMassLow);
  Short_t iBinNormHigh = hInvMassBG->FindFixBin(dNormMassHigh);
  printf("Bin low %g | high %g\n",hInvMassBG->GetBinCenter(iBinNormLow),hInvMassBG->GetBinCenter(iBinNormHigh));

  Double_t dIntSignal = 0;
  Double_t dIntBG = 0;
  for(Short_t iBin(iBinNormLow); iBin < iBinNormHigh+1; iBin++)
  {
    dIntSignal += hInvMass->GetBinContent(iBin);
    dIntBG += hInvMassBG->GetBinContent(iBin);
  }
  Double_t dScaleFact = dIntSignal/dIntBG;
  // printf("integral: sig %g | bg %g | scale %g\n",dIntSignal,dIntBG,dScaleFact);

  // TODO which one?
  TH1D* hInvMassBG_scaled = (TH1D*) hInvMassBG->Clone("hInvMassBG_scaled");
  // hInvMass->Scale(1/dIntSignal);
  // hInvMassBG_scaled->Scale(dScaleFact);
  // hInvMassBG_scaled->Scale(1/dScaleFact);
  hInvMassBG_scaled->SetLineColor(kRed);

  TH1D* hInvMass_subs = (TH1D*) hInvMass->Clone("hInvMass_subs");
  if(task->fFlowPhiSubtLS) { hInvMass_subs->Add(hInvMassBG_scaled,-1); }


  // TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1200,1200);
	canFitInvMass->Divide(3,2);
  // canFitInvMass->cd(1);
  // hInvMass->Draw();
  // if(task->fFlowPhiSubtLS)
  // {
  //   hInvMassBG->Draw("same");
  //   hInvMassBG_scaled->Draw("same");
  // }
  // // canFitInvMass->cd(8);
  // // hInvMass_subs->Draw();

	TH1D* hInvMass_shot = 0x0;
	TH1D* hInvMass_side = 0x0;
	TH1D* hInvMass_residual = 0x0;
	TH1D* hInvMass_ratio = 0x0;

	TH1D* hFlowMass_side = 0x0;

	const Short_t iNumBinsMassFlow = hFlowMass->GetNbinsX();

	// inv mass fitts
	TF1* fitShot = 0x0;
	TF1* fitSide = 0x0;
	TF1* fitRatio = 0x0;

	// flow mass fits
	TF1* fitFlowSide = 0x0;
	TF1* fitFlowTot = 0x0;

	hInvMass_shot = (TH1D*) hInvMass_subs->Clone("hInvMass_shot");
  hInvMass_shot->SetMaximum(hInvMass_shot->GetMaximum()*2);
	hInvMass_side = (TH1D*) hInvMass_subs->Clone("hInvMass_side");
	hInvMass_residual = (TH1D*) hInvMass_subs->Clone("hInvMass_residual");
  hInvMass_residual->GetYaxis()->SetTitle("Residuals");

  // fitting first shot
  // Short_t iBinPeakLow = hInvMass_shot->FindFixBin(dInvMassPeakLow);
  // Short_t iBinPeakHigh = hInvMass_shot->FindFixBin(dInvMassPeakHigh);
  // for(Short_t bin(iBinPeakLow); bin < iBinPeakHigh; bin++)
  // {
  //   // hInvMass_shot->SetBinError(bin,10*hInvMass_subs->GetMaximum());
  // }
  // hInvMass_shot->SetMaximum(hInvMass_subs->GetMaximum());
  // hInvMass_shot->SetMinimum(hInvMass_subs->GetMinimum());

  Debug("====== Fitting InvMass first shot ========","ExtractFlowPhi");
  canFitInvMass->cd(1);
	fitShot = new TF1("fitShot","[0]*TMath::BreitWigner(x,[1],[2])+pol2(3)",dFittingLow,dFittingHigh);
	fitShot->SetNpx(10000);
  fitShot->SetParameter(0,hInvMass_shot->GetMaximum());
  fitShot->SetParameter(1,1.019445);
  fitShot->SetParLimits(1,1.018,1.020);
  fitShot->SetParameter(2,0.005);
  fitShot->SetParLimits(2,0.004,0.006);
	// fitShot->SetParameter(1,0.5);
	// fitShot->SetParLimits(1,0.485,0.515);
	// fitShot->SetParameter(2,0.01);
	// fitShot->SetParLimits(2,0.,0.05);
  // fitShot->SetParameter(3,hInvMass_shot->GetMaximum());
	hInvMass_shot->Fit("fitShot",sFitOpt.Data());

  latex->DrawLatex(0.18,0.81,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitShot->GetChisquare(), fitShot->GetNDF(),fitShot->GetChisquare()/fitShot->GetNDF()));
  latex->DrawLatex(0.18,0.73,Form("#color[8]{M = %.6g#pm%.2g}",fitShot->GetParameter(1),fitShot->GetParError(1)));
  latex->DrawLatex(0.18,0.65,Form("#color[8]{#Gamma = %0.5g#pm%.2g}",fitShot->GetParameter(2),fitShot->GetParError(2)));



	// TODO checking the fitting results

  // extract mean & sigma for sidebands fitting reagion
  dMeanShot = fitShot->GetParameter(4);
  dSigmaShot = fitShot->GetParameter(5);
  // dMassLow = dMeanShot - iNumSigmas*dSigmaShot;
  // dMassHigh = dMeanShot + iNumSigmas*dSigmaShot;
  dMassLow = 1.005;
  // Short_t iBinPeakLow = hInvMass_shot->FindFixBin(dInvMassPeakLow);
  dMassHigh = 1.038;

	Debug(Form("  Fitting region: %f - %f \n", dMassLow,dMassHigh), "ExtractFlowPhi");

  // Excluding mass peak window (setting errors to inf)
	const Short_t iNumBinsMass = hInvMass->GetNbinsX();
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		if(hInvMass_side->GetBinCenter(iMass) > dMassLow && hInvMass_side->GetBinCenter(iMass) < dMassHigh) { hInvMass_side->SetBinError(iMass,9999999999999); }
	}
	hInvMass_side->SetMaximum(0.55*hInvMass_shot->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_side->SetMinimum(0);
	// hInvMass_side->SetMinimum(hInvMass_shot->GetMinimum());


  // checking User defined fitting range
  if(task->fFlowFitRangeLow >= dFittingLow && task->fFlowFitRangeLow < dMassLow) dFittingLow = task->fFlowFitRangeLow;
  if(task->fFlowFitRangeHigh <= dFittingHigh && task->fFlowFitRangeHigh > dMassHigh) dFittingHigh = task->fFlowFitRangeHigh;

  Debug("====== Fitting InvMass sidebands ========","ExtractFlowPhi");
  canFitInvMass->cd(2);
  hInvMass_side->Draw();
	// fitSide = new TF1("fitSide","[0]*TMath::BreitWigner(x,[1],[2])+pol2(3)",dFittingLow,dFittingHigh);
	fitSide = new TF1("fitSide","pol2(0)",dFittingLow,dFittingHigh);
	fitSide->SetNpx(10000);
  // fitSide->SetParameter(0,hInvMass_side->GetMaximum());
  // fitSide->SetParameter(1,fitShot->GetParameter(1));
  // fitSide->SetParLimits(1,1.018,1.020);
  // fitSide->SetParameter(2,fitShot->GetParameter(2));
  // fitSide->SetParLimits(2,0.004,0.006);
  // fitSide->SetParameter(3,fitShot->GetParameter(3));
  // fitSide->SetParameter(4,fitShot->GetParameter(4));
  // fitSide->SetParameter(5,fitShot->GetParameter(5));

  // fitSide->SetParameter(0,hInvMass_side->GetMaximum());
  fitSide->SetParameter(1,fitShot->GetParameter(3));
  // fitSide->SetParLimits(1,1.018,1.020);
  fitSide->SetParameter(2,fitShot->GetParameter(4));
  // fitSide->SetParLimits(2,0.004,0.006);
  fitSide->SetParameter(3,fitShot->GetParameter(5));
  // fitSide->SetParameter(4,fitShot->GetParameter(4));
  // fitSide->SetParameter(5,fitShot->GetParameter(5));
	hInvMass_side->Fit("fitSide",sFitOpt.Data());
  latex->DrawLatex(0.38,0.15,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitSide->GetChisquare(), fitSide->GetNDF(),fitSide->GetChisquare()/fitSide->GetNDF()));


  TF1* fitResBG = new TF1("fitResBG","pol2(0)",dFittingLow,dFittingHigh);
  // fitResBG->SetParameter(0,fitShot->GetParameter(3));
  // fitResBG->SetParameter(1,fitShot->GetParameter(4));
  // fitResBG->SetParameter(2,fitShot->GetParameter(5));
  fitResBG->SetParameter(0,fitSide->GetParameter(0));
  fitResBG->SetParameter(1,fitSide->GetParameter(1));
  fitResBG->SetParameter(2,fitSide->GetParameter(2));
  fitResBG->SetLineColor(kBlue);
  fitResBG->SetLineStyle(kDashed);
  fitResBG->SetLineWidth(2);
  // fitResBG->Draw("same");


  Double_t dContent = 0;
	for(Short_t iMass(1); iMass < hInvMass->GetNbinsX()+1; iMass++)
	{
		dContent = hInvMass_side->GetBinContent(iMass) - fitResBG->Eval(hInvMass_side->GetBinCenter(iMass));
		hInvMass_residual->SetBinContent(iMass,dContent);
	}
  hInvMass_ratio = (TH1D*) hInvMass_residual->Clone("hInvMass_ratio");

	canFitInvMass->cd(3);
  hInvMass_residual->SetMaximum(0.1*hInvMass_residual->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_residual->Draw();

  TLine* vert = new TLine();
  vert->SetLineColor(kRed);
  vert->SetLineStyle(7);
  vert->DrawLine(task->fFlowFitRangeLow,hInvMass_residual->GetMinimum(),task->fFlowFitRangeLow,hInvMass_residual->GetMaximum());
  vert->DrawLine(task->fFlowFitRangeHigh,hInvMass_residual->GetMinimum(),task->fFlowFitRangeHigh,hInvMass_residual->GetMaximum());

	hInvMass_ratio->Sumw2();
	hInvMass_ratio->Divide(hInvMass);
  hInvMass_ratio->GetYaxis()->SetTitle("signal/total ratio");

  Debug("====== Fitting InvMass sig/tot ratio ========","ExtractFlowPhi");
	canFitInvMass->cd(4);

	//hInvMass_ratio->Draw();
	// fitRatio = new TF1("fitRatio","pol2(0)+[3]*TMath::BreitWigner(x,[4],[5])",0.99,1.07);
	// fitRatio->SetNpx(1000);
  // fitRatio->SetParameter(0,fitSide->GetParameter(0));
  // fitRatio->SetParameter(1,fitSide->GetParameter(1));
  // fitRatio->SetParameter(2,fitSide->GetParameter(2));
  // fitRatio->SetParameter(3,1);
  // fitRatio->SetParameter(4,fitSide->GetParameter(4));
  // fitRatio->SetParLimits(4,1.01,1.03);
  // fitRatio->SetParameter(5,fitSide->GetParameter(5));
  // fitRatio->SetParLimits(5,0.004,0.008);
	// hInvMass_ratio->Fit("fitRatio","R");

	fitRatio = new TF1("fitRatio","[0]*TMath::BreitWigner(x,[1],[2])+pol2(3)",dFittingLow,dFittingHigh);
	// fitRatio = new TF1("fitRatio","[0]*TMath::BreitWigner(x,[1],[2])+pol3(3)",dFittingLow,dFittingHigh);
	// fitRatio = new TF1("fitRatio","gaus(0)+pol3(3)",dFittingLow,dFittingHigh);
	// fitRatio = new TF1("fitRatio","gaus(0)+pol2(3)",dFittingLow,dFittingHigh);
	fitRatio->SetNpx(1000);
  fitRatio->SetParameter(0,0.8);
  // fitRatio->SetParLimits(0,0.5,1.);
  fitRatio->SetParameter(1,fitShot->GetParameter(1));
  fitRatio->SetParLimits(1,1.014,1.025);
  fitRatio->SetParameter(2,fitShot->GetParameter(2));
  fitRatio->SetParLimits(2,0.003,0.03);
  // fitRatio->SetParameter(3,10);
  // fitRatio->SetParLimits(3,-300,300);
  // fitRatio->SetParameter(4,10);
  // fitRatio->SetParLimits(4,-300,300);
  // fitRatio->SetParameter(5,10);
  // fitRatio->SetParLimits(5,-300,300);
  // fitRatio->SetParameter(6,10);
  // fitRatio->SetParLimits(6,-300,300);
	hInvMass_ratio->Fit("fitRatio",sFitOpt.Data());

  latex->DrawLatex(0.38,0.15,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitRatio->GetChisquare(), fitRatio->GetNDF(),fitRatio->GetChisquare()/fitRatio->GetNDF()));


  // return kFALSE;

	// flow mass Fitting
	hFlowMass_side = (TH1D*) hFlowMass->Clone("hFlowMass_side");
	hFlowMass_side->SetMaximum(1.2*hFlowMass->GetMaximum());
	hFlowMass_side->SetMinimum(0.8*hFlowMass->GetMinimum());
  hFlowMass_side->GetXaxis()->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  hFlowMass_side->GetYaxis()->SetTitle("v_{2}^{tot}");
  hFlowMass->GetXaxis()->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  hFlowMass->GetYaxis()->SetTitle("v_{2}^{tot}");

	// fitting side bands
  // Excluding mass peak window (setting errors to inf)
  for(Short_t iMass(1); iMass < iNumBinsMassFlow+1; iMass++)
  {
    if(hFlowMass_side->GetBinCenter(iMass) > dMassLow && hFlowMass_side->GetBinCenter(iMass) < dMassHigh)
      hFlowMass_side->SetBinError(iMass,9999999999999);
	}

  Debug("====== Fitting FlowMass sidebands ========","ExtractFlowPhi");
	canFitInvMass->cd(5);

  if(task->fFlowFitFixTerms) { fitFlowSide = new TF1("fitFlowSide","pol1(0)",dFittingLow,dFittingHigh); }
  else { fitFlowSide = new TF1("fitFlowSide","pol1(0)",hFlowMass_side->GetXaxis()->GetXmin(),hFlowMass_side->GetXaxis()->GetXmax()); }
	// fitFlowSide = new TF1("fitFlowSide","pol3(0)",dFittingLow,dFittingHigh);
	hFlowMass_side->Fit("fitFlowSide",sFitOpt.Data());

  Debug("====== Fitting FlowMass total ========","ExtractFlowPhi");
	canFitInvMass->cd(6);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol2(8)",0.99,1.07);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol2(4)) + ( 1-(gaus(1)+pol2(4)) )*pol1(7)",dFittingLow,dFittingHigh);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol1(8)",dFittingLow,dFittingHigh);
	fitFlowTot = new TF1("fitFlowTot","[0]*([1]*TMath::BreitWigner(x,[2],[3])+pol2(4)) + ( 1-([1]*TMath::BreitWigner(x,[2],[3])+pol2(4)) )*pol1(7)",dFittingLow,dFittingHigh);
	// fitFlowTot = new TF1("fitFlowTot","[0]*([1]*TMath::BreitWigner(x,[2],[3])+pol3(4)) + ( 1-([1]*TMath::BreitWigner(x,[2],[3])+pol3(4)) )*pol1(8)",dFittingLow,dFittingHigh);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol2(4)) + ( 1-(gaus(1)+pol2(4)) )*pol3(7)",dFittingLow,dFittingHigh);
	// Inv mass ratio signal/total
	// fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
	// fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
	// fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
	// fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
	// fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
	// fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
	// fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
	// // FlowMass backround / sidebands
	// fitFlowTot->FixParameter(8,fitFlowSide->GetParameter(0));
	// fitFlowTot->FixParameter(9,fitFlowSide->GetParameter(1));
	// fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));

  if(task->fFlowFitFixTerms) // fix terms from background flowe stiamtion
  {
    // Inv mass ratio signal/total

    fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
    fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
    fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
    fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
    fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
    fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
    // FlowMass backround / sidebands
    // fitFlowTot->SetParameter(7,fitFlowSide->GetParameter(0));
    fitFlowTot->FixParameter(7,fitFlowSide->GetParameter(0));
    fitFlowTot->FixParameter(8,fitFlowSide->GetParameter(1));
  	// fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(3));
  }
  else
  {
    // Inv mass ratio signal/total
    fitFlowTot->SetParameter(0,dFlow);
    fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
    fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
    fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
    fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
    fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
    fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
    // fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
    // FlowMass backround / sidebands
    // fitFlowTot->SetParameter(7,fitFlowSide->GetParameter(0));
    fitFlowTot->SetParameter(7,fitFlowSide->GetParameter(0));
    fitFlowTot->SetParameter(8,fitFlowSide->GetParameter(1));
    // fitFlowTot->SetParameter(8,fitFlowSide->GetParameter(0));
    // fitFlowTot->SetParameter(9,fitFlowSide->GetParameter(1));
    // fitFlowTot->SetParameter(9,fitFlowSide->GetParameter(1));
    // fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(3));
  }

  fitFlowTot->SetLineColor(8);
	hFlowMass->Fit("fitFlowTot","RIB");

	dFlow = fitFlowTot->GetParameter(0);
	dFlowError = fitFlowTot->GetParError(0);

  canFitInvMass->cd(6);

  latex->SetNDC();
  latex->DrawLatex(0.2,0.83,Form("#color[8]{v_{2} = %.3g#pm%.2g}",dFlow,dFlowError));
  latex->DrawLatex(0.2,0.75,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitFlowTot->GetChisquare(), fitFlowTot->GetNDF(),fitFlowTot->GetChisquare()/fitFlowTot->GetNDF()));


  fitShot = new TF1("fitShot","[0]*TMath::BreitWigner(x,[1],[2])+pol2(3)",dFittingLow,dFittingHigh);
	fitShot->SetNpx(10000);
	// fitShot->SetParameter(1,0.5);
	// fitShot->SetParLimits(1,0.485,0.515);
	// fitShot->SetParameter(2,0.01);
	// fitShot->SetParLimits(2,0.,0.05);
  // fitShot->SetParameter(3,hInvMass_shot->GetMaximum());
	hInvMass_shot->Fit("fitShot",sFitOpt.Data());



TLegend* leg = new TLegend(0.55,0.3,0.88,0.61);
leg->SetBorderSize(0);
leg->SetFillColorAlpha(0,0);

TLegend* leg2 = new TLegend(0.55,0.3,0.88,0.43);
leg2->SetBorderSize(0);
leg2->SetFillColorAlpha(0,0);


  TF1* fitMass = new TF1("fitMass","[0]*TMath::BreitWigner(x,[1],[2])+pol2(3)",1,1.06);
	// fitMass->SetNpx(50000);
  fitMass->SetParameter(0,hInvMass_shot->GetMaximum());
  fitMass->SetParameter(1,1.019445);
  fitMass->SetParLimits(1,1.018,1.020);
  fitMass->SetParameter(2,0.005);
  fitMass->SetParLimits(2,0.004,0.006);
  fitMass->SetLineColor(kBlue);


  SetStyle();


  TLatex * text = new TLatex();
  text->SetTextFont(42);

  TLatex * text2 = new TLatex();
  text2->SetTextFont(42);
  // text2->SetTextSizePixels(22);



  // TCanvas* canNice = new TCanvas("canNice","canNice",500,1000);
  TCanvas* canNice = new TCanvas("canNice","canNice",600,800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.0);
  pad1->SetLeftMargin(0.15);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.3);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();

  pad1->cd();
  // canNice->Divide(1,2);
  // canNice->cd(1);
  TH1* h1 = pad1->DrawFrame(1.,0.02 ,1.06, 300000);
  h1->SetStats(0);
  // h1->SetTitleSize(0.2,"xy");
  h1->SetTitle(";#it{m}_{inv} (GeV/#it{c}^{2}); Counts");
  h1->GetXaxis()->SetLabelSize(0.);
  hInvMass_subs->SetStats(0);
  hInvMass_subs->SetLineColor(kBlack);
  // fitMass->SetLineStyle(kDashed);
  fitMass->SetNpx(300);
  hInvMass_subs->SetMarkerColor(kBlack);
  hInvMass_subs->Draw("same");
  hInvMass_subs->Fit("fitMass","R");
  fitMass->Draw("same");
    // TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
    // axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    // axis->SetLabelSize(15);
    // axis->Draw();


  TF1* fitSig = new TF1("fitMass","[0]*TMath::BreitWigner(x,[1],[2])",1,1.06);
  fitSig->SetParameter(0,fitMass->GetParameter(0));
  fitSig->SetParameter(1,fitMass->GetParameter(1));
  fitSig->SetParameter(2,fitMass->GetParameter(2));
  fitSig->SetParameter(0,fitMass->GetParameter(0));
  fitSig->SetLineColor(kRed);
  fitSig->SetNpx(300);

  TF1* fitBg = new TF1("fitMass","pol2(3)",1,1.06);
  fitBg->SetParameter(3,fitMass->GetParameter(3));
  fitBg->SetParameter(4,fitMass->GetParameter(4));
  fitBg->SetParameter(5,fitMass->GetParameter(5));
  fitBg->SetLineColor(kBlack);
  fitBg->SetLineStyle(kDashed);

  leg->AddEntry(hInvMass_subs,"data","pl");
  leg->AddEntry(fitMass,"combined fit","l");
  leg->AddEntry(fitSig,"#phi signal (Breit-Wigner)","l");
  leg->AddEntry(fitBg,"background (pol2)","l");
  // leg->AddEntry(hInvMass_subs,"K^{+}+K^{-} yield","pel");
  // leg->AddEntry(fitMass,"Breit-Wigner + pol2","l");
  // leg->AddEntry(fitSig,"#phi #rightarrow K^{+}+K^{-}","l");
  // leg->AddEntry(fitBg,"comb. bg.","l");




  fitSig->Draw("same");
  fitBg->Draw("same");
  leg->Draw();


  leg2->AddEntry(fitFlowTot,"#it{v}_{2}^{ tot} fit","l");

  DrawLogo(1, 0.2, 0.83);
  text->DrawLatexNDC(0.2,0.77,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  text2->DrawLatexNDC(0.2,0.7,Form("0-20%% (V0A)"));
  text2->DrawLatexNDC(0.2,0.63,Form("|#eta| < 0.8, 1.5 < #it{p}_{T} < 2 GeV/#it{c}"));
  text2->DrawLatexNDC(0.65,0.83,Form("#phi #rightarrow K^{+} + K^{-}"));

  pad2->cd();
  // gPad->SetBottomMargin(0.5);
  TH1* h2 = pad2->DrawFrame(1.,0.08 ,1.06,0.199 );
  hFlowMass->SetStats(0);
  hFlowMass->SetLineColor(kBlack);
  h2->SetTickLength(0.07,"X");
  h2->SetTitle(";#it{m}_{K^{+}K^{-}} (GeV/#it{c}^{2}); v_{2}^{tot} {2,|#Delta#eta|>0.8}");
  h2->SetTitleSize(0.09,"X");
  h2->SetTitleSize(0.09,"Y");
  h2->SetTitleOffset(0.7,"Y");
  // h2->SetLabelSize(0.045,"xyz");
  h2->SetLabelSize(0.09,"xy");
  // h2->SetLabelOffset(0.01,"x");
  // h2->SetLabelOffset(0.01,"y");
  // h2->SetLabelOffset(0.01,"x");
  // h2->SetTitleOffset(0.1,"XY");
  hFlowMass->Draw("same");
  leg2->Draw();
  canNice->SaveAs(Form("%s/NICE_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),"Phi",task->fHarmonics,10*task->fEtaGap,1,"pdf"),"pdf");
  canNice->SaveAs(Form("%s/NICE_%s_n%d2_gap%02.2g_cent%d.%s",fsOutputFilePath.Data(),"Phi",task->fHarmonics,10*task->fEtaGap,1,"eps"),"eps");
//

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowK0sAlex(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowK0sAlex"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowK0sAlex"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowK0sAlex"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowK0sAlex"); return kFALSE; }

  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

  // ####### Alex's code
  Info("#### Start of Alex's segment #######################","ExtractFlowK0sAlex");
  TH1D* hInvMassA = (TH1D*) hInvMass;
  TH1D* hVnK0APx = (TH1D*) hFlowMass;

  canFitInvMass->cd(1);

  hInvMassA->Sumw2();
  hInvMassA->GetXaxis()->SetTitle("M_{K^{0}} (GeV/c^{2})");
  hInvMassA->SetMarkerStyle(20);
  hInvMassA->SetMinimum(0);
  hInvMassA->DrawCopy();


  TF1* fitA = new TF1(Form("fitA"), "[0] + [1]*(x-0.4) + [2]*(x-0.4)*sqrt(x-0.4) + [3]*(x-0.4)*(x-0.4) + [4]*(x-0.4)*(x-0.4)*sqrt(x-0.4) + [5]*(x-0.4)*(x-0.4)*(x-0.4) + gaus(6)", 0.4, 0.6);
  fitA->SetParameters(hInvMassA->GetMaximum()/10., -hInvMassA->GetMaximum()/50., hInvMassA->GetMaximum()/30., hInvMassA->GetMaximum()/10., hInvMassA->GetMaximum()/20., hInvMassA->GetMaximum()/20., hInvMassA->GetMaximum(), 0.4976, 0.001);
  fitA->SetParLimits(6, 0, hInvMassA->GetMaximum()*100);
  //fitA->FixParameter(7, 0.4976);
  fitA->SetParLimits(7, 0.48, 0.52);
  fitA->SetParLimits(8, 0, 1);
  hInvMassA->Fit(fitA, "RNL");

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 10)){

      fitA->SetParameters(fitA->GetParameter(0)/nfitsA, fitA->GetParameter(1), fitA->GetParameter(2), fitA->GetParameter(3), fitA->GetParameter(4), fitA->GetParameter(5), fitA->GetParameter(6), 0.4976, fitA->GetParameter(8)*nfitsA);
      //fitA->FixParameter(7, 0.4976);
      fitA->SetParLimits(7, 0.48, 0.52);
      fitA->SetParLimits(8, 0, 1);
      hInvMassA->Fit(fitA, "RNL");

      statusA = gMinuit->fCstatu.Data();
      nfitsA++;
  }

  cout<<"Fit iteration: "<<nfitsA<<endl;

  fitA->DrawCopy("same");

  latex->DrawLatex(0.17,0.81,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitA->GetChisquare(), fitA->GetNDF(),fitA->GetChisquare()/fitA->GetNDF()));
  latex->DrawLatex(0.17,0.73,Form("#color[8]{#mu = %.6g#pm%.2g}",fitA->GetParameter(7),fitA->GetParError(7)));
  latex->DrawLatex(0.17,0.65,Form("#color[8]{#sigma = %0.5g#pm%.2g}",fitA->GetParameter(8),fitA->GetParError(8)));

  TF1* poliA = new TF1(Form("poliA"), "[0] + [1]*(x-0.4) + [2]*(x-0.4)*sqrt(x-0.4) + [3]*(x-0.4)*(x-0.4) + [4]*(x-0.4)*(x-0.4)*sqrt(x-0.4) + [5]*(x-0.4)*(x-0.4)*(x-0.4)", 0.4, 0.6);
  poliA->SetParameters(fitA->GetParameter(0), fitA->GetParameter(1), fitA->GetParameter(2), fitA->GetParameter(3), fitA->GetParameter(4), fitA->GetParameter(5));
  poliA->SetLineColor(2);
  poliA->SetLineStyle(2);
  poliA->DrawCopy("same");

  TF1* gA = new TF1(Form("gausA"), "gaus", 0.4, 0.6);
  gA->SetParameters(fitA->GetParameter(6), fitA->GetParameter(7), fitA->GetParameter(8));
  gA->SetLineColor(4);
  gA->SetLineStyle(2);
  gA->DrawCopy("same");

  cout<<fitA->GetParameter(6)<<"     "<<fitA->GetParameter(7)<<"    "<<fitA->GetParameter(8)<<endl;

  canFitInvMass->cd(2);

  hVnK0APx->GetXaxis()->SetTitle("M_{K^{0}} (GeV/c^{2})");
  hVnK0APx->DrawCopy();

  TF1* fvna = new TF1(Form("fvna"), "[0]*gaus(1)/(gaus(4) + [7] + [8]*(x-0.4) + [9]*(x-0.4)*sqrt(x-0.4) + [10]*(x-0.4)*(x-0.4) + [11]*(x-0.4)*(x-0.4)*sqrt(x-0.4) + [12]*(x-0.4)*(x-0.4)*(x-0.4)) + ([13] + [14]*(x-0.4) + [15]*(x-0.4)*sqrt(x-0.4) + [16]*(x-0.4)*(x-0.4) + [17]*(x-0.4)*(x-0.4)*sqrt(x-0.4) + [18]*(x-0.4)*(x-0.4)*(x-0.4)) / ( [19] + [20]*(x-0.4) + [21]*(x-0.4)*sqrt(x-0.4) + [22]*(x-0.4)*(x-0.4) + [23]*(x-0.4)*(x-0.4)*sqrt(x-0.4) + [24]*(x-0.4)*(x-0.4)*(x-0.4) + gaus(25))*([28]+[29]*x)", 0.4, 0.6);

 Double_t parama[30] = {0.1, fitA->GetParameter(6), fitA->GetParameter(7), fitA->GetParameter(8), fitA->GetParameter(6), fitA->GetParameter(7), fitA->GetParameter(8), fitA->GetParameter(0), fitA->GetParameter(1), fitA->GetParameter(2), fitA->GetParameter(3), fitA->GetParameter(4), fitA->GetParameter(5), fitA->GetParameter(0), fitA->GetParameter(1), fitA->GetParameter(2), fitA->GetParameter(3), fitA->GetParameter(4), fitA->GetParameter(5), fitA->GetParameter(0), fitA->GetParameter(1), fitA->GetParameter(2), fitA->GetParameter(3), fitA->GetParameter(4), fitA->GetParameter(5), fitA->GetParameter(6), fitA->GetParameter(7), fitA->GetParameter(8), 4.65731e-01, -5.44560e+00};


 fvna->SetParName(0, "v_{2}");
 //fvna->SetParLimits(0, 0.001, 0.4);
 fvna->SetLineColor(2);
 fvna->SetParameters(parama);

 fvna->FixParameter(1, fitA->GetParameter(6));
 fvna->FixParameter(2, fitA->GetParameter(7));
 fvna->FixParameter(3, fitA->GetParameter(8));

 fvna->FixParameter(4, fitA->GetParameter(6));
 fvna->FixParameter(5, fitA->GetParameter(7));
 fvna->FixParameter(6, fitA->GetParameter(8));

 fvna->FixParameter(7, fitA->GetParameter(0));
 fvna->FixParameter(8, fitA->GetParameter(1));
 fvna->FixParameter(9, fitA->GetParameter(2));
 fvna->FixParameter(10, fitA->GetParameter(3));
 fvna->FixParameter(11, fitA->GetParameter(4));
 fvna->FixParameter(12, fitA->GetParameter(5));

 fvna->FixParameter(13, fitA->GetParameter(0));
 fvna->FixParameter(14, fitA->GetParameter(1));
 fvna->FixParameter(15, fitA->GetParameter(2));
 fvna->FixParameter(16, fitA->GetParameter(3));
 fvna->FixParameter(17, fitA->GetParameter(4));
 fvna->FixParameter(18, fitA->GetParameter(5));

 fvna->FixParameter(19, fitA->GetParameter(0));
 fvna->FixParameter(20, fitA->GetParameter(1));
 fvna->FixParameter(21, fitA->GetParameter(2));
 fvna->FixParameter(22, fitA->GetParameter(3));
 fvna->FixParameter(23, fitA->GetParameter(4));
 fvna->FixParameter(24, fitA->GetParameter(5));

 fvna->FixParameter(25, fitA->GetParameter(6));
 fvna->FixParameter(26, fitA->GetParameter(7));
 fvna->FixParameter(27, fitA->GetParameter(8));


  hVnK0APx->Fit(fvna, "RN");
  fvna->DrawCopy("same");

  dFlow = fvna->GetParameter(0);
  dFlowError = fvna->GetParError(0);

  latex->DrawLatex(0.2,0.83,Form("#color[8]{v_{2} = %.3g#pm%.2g}",dFlow,dFlowError));
  latex->DrawLatex(0.2,0.75,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fvna->GetChisquare(), fvna->GetNDF(),fvna->GetChisquare()/fvna->GetNDF()));

  Info("#### End of Alex's segment #######################","ExtractFlowK0sAlex");
  // ####### END of Alex's code

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowK0s(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass)
{
  //gStyle->SetOptFit(1110);

	if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowK0s"); return kFALSE; }
	if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowK0s"); return kFALSE; }
	if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowK0s"); return kFALSE; }
	if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowK0s"); return kFALSE; }

  // setting fitting options
  TString sFitOpt = TString("R");
  if(!fbDebug) sFitOpt.Append("Q");

  hFlowMass->SetMinimum(-0.1);
  if(hFlowMass->GetMaximum() < 0.5) hFlowMass->SetMaximum(0.5);

	// Reseting the canvas (removing drawn things)
	canFitInvMass->Clear();

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

	// Fitting K0s
	const TString sOutputFormat = fsOutputFileFormat;
  // Overall fitting region
  Double_t dFittingLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dFittingHigh = hInvMass->GetXaxis()->GetXmax();

  // InvMass rejection region for sidebands fitting (based on mean ± iNumSigmas*sigma)
  Double_t dMassLow = 0.;
  Double_t dMassHigh = 0.;

	UShort_t iNumSigmas = 7; if(task->fFlowFitRejectNumSigmas > 0) iNumSigmas = task->fFlowFitRejectNumSigmas;
	Double_t dMeanShot = 0;
	Double_t dSigmaShot = 0;

	//TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1200,1200);
	canFitInvMass->Divide(3,2);

	TH1D* hInvMass_side = 0x0;
	TH1D* hInvMass_residual = 0x0;
	TH1D* hInvMass_ratio = 0x0;

	TH1D* hFlowMass_side = 0x0;

	const Short_t iNumBinsMassFlow = hFlowMass->GetNbinsX();

	// inv mass fitts
	TF1* fitShot = 0x0;
	TF1* fitSide = 0x0;
	TF1* fitRatio = 0x0;

	// flow mass fits
	TF1* fitFlowSide = 0x0;
	TF1* fitFlowTot = 0x0;


	hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side");
	hInvMass_residual = (TH1D*) hInvMass->Clone("hInvMass_residual");

	Debug("====== Fitting InvMass first shot ========","ExtractFlowK0s");
	canFitInvMass->cd(1);
	fitShot = new TF1("fitShot","gaus(0)+pol2(3)",dFittingLow,dFittingHigh);
	fitShot->SetNpx(10000);
	fitShot->SetParameter(1,0.5);
	fitShot->SetParLimits(1,0.485,0.515);
	fitShot->SetParameter(2,0.01);
	fitShot->SetParLimits(2,0.,0.05);
	hInvMass->Fit("fitShot",Form("%s N",sFitOpt.Data()));
	hInvMass->Draw();
  fitShot->Draw("same");


  latex->DrawLatex(0.17,0.81,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitShot->GetChisquare(), fitShot->GetNDF(),fitShot->GetChisquare()/fitShot->GetNDF()));
  latex->DrawLatex(0.17,0.73,Form("#color[8]{#mu = %.6g#pm%.2g}",fitShot->GetParameter(1),fitShot->GetParError(1)));
  latex->DrawLatex(0.17,0.65,Form("#color[8]{#sigma = %0.5g#pm%.2g}",fitShot->GetParameter(2),fitShot->GetParError(2)));



	// TODO checking the fitting results

  // extract mean & sigma for sidebands fitting reagion
  dMeanShot = fitShot->GetParameter(1);
  dSigmaShot = fitShot->GetParameter(2);
  dMassLow = dMeanShot - iNumSigmas*dSigmaShot;
  dMassHigh = dMeanShot + iNumSigmas*dSigmaShot;
	Debug(Form("  Fitting region: %f - %f \n", dMassLow,dMassHigh), "ExtractFlowK0s");

  // Excluding mass peak window (setting errors to inf)
	const Short_t iNumBinsMass = hInvMass->GetNbinsX();
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		if(hInvMass_side->GetBinCenter(iMass) > dMassLow && hInvMass_side->GetBinCenter(iMass) < dMassHigh) { hInvMass_side->SetBinError(iMass,9999999999999); }
	}
	canFitInvMass->cd(2);
	hInvMass_side->SetMaximum(0.2*hInvMass->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_side->SetMinimum(0);
	hInvMass_side->Draw();

  // checking User defined fitting range
  if(task->fFlowFitRangeLow >= dFittingLow && task->fFlowFitRangeLow < dMassLow) dFittingLow = task->fFlowFitRangeLow;
  if(task->fFlowFitRangeHigh <= dFittingHigh && task->fFlowFitRangeHigh > dMassHigh) dFittingHigh = task->fFlowFitRangeHigh;

	// fitting background in sidebands
  Debug("====== Fitting InvMass side bands ========","ExtractFlowK0s");
	fitSide = new TF1("fitSide","pol3(0)",dFittingLow,dFittingHigh);
	fitSide->SetNpx(10000);
	hInvMass_side->Fit("fitSide",sFitOpt.Data());

  // latex->DrawLatex(0.2,0.83,Form("#color[8]{v_{2} = %.3g#pm%.2g}",dFlow,dFlowError));
  latex->DrawLatex(0.2,0.75,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitSide->GetChisquare(), fitSide->GetNDF(),fitSide->GetChisquare()/fitSide->GetNDF()));


	Double_t dContent = 0;
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		dContent = hInvMass_side->GetBinContent(iMass) - fitSide->Eval(hInvMass_side->GetBinCenter(iMass));
		hInvMass_residual->SetBinContent(iMass,dContent);
	}

  hInvMass_ratio = (TH1D*) hInvMass_residual->Clone("hInvMass_ratio");
  hInvMass_ratio->GetYaxis()->SetTitle("signal/total ratio");

	canFitInvMass->cd(3);
  hInvMass_residual->SetMaximum(0.01*hInvMass_residual->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_residual->Draw();
  hInvMass_residual->GetYaxis()->SetTitle("Residuals");


  TLine* vert = new TLine();
  vert->SetLineColor(kRed);
  vert->SetLineStyle(7);
  vert->DrawLine(task->fFlowFitRangeLow,hInvMass_residual->GetMinimum(),task->fFlowFitRangeLow,hInvMass_residual->GetMaximum());
  vert->DrawLine(task->fFlowFitRangeHigh,hInvMass_residual->GetMinimum(),task->fFlowFitRangeHigh,hInvMass_residual->GetMaximum());

  hInvMass_ratio->Sumw2();
  hInvMass_ratio->Divide(hInvMass);
	// hInvMass_ratio->Divide(hInvMass_ratio,hInvMass,1.,1.,"b");

  Debug("====== Fitting InvMass sig/tot ratio ========","ExtractFlowK0s");

	canFitInvMass->cd(4);

	//hInvMass_ratio->Draw();
	fitRatio = new TF1("fitRatio","gaus(0)+pol3(3)",dFittingLow,dFittingHigh);
	// fitRatio = new TF1("fitRatio","gaus(0)+pol2(3)",dFittingLow,dFittingHigh);
	// fitRatio = new TF1("fitRatio","gaus(0)+gaus(3)+pol2(6)",dFittingLow,dFittingHigh);
	fitRatio->SetNpx(1000);
	fitRatio->SetParameter(0,0.98);
	fitRatio->SetParLimits(0,0.1,1.);
	fitRatio->SetParameter(1,0.5);
	fitRatio->SetParLimits(1,0.48,0.51);
	fitRatio->SetParameter(2,0.01);
	fitRatio->SetParLimits(2,0.,0.05);
	// fitRatio->SetParameter(3,0.3);
	// fitRatio->SetParLimits(3,0.1,0.9);
	// fitRatio->SetParameter(4,0.5);
	// fitRatio->SetParLimits(4,0.48,0.51);
	// fitRatio->SetParameter(5,0.05);
	// fitRatio->SetParLimits(5,0.003,0.01);
	hInvMass_ratio->Fit("fitRatio",sFitOpt.Data());

  // latex->DrawLatex(0.2,0.83,Form("#color[8]{v_{2} = %.3g#pm%.2g}",dFlow,dFlowError));
  latex->DrawLatex(0.2,0.75,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitRatio->GetChisquare(), fitRatio->GetNDF(),fitRatio->GetChisquare()/fitRatio->GetNDF()));


	// flow mass Fitting
	hFlowMass_side = (TH1D*) hFlowMass->Clone("hFlowMass_side");
	hFlowMass_side->SetMaximum(1.5*hFlowMass->GetMaximum());
	hFlowMass_side->SetMinimum(0.5*hFlowMass->GetMinimum());
  hFlowMass_side->GetXaxis()->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  hFlowMass_side->GetYaxis()->SetTitle("v_{2}^{tot}");
  hFlowMass->GetXaxis()->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  hFlowMass->GetYaxis()->SetTitle("v_{2}^{tot}");

	// fitting side bands
	for(Short_t iMass(1); iMass < iNumBinsMassFlow+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hFlowMass_side->GetBinCenter(iMass) > dMassLow && hFlowMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hFlowMass_side->SetBinError(iMass,9999999999999);
		}
	}

  Debug("====== Fitting FlowMass sidebands ========","ExtractFlowK0s");
	canFitInvMass->cd(5);

  if(task->fFlowFitFixTerms) { fitFlowSide = new TF1("fitFlowSide","pol1(0)",dFittingLow,dFittingHigh); }
  else { fitFlowSide = new TF1("fitFlowSide","pol1(0)",hFlowMass_side->GetXaxis()->GetXmin(),hFlowMass_side->GetXaxis()->GetXmax()); }

	// fitFlowSide = new TF1("fitFlowSide","pol1(0)",dFittingLow,dFittingHigh);
	hFlowMass_side->Fit("fitFlowSide",sFitOpt.Data());

  Debug("====== Fitting FlowMass total flow ========","ExtractFlowK0s");
  canFitInvMass->cd(6);

	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+gaus(4)+pol2(7)) + ( 1-(gaus(1)+gaus(4)+pol3(7)) )*pol1(10)",dFittingLow,dFittingHigh);
	fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol1(8)",dFittingLow,dFittingHigh);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol2(4)) + ( 1-(gaus(1)+pol2(4)) )*pol1(7)",dFittingLow,dFittingHigh);

  if(task->fFlowFitFixTerms) // fix terms from background flowe stiamtion
  {
    // Inv mass ratio signal/total
    fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
    fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
    fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
    fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
    fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
    fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
    fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
    // fitFlowTot->FixParameter(8,fitRatio->GetParameter(7));
    // fitFlowTot->FixParameter(9,fitRatio->GetParameter(8));
    // fitFlowTot->FixParameter(10,fitRatio->GetParameter(6));
    // FlowMass backround / sidebands
    fitFlowTot->FixParameter(8,fitFlowSide->GetParameter(0));
    fitFlowTot->FixParameter(9,fitFlowSide->GetParameter(1));
    // fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));
  }
  else
  {
    // Inv mass ratio signal/total
    fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
    fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
    fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
    fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
    fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
    fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
    fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
    // fitFlowTot->FixParameter(8,fitRatio->GetParameter(7));
    // fitFlowTot->FixParameter(9,fitRatio->GetParameter(8));
    // FlowMass backround / sidebands
    // fitFlowTot->SetParameter(8,fitFlowSide->GetParameter(0));
    // fitFlowTot->SetParameter(9,fitFlowSide->GetParameter(1));
    // fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));
  }
  fitFlowTot->SetParLimits(0,0.,1.);
	hFlowMass->Fit("fitFlowTot","R");

	dFlow = fitFlowTot->GetParameter(0);
	dFlowError = fitFlowTot->GetParError(0);

  latex->DrawLatex(0.2,0.83,Form("#color[8]{v_{2} = %.3g#pm%.2g}",dFlow,dFlowError));
  latex->DrawLatex(0.2,0.75,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitFlowTot->GetChisquare(), fitFlowTot->GetNDF(),fitFlowTot->GetChisquare()/fitFlowTot->GetNDF()));

	return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowLambda(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowLambda"); return kFALSE; }
	if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowLambda"); return kFALSE; }
	if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowLambda");	return kFALSE; }
	if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowLambda"); return kFALSE; }

	// Reseting the canvas (removing drawn things)
	canFitInvMass->Clear();

  TString sFitOpt = TString("RI");
  if(!fbDebug) sFitOpt.Append("Q");

  //gStyle->SetOptFit(1110);

    hFlowMass->SetMinimum(-0.1);
    if(hFlowMass->GetMaximum() < 0.5) hFlowMass->SetMaximum(0.5);


  // Fitting Lambda
  const TString sOutputFormat = fsOutputFileFormat;
  // Overall fitting region
  Double_t dFittingLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dFittingHigh = hInvMass->GetXaxis()->GetXmax();

  // InvMass rejection region for sidebands fitting (based on mean ± iNumSigmas*sigma)
  Double_t dMassLow = 0.;
  Double_t dMassHigh = 0.;

  UShort_t iNumSigmas = 6; if(task->fFlowFitRejectNumSigmas > 0) iNumSigmas = task->fFlowFitRejectNumSigmas;
  Double_t dMeanShot = 0;
  Double_t dSigmaShot = 0;

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

	//TCanvas* canFitInvMass = new TCanvas("canFitInvMass","FitInvMass",1200,1200);
	canFitInvMass->Divide(3,2);

	TH1D* hInvMass_side = 0x0;
	TH1D* hInvMass_residual = 0x0;
	TH1D* hInvMass_ratio = 0x0;

	TH1D* hFlowMass_side = 0x0;

	const Short_t iNumBinsMassFlow = hFlowMass->GetNbinsX();

	// inv mass fitts
	TF1* fitShot = 0x0;
	TF1* fitSide = 0x0;
	TF1* fitRatio = 0x0;

	// flow mass fits
	TF1* fitFlowSide = 0x0;
	TF1* fitFlowTot = 0x0;


	hInvMass_side = (TH1D*) hInvMass->Clone("hInvMass_side");
	hInvMass_residual = (TH1D*) hInvMass->Clone("hInvMass_residual");
	hInvMass_residual->GetYaxis()->SetTitle("Residuals");

  Debug("====== Fitting InvMass first shot ========","ExtractFlowLambda");
  canFitInvMass->cd(1);
	// fitShot = new TF1("fitShot","gaus(0)+pol2(3)",dFittingLow,dFittingHigh);
	fitShot = new TF1("fitShot","gaus(0)+pol2(3)",1.095,1.135);
	fitShot->SetNpx(10000);
	fitShot->SetParameter(1,1.115);
	fitShot->SetParLimits(1,1.113,1.12);
	fitShot->SetParameter(2,0.001);
  fitShot->SetParLimits(2,0.,0.0025);
  hInvMass->Fit("fitShot",Form("%s N",sFitOpt.Data()));
  hInvMass->Draw();
  fitShot->Draw("same");

  latex->DrawLatex(0.17,0.81,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitShot->GetChisquare(), fitShot->GetNDF(),fitShot->GetChisquare()/fitShot->GetNDF()));
  latex->DrawLatex(0.17,0.73,Form("#color[8]{#mu = %.6g#pm%.2g}",fitShot->GetParameter(1),fitShot->GetParError(1)));
  latex->DrawLatex(0.17,0.65,Form("#color[8]{#sigma = %0.5g#pm%.2g}",fitShot->GetParameter(2),fitShot->GetParError(2)));




	// TODO checking the fitting results

	// extract mean & sigma for sidebands fitting reagion
	dMeanShot = fitShot->GetParameter(1);
	dSigmaShot = fitShot->GetParameter(2);
	dMassLow = dMeanShot - iNumSigmas*dSigmaShot;
	dMassHigh = dMeanShot + iNumSigmas*dSigmaShot;
  Debug(Form(" Fitting region: %f - %f\n", dMassLow,dMassHigh),"ExtractFlowLambda");

  dMassLow = 1.105;
  dMassHigh = 1.131;

	const Short_t iNumBinsMass = hInvMass->GetNbinsX();
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hInvMass_side->GetBinCenter(iMass) > dMassLow && hInvMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hInvMass_side->SetBinError(iMass,9999999999999);
		}
	}
	canFitInvMass->cd(2);
	hInvMass_side->SetMaximum(0.55*hInvMass->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_side->SetMinimum(0);
	hInvMass_side->Draw();

  // checking User defined fitting range
  if(task->fFlowFitRangeLow >= dFittingLow && task->fFlowFitRangeLow < dMassLow) dFittingLow = task->fFlowFitRangeLow;
  if(task->fFlowFitRangeHigh <= dFittingHigh && task->fFlowFitRangeHigh > dMassHigh) dFittingHigh = task->fFlowFitRangeHigh;

  Debug("====== Fitting InvMass sidebands ========","ExtractFlowLambda");
  	fitSide = new TF1("fitSide","pol3(0)",dFittingLow,dFittingHigh);
	fitSide->SetNpx(10000);
	hInvMass_side->Fit("fitSide",sFitOpt.Data());
  latex->DrawLatex(0.17,0.81,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitSide->GetChisquare(), fitSide->GetNDF(),fitSide->GetChisquare()/fitSide->GetNDF()));

	Double_t dContent = 0;
	for(Short_t iMass(1); iMass < iNumBinsMass+1; iMass++)
	{
		dContent = hInvMass_side->GetBinContent(iMass) - fitSide->Eval(hInvMass_side->GetBinCenter(iMass));
		hInvMass_residual->SetBinContent(iMass,dContent);
	}

  hInvMass_ratio = (TH1D*) hInvMass_residual->Clone("hInvMass_ratio");
  hInvMass_ratio->GetYaxis()->SetTitle("signal/total ratio");


  canFitInvMass->cd(3);
  hInvMass_residual->SetMaximum(0.1*hInvMass_residual->GetMaximum()); // setting maximum & minimum (otherwise overshoteed with errors)
	hInvMass_residual->Draw();

  TLine* vert = new TLine();
  vert->SetLineColor(kRed);
  vert->SetLineStyle(7);
  vert->DrawLine(task->fFlowFitRangeLow,hInvMass_residual->GetMinimum(),task->fFlowFitRangeLow,hInvMass_residual->GetMaximum());
  vert->DrawLine(task->fFlowFitRangeHigh,hInvMass_residual->GetMinimum(),task->fFlowFitRangeHigh,hInvMass_residual->GetMaximum());


	hInvMass_ratio->Sumw2();
	hInvMass_ratio->Divide(hInvMass);
	hInvMass_ratio->SetMinimum(-0.5);
	hInvMass_ratio->SetMaximum(1.);

  Debug("====== Fitting InvMass sig/tot ratio ========","ExtractFlowLambda");
  canFitInvMass->cd(4);
	//hInvMass_ratio->Draw();
	fitRatio = new TF1("fitRatio","gaus(0)+pol2(3)",dFittingLow,dFittingHigh);
	// fitRatio = new TF1("fitRatio","gaus(0)+pol3(3)",dFittingLow,dFittingHigh);
	// fitRatio = new TF1("fitRatio","gaus(0)+gaus(3)+pol2(6)",dFittingLow,dFittingHigh);
	fitRatio->SetNpx(1000);
	fitRatio->SetParameter(0,hInvMass_ratio->GetMaximum());
  fitRatio->SetParLimits(0,0.7,1.);
	fitRatio->SetParameter(1,1.115);
	fitRatio->SetParameter(2,0.001);
	// fitRatio->SetParameter(3,0.98);
  // fitRatio->SetParLimits(3,0.,5.);
	// fitRatio->SetParameter(4,1.115);
	// fitRatio->SetParameter(5,0.005);
	// fitRatio->SetParLimits(5,0.01,0.3);
	hInvMass_ratio->Fit("fitRatio",sFitOpt.Data());
  latex->DrawLatex(0.17,0.81,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitRatio->GetChisquare(), fitRatio->GetNDF(),fitRatio->GetChisquare()/fitRatio->GetNDF()));

	// flow mass Fitting
	hFlowMass_side = (TH1D*) hFlowMass->Clone("hFlowMass_side");
	hFlowMass_side->SetMaximum(1.5*hFlowMass->GetMaximum());
	hFlowMass_side->SetMinimum(0.5*hFlowMass->GetMinimum());
  hFlowMass_side->GetXaxis()->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  hFlowMass_side->GetYaxis()->SetTitle("v_{2}^{tot}");
  hFlowMass->GetXaxis()->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
  hFlowMass->GetYaxis()->SetTitle("v_{2}^{tot}");

	// fitting side bands
	for(Short_t iMass(1); iMass < iNumBinsMassFlow+1; iMass++)
	{
		// Excluding mass peak window (setting errors to inf)
		if(hFlowMass_side->GetBinCenter(iMass) > dMassLow && hFlowMass_side->GetBinCenter(iMass) < dMassHigh)
		{
			hFlowMass_side->SetBinError(iMass,9999999999999);
		}
	}
  Debug("====== Fitting FlowMass sidebands ========","ExtractFlowLambda");
  canFitInvMass->cd(5);

  if(task->fFlowFitFixTerms) { fitFlowSide = new TF1("fitFlowSide","pol1(0)",dFittingLow,dFittingHigh); }
  else { fitFlowSide = new TF1("fitFlowSide","pol1(0)",hFlowMass_side->GetXaxis()->GetXmin(),hFlowMass_side->GetXaxis()->GetXmax()); }

  // fitFlowSide = new TF1("fitFlowSide","pol1(0)",dFittingLow,dFittingHigh);
	// fitFlowSide = new TF1("fitFlowSide","pol2(0)",dFittingLow,dFittingHigh);
	hFlowMass_side->Fit("fitFlowSide",sFitOpt.Data());

  Debug("====== Fitting FlowMass total flow ========","ExtractFlowLambda");
	canFitInvMass->cd(6);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+gaus(4)+pol2(7)) + ( 1-(gaus(1)+gaus(4)+pol2(7)) )*pol1(9)",dFittingLow,dFittingHigh);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+gaus(4)+pol3(7)) + ( 1-(gaus(1)+gaus(4)+pol3(7)) )*pol1(11)",dFittingLow,dFittingHigh);
	fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol2(4)) + ( 1-(gaus(1)+pol2(4)) )*pol1(7)",dFittingLow,dFittingHigh);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol1(8)",dFittingLow,dFittingHigh);
	// fitFlowTot = new TF1("fitFlowTot","[0]*(gaus(1)+pol3(4)) + ( 1-(gaus(1)+pol3(4)) )*pol2(8)",dFittingLow,dFittingHigh);

  if(task->fFlowFitFixTerms) // fix terms from background flowe stiamtion
  {
    // Inv mass ratio signal/total
    fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
    fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
    fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
    fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
    fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
    fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
    // fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
    // fitFlowTot->FixParameter(8,fitRatio->GetParameter(7));
    // fitFlowTot->FixParameter(9,fitRatio->GetParameter(8));
    // fitFlowTot->FixParameter(10,fitRatio->GetParameter(9));
    // FlowMass backround / sidebands
    fitFlowTot->FixParameter(7,fitFlowSide->GetParameter(0));
    fitFlowTot->FixParameter(8,fitFlowSide->GetParameter(1));
    // fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));
  }
  else
  {
    // Inv mass ratio signal/total
    fitFlowTot->FixParameter(1,fitRatio->GetParameter(0));
    fitFlowTot->FixParameter(2,fitRatio->GetParameter(1));
    fitFlowTot->FixParameter(3,fitRatio->GetParameter(2));
    fitFlowTot->FixParameter(4,fitRatio->GetParameter(3));
    fitFlowTot->FixParameter(5,fitRatio->GetParameter(4));
    fitFlowTot->FixParameter(6,fitRatio->GetParameter(5));
    // fitFlowTot->FixParameter(7,fitRatio->GetParameter(6));
    // fitFlowTot->FixParameter(8,fitRatio->GetParameter(7));
    // fitFlowTot->FixParameter(9,fitRatio->GetParameter(8));
    // FlowMass backround / sidebands
    fitFlowTot->SetParameter(0,dFlow);
    // fitFlowTot->SetParLimits(0,0.01,0.5);
    // if(dFlow > 0.01) fitFlowTot->SetParLimits(0,0.5*dFlow,2*dFlow);


    fitFlowTot->SetParameter(7,fitFlowSide->GetParameter(0));
    fitFlowTot->SetParameter(8,fitFlowSide->GetParameter(1));
    // fitFlowTot->SetParameter(8,fitFlowSide->GetParameter(0));
    // fitFlowTot->SetParameter(9,fitFlowSide->GetParameter(1));
    // fitFlowTot->SetParameter(12,fitFlowSide->GetParameter(1));
    // fitFlowTot->FixParameter(10,fitFlowSide->GetParameter(2));
  }
  hFlowMass->SetMinimum(-0.1);
  hFlowMass->SetMaximum(0.5);
	hFlowMass->Fit("fitFlowTot","RIB");

	dFlow = fitFlowTot->GetParameter(0);
	dFlowError = fitFlowTot->GetParError(0);

  canFitInvMass->cd(6);

  latex->DrawLatex(0.2,0.83,Form("#color[8]{v_{2} = %.3g#pm%.2g}",dFlow,dFlowError));
  latex->DrawLatex(0.2,0.75,Form("#color[8]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitFlowTot->GetChisquare(), fitFlowTot->GetNDF(),fitFlowTot->GetChisquare()/fitFlowTot->GetNDF()));

	return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::Fatal(TString sMsg, TString sMethod)
{
	printf("Fatal::%s  %s. Terminating!\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Error(TString sMsg, TString sMethod)
{
	printf("Error::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Info(TString sMsg, TString sMethod)
{
	printf("Info::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Warning(TString sMsg, TString sMethod)
{
	printf("Warning::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Debug(TString sMsg, TString sMethod)
{
	if(fbDebug) printf("Debug::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}
//_____________________________________________________________________________
void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}
//_____________________________________________________________________________
void FakeHistosOnlyForExample(TH1* &hstat, TH1* &hsyst, TH1*&hsystCorr) {

  TF1 * fExpo = new TF1 ("fExpo", "expo");
  fExpo->SetParameters(10, -0.3);
  hstat     = new TH1F("hstat", "hstat", 100, 0, 10);
  hsyst     = new TH1F("hsyst", "hsyst", 100, 0, 10);
  hsystCorr = new TH1F("hsystCorr", "hsystCorr", 100, 0, 10);
  hstat     ->Sumw2();
  hsyst     ->Sumw2();
  hsystCorr ->Sumw2();
  hstat->FillRandom("fExpo",20000);
  Int_t nbinx = hstat->GetNbinsX();

  for(Int_t ibinx = 1; ibinx <= nbinx; ibinx++){
    hsyst->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsyst->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.08);
			  hsystCorr->SetBinContent (ibinx, hstat->GetBinContent(ibinx));
			  hsystCorr->SetBinError   (ibinx, hstat->GetBinContent(ibinx)*0.15);
  }

}
//_____________________________________________________________________________
void LoadLibs() {
  gSystem->Load("libCore.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
}
//_____________________________________________________________________________
void DrawLogo (Int_t logo, Double_t xmin, Double_t ymin) {

  // Logo is not needed anymore, now we only write alice preliminary
  // Logo:
  // 0: Justr writes "ALICE" (for final data)
  // Anything eles: writes "ALICE Preliminary"
  TString sText;

  switch (logo)
  {
    case 0:
      sText = TString("ALICE");
      break;

    case 1:
      sText = TString("ALICE Preliminary");
      break;

    case 2:
      sText = TString("ALICE This thesis");
      break;
  }


  TLatex *   tex = new TLatex(xmin,ymin, sText.Data());
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->Draw();

  // OLD logo
  //  TPad * currentPad = gPad;
  // Double_t AliLogo_LowX =xmin;
  // Double_t AliLogo_LowY = ymin;
  // Double_t AliLogo_Height = size;
  // //ALICE logo is a  file that is 821x798 pixels->should be wider than a square
  // Double_t AliLogo_Width  = (821./798.) * AliLogo_Height * gPad->GetWh() / gPad->GetWw();

  // TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",AliLogo_LowX,AliLogo_LowY,AliLogo_LowX+AliLogo_Width,AliLogo_LowY+AliLogo_Height);
  // myPadSetUp(myPadLogo,0,0,0,0);
  // //myPadLogo->SetFixedAspectRatio(1);
  // myPadLogo->Draw();
  // myPadLogo->cd();
  // if (logo == 0) {
  //   myPadLogo->SetFillColor(2); // color to first figure out where is the pad then comment !
  // } else if (logo == 1){
  //   TASImage *myAliceLogo = new TASImage(performanceLogoPath);
  //   myAliceLogo->Draw();
  // } else if (logo == 2) {
  //   TASImage *myAliceLogo = new TASImage(preliminaryLogoPath);
  //   myAliceLogo->Draw();
  // }
  // // go back to the old pad
  // currentPad->cd();

}//_____________________________________________________________________________
