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

    TString     GetSpeciesName(); // system species name (Charged, K0s, ...)
    TString     GetSpeciesLabel(); // readable species name (h^{#pm}, K^{0}_{S}, ...)
    TString     GetEtaGapString() { return TString(Form("%02.2g",10*fEtaGap)); } // used for "character-safe" names

    void        SetHarmonics(Int_t harm) { fHarmonics = harm; }
    void        SetEtaGap(Float_t eta) { fEtaGap = eta; }
    void        SetNumSamples(Short_t num) { fNumSamples = num; }
    void        SetInputTag(const char* name) { fInputTag = name; }
    void        SetPtBins(std::vector<Double_t> array) { fPtBinsEdges = array; fNumPtBins = (Int_t) array.size() - 1; } // setup the pt binning for this task using std::vectors. NB: possible with {}
    void        SetShowMultDist(Bool_t show) { fShowMult = show; }
    void        SetConsiderCorrelations(Bool_t cor = kTRUE) { fConsCorr = cor; }
    void        SetDoFourCorrelations(Bool_t four = kTRUE) { fDoFour = four; }
    void        SetRebinning(Bool_t rebin = kTRUE) { fRebinning = rebin; }
    void        SetMergePosNeg(Bool_t merge = kTRUE) { fMergePosNeg = merge; }
    void        SetDesamplingUseRMS(Bool_t use = kTRUE) { fDesampleUseRMS = use; }
    void        SetProcessMixedHarmonics(TString nameDiff, TString nameRefs) { fProcessMixed = kTRUE; fMixedDiff = nameDiff; fMixedRefs = nameRefs; }

    // fitting
    void        SetInvMassRebin(Short_t rebin = 2) { fRebinInvMass = rebin; }
    void        SetFlowMassRebin(Short_t rebin = 2) { fRebinFlowMass = rebin; }
    // fitting parameters
    void        SetFitPhiSubtLS(Bool_t sub = kTRUE) { fFlowFitPhiSubtLS = sub; }
    void        SetFitMassRange(Double_t dLow, Double_t dHigh) { fFlowFitRangeLow = dLow; fFlowFitRangeHigh = dHigh; }
    void        SetFitMassSig(TString func, Int_t pars) { fFlowFitMassSig = func; fNumParMassSig = pars; }
    void        SetFitMassBG(TString func, Int_t pars) { fFlowFitMassBG = func; fNumParMassBG = pars; }
    void        SetFitFlowBG(TString func, Int_t pars) { fFlowFitFlowBG = func; fNumParFlowBG = pars; }
    void        SetFitParDefaults(Double_t* array, Int_t size);
    void        SetFitParLimitsLow(Double_t* array, Int_t size);
    void        SetFitParLimitsHigh(Double_t* array, Int_t size);

    TString     fTaskTag; // "unique" tag used primarily for storing output
    TString     fName; // task name
    PartSpecies fSpecies; // species involved
    TString     fInputTag; // alterinative tag appended to name of input histos & profiles
    Int_t       fHarmonics; // harmonics
    Double_t    fEtaGap; // eta gap
    Bool_t      fDoFour; // process 4-particle correlations
    Bool_t      fConsCorr; // consider correlations in cumulant / flow calculations
    Int_t       fNumSamples; // [10] number of samples
    Int_t       fNumPtBins; // actual number of pT bins (not size of array) for rebinning
    std::vector<Double_t>   fPtBinsEdges; // pt binning
    Bool_t      fShowMult; // show multiplicity distribution
    Bool_t      fSampleMerging; // [kFALSE] flag for merging TProfiles (good for refs)
    Bool_t      fRebinning; // [kTRUE] flag for rebinning prior to desampling
    Bool_t      fDesampleUseRMS; // [kFALSE] flag for using RMS as uncertainty during desampling
    Bool_t      fMergePosNeg; // [kFALSE] flag for merging results corresponding to positive and negative POIs
    Bool_t      fProcessMixed; // [kFALSE] flag for processing mixed harmonics (non-linear flow modes)
    TString     fMixedDiff; // name (tag) for diff. profile for mixed harmonics
    TString     fMixedRefs; // name (tag) for reference profile for mixed harmonics
    // Reconstructed fitting
    Bool_t      fFlowFitPhiSubtLS; // [kFALSE] flag for subtraction of like-sign background from the unlike-sign one
    Short_t     fRebinInvMass; // flag for rebinning inv-mass (and BG) histo
    Short_t     fRebinFlowMass; // flag for rebinning flow-mass profile

    Double_t    fFlowFitRangeLow; // lower edge for fitting during flow extraction
    Double_t    fFlowFitRangeHigh; // high edge for fitting during flow extraction
    TString     fFlowFitMassSig; // formula for signal contribution of inv.mass fit
    TString     fFlowFitMassBG; // formula for bg contribution of inv.mass fit
    TString     fFlowFitFlowBG; // formula for bg contribution of flow-mass fit
    Int_t       fNumParMassSig; // number of parameters in 'fFlowFitMassSig'
    Int_t       fNumParMassBG; // number of parameters in 'fFlowFitMassBG'
    Int_t       fNumParFlowBG; // number of parameters in 'fFlowFitFlowBG'
    static const Int_t fNumParsMax = 20;
    Double_t    fFitParDefaults[fNumParsMax]; // default values for all parameters in all formulas (i.e. master flow-mass formula)
    Double_t    fFitParLimLow[fNumParsMax]; // low limit values for all parameters in all formulas (i.e. master flow-mass formula)
    Double_t    fFitParLimHigh[fNumParsMax]; // high limit values for all parameters in all formulas (i.e. master flow-mass formula)
    TList*      fListProfiles; // TList with profiles slices
    TList*      fListHistos; // TList with histos slices


    std::vector<TH1D*>* fVecHistInvMass; // container for sliced inv. mass projections
    std::vector<TH1D*>* fVecHistInvMassBG; // container for sliced inv. mass projections for BG candidates (phi)
    std::vector<TH1D*>* fVecHistFlowMass; // container for sliced flow-mass projections
    std::vector<TH1D*>* fVecHistFlowMassFour; // container for sliced flow-mass projections
    TCanvas*     fCanvas; // temporary canvas for mass plotting

  protected:
  private:

};

//_____________________________________________________________________________
FlowTask::FlowTask(PartSpecies species, const char* name) :
  fHarmonics(0),
  fEtaGap(0),
  fNumPtBins(-1),
  fPtBinsEdges()
{
  fName = name;
  fSpecies = species;
  fInputTag = "";
  fNumSamples = 10;
  fDoFour = kFALSE;
  fConsCorr = kFALSE;
  fShowMult = kFALSE;
  fRebinning = kTRUE;
  fSampleMerging = kFALSE;
  fDesampleUseRMS = kFALSE;
  fProcessMixed = kFALSE;
  fMixedDiff = "";
  fMixedRefs = "";
  fMergePosNeg = kFALSE;
  fFlowFitPhiSubtLS = kFALSE;
  fRebinFlowMass = 0;
  fRebinInvMass = 0;
  fFlowFitRangeLow = -1.0;
  fFlowFitRangeHigh = -1.0;
  fFlowFitMassSig = TString();
  fFlowFitMassBG = TString();
  fFlowFitFlowBG = TString();
  fNumParMassSig = 0;
  fNumParMassBG = 0;
  fNumParFlowBG = 0;
  fListProfiles = new TList();
  fListProfiles->SetOwner(kTRUE);
  fListHistos = new TList();
  fListHistos->SetOwner(kTRUE);
  fVecHistInvMass = new std::vector<TH1D*>;
  fVecHistInvMassBG = new std::vector<TH1D*>;
  fVecHistFlowMass = new std::vector<TH1D*>;
  fVecHistFlowMassFour = new std::vector<TH1D*>;

  fTaskTag = this->GetSpeciesName();
  if(!fName.EqualTo("")) fTaskTag.Append(Form("_%s",fName.Data()));
}
//_____________________________________________________________________________
FlowTask::~FlowTask()
{
  if(fListProfiles) { fListProfiles->Clear(); delete fListProfiles; }
  if(fListHistos) { fListHistos->Clear(); delete fListHistos; }
  if(fVecHistFlowMass) delete fVecHistFlowMass;
  if(fVecHistFlowMassFour) delete fVecHistFlowMassFour;
  if(fVecHistInvMass) delete fVecHistInvMass;
  if(fVecHistInvMassBG) delete fVecHistInvMassBG;
  if(fCanvas) delete fCanvas;
}
//_____________________________________________________________________________
void FlowTask::SetFitParDefaults(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { Error("Wrong size of parameters array.","SetFitParDefaults"); return; }
  if(!array) { Error("Wrong array.","SetFitParDefaults"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParDefaults[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
void FlowTask::SetFitParLimitsLow(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { Error("Wrong size of parameters array.","SetFitParLimitsLow"); return; }
  if(!array) { Error("Wrong array.","SetFitParLimitsLow"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParLimLow[i] = array[i]; }

  return;
}
//_____________________________________________________________________________
void FlowTask::SetFitParLimitsHigh(Double_t* array, Int_t size)
{
  if(size < 0 || size > fNumParsMax) { Error("Wrong size of parameters array.","SetFitParLimitsHigh"); return; }
  if(!array) { Error("Wrong array.","SetFitParLimitsHigh"); return; }

  for(Int_t i(0); i < size; ++i) { fFitParLimHigh[i] = array[i]; }

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
TString FlowTask::GetSpeciesLabel()
{
  TString label = TString();
  switch (fSpecies)
  {
    case kRefs : label.Append("RFP"); break;
    case kCharged : label.Append("h^{#pm}"); break;
    case kPion : label.Append("#pi^{#pm}"); break;
    case kKaon : label.Append("K^{#pm}"); break;
    case kProton : label.Append("p/{#bar{p}}"); break;
    case kPhi : label.Append("#phi"); break;
    case kK0s : label.Append("K_{S}^{0}"); break;
    case kLambda : label.Append("#Lambda/#bar{#Lambda}"); break;
    default: label.Append("Non");
  }
  return label;
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
  printf("   fProcessMixed: %s\n", fProcessMixed ? "true" : "false");
  printf("   fDoFour: %s\n", fDoFour ? "true" : "false");
  printf("   fShowMult: %s\n", fShowMult ? "true" : "false");
  printf("   fMergePosNeg: %s\n", fMergePosNeg ? "true" : "false");
  printf("   fFlowFitRangeLow: %g\n",fFlowFitRangeLow);
  printf("   fFlowFitRangeHigh: %g\n",fFlowFitRangeHigh);
  printf("   fNumPtBins: %d\n",fNumPtBins);
  if(fNumPtBins > 1) { printf("   fPtBinsEdges: "); for(Int_t i(0); i < (Int_t) fPtBinsEdges.size() ; ++i) printf("%g ",fPtBinsEdges[i]); printf("\n"); }
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
    void        SetSaveSlices(Bool_t save = kTRUE) { fbSaveSlices = save; } // save prepared slices to output file
    void        SetNoFitting(Bool_t save = kTRUE) { fbNoFitting = save; } //
    void        SetGlobalProfNameLabel(const char* label = "") { fsGlobalProfNameLabel = label; } // add global profile label for all tasks NOTE: for the purpose of Flow sub
    void        SetSaveMult(Bool_t bSave = kTRUE) { fbSaveMult = bSave; } // save reference multiplicity
    void        SetMultiplicityBins(std::vector<Double_t> array) { fdMultBins = array; fiNumMultBins = (Int_t) array.size() - 1; } // setup the global multiplicity binning, where size is number of elements in array
    void        SetFitCumulants(Bool_t cum = kTRUE) { fFlowFitCumulants = cum; } // use cn{2} vs m_inv instead of vn{2} vs. m_inv
    void        SetSaveInterSteps(Bool_t save = kTRUE) { fSaveInterSteps = save; }
    void        SetDebug(Bool_t debug = kTRUE) { fbDebug = debug; }
    void        AddTask(FlowTask* task = 0x0); // add task to internal lists of all tasks
    void        Run(); // running the task (main body of the class)
    void        Clear(); // clearing (removing tasks, etc.) after running
  protected:

  private:
    Bool_t      Initialize(); // initialization task
    Bool_t      LoadLists(); // loading flow lists from input file

    Bool_t      ProcessMixed(FlowTask* task); // prepare FlowTask input for mixed harmonics
    Bool_t      ProcessTask(FlowTask* task); // process FlowTask according to it setting
    Bool_t      ProcessRefs(FlowTask* task); // process reference flow task
    Bool_t      ProcessDirect(FlowTask* task, Short_t iMultBin = 0); // process PID (pion,kaon,proton) flow task
    Bool_t      ProcessReconstructed(FlowTask* task, Short_t iMultBin = 0); // process  V0s flow
    Bool_t      PrepareSlices(const Short_t multBin, FlowTask* task, TProfile3D* p3Cor = 0x0, TH3D* h3Entries = 0x0, TH3D* h3EntriesBG = 0x0, TProfile3D* p3CorFour = 0x0); // prepare
    Bool_t      PrepareSlicesNew(FlowTask* task, TString histName, Bool_t bDoCand = kTRUE); // wrapper for making/preparing per-task slices
    Bool_t      MakeProfileSlices(FlowTask* task, TH1* inputProf, TList* outList); // prepare slices out of inputHist
    Bool_t      MakeSparseSlices(FlowTask* task, THnSparse* inputSparse, TList* outList, const char* outName = "hInvMass"); // prepare slices out of 'inputSparse'
    Bool_t      SaveSlices(FlowTask* task); // save prepared slices into a (sepated) output file

    TH1D*       CalcRefCumTwo(TProfile* hTwoRef, FlowTask* task); // calculate cn{2} out of correlation
    TH1D*       CalcRefCumFour(TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate cn{4} out of correlation
    TH1D*       CalcDifCumTwo(TH1D* hTwoDif, FlowTask* task); // calculate dn{2} out of correlation
    TH1D*       CalcDifCumTwo(TProfile* hTwoDif, FlowTask* task); // calculate dn{2} out of correlation
    TH1D*       CalcDifCumFour(TH1D* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate dn{4} out of correlation
    TH1D*       CalcDifCumFour(TProfile* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate dn{4} out of correlation

    TH1D*       CalcRefFlowTwo(TH1D* hTwoRef, FlowTask* task); // calculate vn{2} out of cn{2}
    TH1D*       CalcRefFlowFour(TH1D* hFourRef, FlowTask* task); // calculate vn{4} out of cn{4}
    TH1D*       CalcDifFlowTwo(TH1D* hTwoDif, TH1D* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate vn'{2} out of dn{2} & vn{2}
    TH1D*       CalcDifFlowFour(TH1D* hFourDif, TH1D* hFourRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel = kFALSE); // calculate vn'{4} out of dn{4} and vn{4}


    Bool_t      FitInvMass(TH1* hist, FlowTask* task, TF1& fitOutSig, TF1& fitOutBg);
    Bool_t      FitCorrelations(TH1* hist, FlowTask* task, TF1& fitOutSig, TF1& fitOutBg, TF1& fitInSig, TF1& fitInBg);

    Bool_t 	    ExtractFlowOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowPhiOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowK0sOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for K0s candidates
    Bool_t 	    ExtractFlowLambdaOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits); // extract flow via flow-mass method for Lambda candidates

    TH1*        MergeListProfiles(TList* list); // merge list of TProfiles into single TProfile
    TH1D*       DesampleList(TList* list, TH1D* merged, FlowTask* task, TString name, Bool_t bSkipDesampling = kFALSE); // Desample list of samples for estimating the uncertanity
    Bool_t      PlotDesamplingQA(TList* list, TH1D* hDesampled, FlowTask* task); // produce QA plots for result of desampling procedure
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

    std::vector<Double_t>    fdMultBins; // global multiplicity/centrality binning
    Int_t     fiNumMultBins; // number of multiplicity bins (not size of array)

    TString     fsInputFilePath; // path to the input folder with input file
    TString     fsInputFileName; // name of input file
    TString     fsOutputFilePath; // current (working) path to the output folder (process/task dependent)
    TString     fsOutputFilePathRoot; // root (global/top-level) path to the output folder
    TString     fsOutputFileName; // name of output file
    TString     fsOutputFileMode; // [RECREATE] mode of output file
    TString     fsTaskName; // name of task (inchluded in data structure names)
    TString     fsOutputFileFormat; // [pdf] format of output files (pictures)
    TString     fsGlobalProfNameLabel; // global profile label for all task
    Bool_t      fbSaveSlices; // [kTRUE] save slices of reconstructed particles
    Bool_t      fbNoFitting; // [kFALSE] skip fitting part for reconstructed for debugging processes or for Slicing only
    Bool_t      fbSaveMult; // [kFALSE]
    Bool_t      fFlowFitCumulants; // [kFALSE]

    Bool_t      fbInit; // flag for initialization status
    Bool_t      fbDebug; // flag for debugging : if kTRUE Debug() messages are displayed
    Bool_t      fSaveInterSteps; // flag for saving intermediate steps (correlations, cumulants) into a main output file
    TFile*      ffInputFile; //! input file container
    TFile*      ffOutputFile; //! output file container
    TFile*      ffDesampleFile; //! output file for results of desampling
    TFile*      ffSlicesFile; //! output file for slices
    TFile*      ffFitsFile; //! output file for fitting procedure
    TList*      flFlowRefs; //! TList from input file with RFPs flow profiles
    TList*      flFlowCharged; //! TList from input file with Charged flow profiles
    TList*      flFlowPion; //! TList from input file with pion flow profiles
    TList*      flFlowKaon; //! TList from input file with kaon flow profiles
    TList*      flFlowProton; //! TList from input file with proton flow profiles
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
  fbSaveSlices(kTRUE),
  fbNoFitting(kFALSE),
  fbSaveMult(kFALSE),
  fFlowFitCumulants(kFALSE),
  fSaveInterSteps(kFALSE),
  fdMultBins(),
  fiNumMultBins(0),
  ffInputFile(0x0),
  ffOutputFile(0x0),
  ffSlicesFile(0x0),
  ffFitsFile(0x0),
  ffDesampleFile(0x0),
  flFlowRefs(0x0),
  flFlowCharged(0x0),
  flFlowPion(0x0),
  flFlowKaon(0x0),
  flFlowProton(0x0),
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
}
//_____________________________________________________________________________
ProcessUniFlow::~ProcessUniFlow()
{
  // default destructor
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;
  if(ffFitsFile) delete ffFitsFile;

  if(flFlowRefs) delete flFlowRefs;
  if(flFlowCharged) delete flFlowCharged;
  if(flFlowPion) delete flFlowPion;
  if(flFlowKaon) delete flFlowKaon;
  if(flFlowProton) delete flFlowProton;
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
  flFlowPion = 0x0;
  flFlowKaon = 0x0;
  flFlowProton = 0x0;
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
    if(!ProcessTask(currentTask)) {return;}
  }

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::Initialize()
{
  // initialization of all necessery prerequisits
  if(fiNumMultBins < 1) { Error(Form("Not enough mult. bins: %d (at least 1 needed)!",fiNumMultBins),"Initialize"); return kFALSE; }

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
  ffDesampleFile = TFile::Open(Form("%s/desampling.root",fsOutputFilePath.Data()),fsOutputFileMode.Data());
  if(!ffDesampleFile) { Fatal(Form("Output desampling file '%s/desampling.root' not open!","Initialize")); return kFALSE; }

  if(fbSaveSlices) {
    ffSlicesFile = TFile::Open(Form("%s/slices.root",fsOutputFilePath.Data()),fsOutputFileMode.Data());
    if(!ffSlicesFile)  { Fatal(Form("Output slices file '%s/slices.root' not open!","Initialize")); return kFALSE; }
  }

  // creating output file for fits
  ffFitsFile = TFile::Open(Form("%s/fits.root",fsOutputFilePath.Data()),fsOutputFileMode.Data());
  if(!ffFitsFile) { Fatal(Form("Output desampling file '%s/fits.root' not open!","Initialize")); return kFALSE; }

  Info("Files loaded","Initialize");

  if(!LoadLists()) return kFALSE;
  Info("Flow lists loaded","Initialize");

  // initialization succesfull
  Info("Initialization succesfull","Initialize");
  return kTRUE;
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
  flFlowPion = (TList*) gDirectory->Get(Form("Flow_Pion_%s",fsTaskName.Data()));
  if(!flFlowPion) { Fatal("'flFlowPion' list does not exists!","LoadLists"); return kFALSE; }
  flFlowKaon = (TList*) gDirectory->Get(Form("Flow_Kaon_%s",fsTaskName.Data()));
  if(!flFlowKaon) { Fatal("'flFlowKaon' list does not exists!","LoadLists"); return kFALSE; }
  flFlowProton = (TList*) gDirectory->Get(Form("Flow_Proton_%s",fsTaskName.Data()));
  if(!flFlowProton) { Fatal("'flFlowProton' list does not exists!","LoadLists"); return kFALSE; }
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
Bool_t ProcessUniFlow::ProcessTask(FlowTask* task)
{
  fsOutputFilePath = Form("%s/%s",fsOutputFilePathRoot.Data(),task->fTaskTag.Data());
  gSystem->mkdir(fsOutputFilePath.Data(),kTRUE);

  Info(Form("Processing task: %s",task->fName.Data()),"ProcessTask");
  if(!task) { Error("Task not valid!","ProcessTask"); return kFALSE; }

  task->PrintTask();

  // task checks & initialization
  if(task->fNumPtBins < 1) { Error(Form("Too small number of bins: %d (at least 1 needed)!",task->fNumPtBins),"ProcesTask"); return kFALSE; }
  if(task->fEtaGap < 0.0 && task->fMergePosNeg) { task->fMergePosNeg = kFALSE; Warning("Merging Pos&Neg 'fMergePosNeg' switch off (no gap)","ProcessTask"); }

  // processing mixed
  if(task->fProcessMixed)
  {
    TList* listSlicesProfiles = task->fListProfiles;
    TList* listSlicesHistos = task->fListHistos;

    if(!PrepareSlicesNew(task,task->fMixedDiff)) { Error("Preparing slices failed!","ProcessTask"); return kFALSE; }

    ffOutputFile->cd();
    if(fSaveInterSteps) {
      listSlicesProfiles->Write("MakeProfileSlices",TObject::kSingleKey);
      listSlicesHistos->Write("MakeHistosSlices",TObject::kSingleKey);
    }
    if(fbNoFitting) { Warning("No fitting flag is ON, skipping Processing!","ProcessTask"); return kTRUE; }
    listSlicesProfiles->Write("MakeProfileSlices",TObject::kSingleKey);
    listSlicesHistos->Write("MakeHistosSlices",TObject::kSingleKey);

    if(!ProcessMixed(task)) { Error("ProcessMixed failed!","ProcessTask"); return kFALSE; }

    return kTRUE;
  }

  // processing standard cumulants
  switch (task->fSpecies)
  {
    case FlowTask::kRefs:
      if(!ProcessRefs(task)) { Error(Form("Task '%s' (%s) not processed correctly!",task->fName.Data(), task->GetSpeciesName().Data()),"ProcessTask"); return kFALSE; }
    break;

    case FlowTask::kCharged:
    case FlowTask::kPion:
    case FlowTask::kKaon:
    case FlowTask::kProton:
      for(Short_t binMult(0); binMult < fiNumMultBins; binMult++) { if(!ProcessDirect(task,binMult)) { Error(Form("Task '%s' (%s; mult. bin %d) not processed correctly!",task->fName.Data(),task->GetSpeciesName().Data(),binMult),"ProcessTask"); return kFALSE; }}
    break;

    case FlowTask::kPhi:
    case FlowTask::kK0s:
    case FlowTask::kLambda:
      if(!ProcessReconstructed(task,0)) { Error(Form("Task '%s' (%s) not processed correctly!",task->fName.Data(),task->GetSpeciesName().Data()),"ProcessTask"); return kFALSE; }
    break;

    default:
    break;
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessMixed(FlowTask* task)
{
  if(!task) { Error("Task not valid!","ProcessMixed"); return kFALSE; }
  Debug("Processing mixed","ProcessMixed");

  FlowTask::PartSpecies species = task->fSpecies;
  Bool_t bReco = kFALSE;
  if(species == FlowTask::kK0s || species == FlowTask::kLambda || species == FlowTask::kPhi) { bReco = kTRUE; }
  // if(bReco) { Error("Currently fully hardcoded for validation purposes for Direct only","ProcessMixed"); return kFALSE; }

  TList* listPOIs = 0x0;
  switch (species)
  {
    case FlowTask::kCharged: listPOIs = flFlowCharged; break;
    case FlowTask::kPion: listPOIs = flFlowPion; break;
    case FlowTask::kKaon: listPOIs = flFlowKaon; break;
    case FlowTask::kProton: listPOIs = flFlowProton; break;
    case FlowTask::kK0s: listPOIs = flFlowK0s; break;
    case FlowTask::kLambda: listPOIs = flFlowLambda; break;
    case FlowTask::kPhi: listPOIs = flFlowPhi; break;
    default: Error("Invalid species!","ProcessMixed"); return kFALSE;
  }

  Int_t iSample = 0;
  TString sNameRefs = Form("%s_Pos_sample%d", task->fMixedRefs.Data(),iSample);
  TString sNamePOIs = Form("%s_Pos_sample%d", task->fMixedDiff.Data(), iSample);
  TString sNamePOIsNeg = Form("%s_Neg_sample%d", task->fMixedDiff.Data(), iSample);

  TList trashCol;
  trashCol.SetOwner(kTRUE);

  // ### Preparing Refs ###
  TProfile* profRef_preRebin = (TProfile*) flFlowRefs->FindObject(sNameRefs.Data());
  if(!profRef_preRebin) { Error(Form("Refs profile '%s' pre-rebin not found!",sNameRefs.Data()),"ProcessMixed"); flFlowRefs->ls(); return kFALSE; }
  TProfile* profRef = (TProfile*) profRef_preRebin->Rebin(fiNumMultBins,Form("%s_rebin",sNameRefs.Data()),fdMultBins.data());
  if(!profRef) { Error("Refs profile rebinning failed!","ProcessMixed"); return kFALSE; }
  trashCol.Add(profRef);

  // ### Preparing POIs ###
  for(Int_t iMultBin(0); iMultBin < fiNumMultBins; ++iMultBin) {

    TH1D* histFlow = 0x0; // histo with final results
    // TString sHistFlowName = Form("%s_%s_mult%d",task->GetSpeciesName().Data(),sNamePOIs.Data(),iMultBin);
    TString sHistFlowName = Form("%s_%s_mult%d",task->GetSpeciesName().Data(),task->fMixedDiff.Data(),iMultBin);

    // direct species
    if(!bReco) {
      TString sName = Form("%s_Pos_sample0_mult%d",task->fMixedDiff.Data(),iMultBin);
      TProfile* profVn = (TProfile*) task->fListProfiles->FindObject(sName.Data());
      if(!profVn) { Error("Loading slice failed!","ProcessMixed"); task->fListProfiles->ls(); return kFALSE; }

      // Making vn out of cn,dn
      TH1D* histVn = (TH1D*) profVn->ProjectionX();
      trashCol.Add(histVn);

      // dividing POIS / sqrt(refs)
      Double_t dRefCont = profRef->GetBinContent(iMultBin+1);
      Double_t dRefErr = profRef->GetBinError(iMultBin+1);

      for(Int_t bin(0); bin < profVn->GetNbinsX()+1; ++bin) {
        if(dRefCont < 0.0) {
          histVn->SetBinContent(bin, 9999.9);
          histVn->SetBinError(bin, 9999.9);
          continue;
        }

        Double_t dOldCont = histVn->GetBinContent(bin);
        Double_t dOldErr = histVn->GetBinError(bin);

        Double_t dNewCont = dOldCont / TMath::Sqrt(dRefCont);
        Double_t dNewErrSq = dOldErr*dOldErr/dRefCont + 0.25*TMath::Power(dRefCont,-3.0)*dOldCont*dOldCont*dRefErr*dRefErr;

        histVn->SetBinContent(bin, dNewCont);
        histVn->SetBinError(bin, TMath::Sqrt(dNewErrSq));
      }

      histFlow = histVn;
    } else { // end-if {!bReco}
      // recnstructed
      gSystem->mkdir(Form("%s/fits/",fsOutputFilePath.Data()));

      histFlow = new TH1D(sHistFlowName.Data(),Form("%s: %s; #it{p}_{T} (GeV/#it{c});",task->GetSpeciesLabel().Data(),sHistFlowName.Data()), task->fNumPtBins,task->fPtBinsEdges.data());
      if(!histFlow) { Error("Creation of 'histFlow' failed!","ProcessMixed"); return kFALSE; }
      trashCol.Add(histFlow);

      // dividing POIS / sqrt(refs)
      Double_t dRefCont = profRef->GetBinContent(iMultBin+1);
      Double_t dRefErr = profRef->GetBinError(iMultBin+1);

      for(Int_t iPtBin(0); iPtBin < task->fNumPtBins; ++iPtBin) {

        TH1D* hInvMass = (TH1D*) task->fListHistos->FindObject(Form("hInvMass_mult%d_pt%d",iMultBin,iPtBin));
        if(!hInvMass) { Error("Loading inv. mass slice failed!","ProcessMixed"); task->fListHistos->ls(); return kFALSE; }

        TH1D* hInvMassBg = 0x0;
        if(species == FlowTask::kPhi) {
          hInvMassBg = (TH1D*) task->fListHistos->FindObject(Form("hInvMassBg_mult%d_pt%d",iMultBin,iPtBin));
          if(!hInvMassBg) { Error("Loading inv. mass (Bg) slice failed!","ProcessMixed"); task->fListHistos->ls(); return kFALSE; }
        }

        TString sName = Form("%s_Pos_sample0_mult%d_pt%d",task->fMixedDiff.Data(),iMultBin,iPtBin);
        TProfile* profVn = (TProfile*) task->fListProfiles->FindObject(sName.Data());
        if(!profVn) { Error("Loading correlation slice failed!","ProcessMixed"); task->fListProfiles->ls(); return kFALSE; }

        // Making vn out of cn,dn
        TH1D* histVn = (TH1D*) profVn->ProjectionX();
        trashCol.Add(histVn);

        for(Int_t bin(0); bin < profVn->GetNbinsX()+1; ++bin) {
          if(dRefCont < 0.0) {
            histVn->SetBinContent(bin, 9999.9);
            histVn->SetBinError(bin, 9999.9);
            continue;
          }

          Double_t dOldCont = histVn->GetBinContent(bin);
          Double_t dOldErr = histVn->GetBinError(bin);

          Double_t dNewCont = dOldCont / TMath::Sqrt(dRefCont);
          Double_t dNewErrSq = dOldErr*dOldErr/dRefCont + 0.25*TMath::Power(dRefCont,-3.0)*dOldCont*dOldCont*dRefErr*dRefErr;

          histVn->SetBinContent(bin, dNewCont);
          histVn->SetBinError(bin, TMath::Sqrt(dNewErrSq));
        }

        // Here ready for fitting

        TCanvas* canFitInvMass = new TCanvas("canFitInvMass","canFitInvMass",1600,1200); // canvas for fitting results

        TList* listFits = new TList();
        // listFits->SetOwner(kTRUE); // NB: when on, seg fault happen

        Double_t dFlow = 0.0;
        Double_t dFlowError = 0.0;

        // Bool_t bExtracted = ExtractFlowOneGo(task,hInvMass,hInvMassBg,histVn,dFlow,dFlowError,canFitInvMass,listFits);

        TF1 fitMassSig, fitMassBg, fitFlowSig, fitFlowBg;

        Bool_t bFitMass = FitInvMass(hInvMass, task, fitMassSig, fitMassBg);

        if(!bFitMass) {
          Warning("Fitting inv.mass unsuccesfull","ProcessMixed");
          delete canFitInvMass;
          delete listFits;
          return kFALSE;
        }

        listFits->Add(hInvMass);
        listFits->Add(&fitMassSig);
        listFits->Add(&fitMassBg);

        Bool_t bFitFlow = FitCorrelations(histVn, task, fitFlowSig, fitFlowBg, fitMassSig, fitMassBg);


        if(!bFitFlow) {
          Warning("Fitting flow unsuccesfull","ProcessMixed");
          delete canFitInvMass;
          delete listFits;
          return kFALSE;
        }

        listFits->Add(histVn);
        listFits->Add(&fitFlowSig);
        listFits->Add(&fitFlowBg);

        Int_t iParFlow = fitFlowSig.GetNpar() - 1;
        dFlow = fitFlowSig.GetParameter(iParFlow);
        dFlowError = fitFlowSig.GetParError(iParFlow);

        Double_t dFlowRel = -999.9; if(TMath::Abs(dFlow) > 0.0) { dFlowRel = dFlowError / dFlow; }
        Info(Form("Final v(n,m,k): (mult %d | pt %d) %g +- %g (rel. %.3f)",iMultBin,iPtBin,dFlow,dFlowError,dFlowRel), "ProcessMixed");

        histFlow->SetBinContent(iPtBin+1,dFlow);
        histFlow->SetBinError(iPtBin+1,dFlowError);

        ffFitsFile->cd();
        listFits->Write(Form("fits_%s_cent%d_pt%d",task->GetSpeciesName().Data(),iMultBin,iPtBin),TObject::kSingleKey);

        // === Plotting fits ===
        TLatex latex2;
        // latex2.SetTextFont(43);
        // latex2.SetTextSize(40);
        latex2.SetNDC();

        canFitInvMass->cd(1);
        // if(task->fSpecies == FlowTask::kPhi) canFitInvMass->cd(2);
        latex2.DrawLatex(0.17,0.85,Form("#color[9]{pt %g-%g GeV/c (%g-%g%%)}",task->fPtBinsEdges[iPtBin],task->fPtBinsEdges[iPtBin+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));
        canFitInvMass->cd(2);
        latex2.DrawLatex(0.17,0.85,Form("#color[9]{pt %g-%g GeV/c (%g-%g%%)}",task->fPtBinsEdges[iPtBin],task->fPtBinsEdges[iPtBin+1],fdMultBins[iMultBin],fdMultBins[iMultBin+1]));
        canFitInvMass->SaveAs(Form("%s/fits/%s_%s_mult%d_pt%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),sNamePOIs.Data(),iMultBin,iPtBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

        delete canFitInvMass;
        delete listFits;
      } // end-for {iPtBin}
    } // end-else {!bReco}

    histFlow->SetName(sHistFlowName.Data());
    histFlow->SetTitle(sHistFlowName.Data());

    ffOutputFile->cd();
    histFlow->Write();

    TCanvas* cFlow = new TCanvas("cFlow","cFlow");
    cFlow->cd();
    histFlow->SetStats(0);
    histFlow->DrawCopy();
    cFlow->SaveAs(Form("%s/Flow_%s_%s_mult%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),sNamePOIs.Data(),iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
    delete cFlow;
  } // end-for {iMultBin}

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessRefs(FlowTask* task)
{
  Info("Processing Refs task","ProcesRefs");
  if(!task) { Error("Task not valid!","ProcessRefs"); return kFALSE; }
  if(task->fSpecies != FlowTask::kRefs) { Error("Task species not kRefs!","ProcessRefs"); return kFALSE; }

  Bool_t bDoFour = task->fDoFour; // check if cn{4} should be processed
  Bool_t bCorrelated = task->fConsCorr; // check if correlated uncrt. are considered

  // for saving profiles into TList for merging (into single one) -> estimation for central values
  TList* listCorTwo = new TList(); TString nameCorTwo = Form("pCor2_Refs_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());
  TList* listCorFour = new TList(); TString nameCorFour = Form("pCor4_Refs_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());

  // for cumulants (applicable for diff. flow)
  TList* listCumTwo = new TList(); TString nameCumTwo = Form("hCum2_Refs_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());
  TList* listCumFour = new TList(); TString nameCumFour = Form("hCum4_Refs_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());

  // for vns desampling
  TList* listFlowTwo = new TList(); TString nameFlowTwo = Form("hFlow2_Refs_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());
  TList* listFlowFour = new TList(); TString nameFlowFour = Form("hFlow4_Refs_harm%d_gap%s",task->fHarmonics, task->GetEtaGapString().Data());

  // estimating <multiplicity>
  if(fbSaveMult)
  {
    TProfile* profMult = (TProfile*) flQACharged->FindObject(Form("fpRefsMult"));
    if(!profMult) { Error("MeanMult profile not found!"); flFlowRefs->ls(); return kFALSE; }
    TProfile* profMult_rebin = (TProfile*) profMult->Rebin(fiNumMultBins,Form("%s_rebin",profMult->GetName()),fdMultBins.data());

    ffOutputFile->cd();
    profMult_rebin->Write(profMult_rebin->GetName());
  }

  // new naming convention for input histos (from FlowTask)
  TString sProfTwoName = Form("<<2>>(%d,-%d)",task->fHarmonics, task->fHarmonics);
  TString sProfFourName = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarmonics, task->fHarmonics, task->fHarmonics, task->fHarmonics);
  if(task->fEtaGap > -1.0) {
    sProfTwoName += Form("_2sub(%.2g)",task->fEtaGap);
    sProfFourName += Form("_2sub(%.2g)",task->fEtaGap);
  }

  Debug("Processing samples","ProcessRefs");
  for(Short_t iSample(0); iSample < task->fNumSamples; ++iSample)
  {
    TProfile* pCorTwo = (TProfile*) flFlowRefs->FindObject(Form("%s_Pos_sample%d",sProfTwoName.Data(), iSample));
    if(!pCorTwo) { Warning(Form("Profile '%s' not valid",Form("%s_Pos_sample%d",sProfTwoName.Data(), iSample)),"ProcesRefs"); flFlowRefs->ls(); return kFALSE; }
    // TProfile* pCorTwo = (TProfile*) flFlowRefs->FindObject(Form("fpRefs_%s<2>_harm%d_gap%s_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iSample));
    // if(!pCorTwo) { Warning(Form("Profile 'pCorTwo' (sample %d) not valid",iSample),"ProcesRefs"); flFlowRefs->ls(); return kFALSE; }

    // Process 4-particle correlations
    TProfile* pCorFour = 0x0;
    if(bDoFour)
    {
      pCorFour = (TProfile*) flFlowRefs->FindObject(Form("%s_Pos_sample%d",sProfFourName.Data(),iSample));
      if(!pCorFour) { Warning(Form("Profile '%s' not valid!",Form("%s_Pos_sample%d",sProfFourName.Data(),iSample)),"ProcesRefs"); flFlowRefs->ls(); return kFALSE; }
      // pCorFour = (TProfile*) flFlowRefs->FindObject(Form("fpRefs_%s<4>_harm%d_gap%s_sample%d",fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iSample));
      // if(!pCorFour) { Warning(Form("Profile 'pCorFour' (sample %d) not valid!",iSample),"ProcesRefs"); flFlowRefs->ls(); return kFALSE; }
    }

    // rebinning the profiles
    if(task->fRebinning)
    {
      pCorTwo = (TProfile*) pCorTwo->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorTwo.Data(), iSample),fdMultBins.data());
      if(bDoFour) { pCorFour = (TProfile*) pCorFour->Rebin(fiNumMultBins,Form("%s_sample%d_rebin", nameCorFour.Data(), iSample),fdMultBins.data()); }
    }

    // naming <<X>>
    TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form("{|#Delta#eta| > %g}",task->fEtaGap)); }
    pCorTwo->SetTitle(Form("%s: <<2>>_{%d} %s",task->GetSpeciesName().Data(), task->fHarmonics, sGap.Data()));
    listCorTwo->Add(pCorTwo);
    if(bDoFour)
    {
      pCorFour->SetTitle(Form("%s: <<4>>_{%d} %s",task->GetSpeciesName().Data(), task->fHarmonics, sGap.Data()));
      listCorFour->Add(pCorFour);
    }

    // Making cumulants out of correlations : <<N>>_n -> c_n{N} -> v_n{N}
    // cn{2}
    TH1D* hCumTwo = CalcRefCumTwo(pCorTwo,task);
    if(!hCumTwo) { Error(Form("cn{2} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
    hCumTwo->SetName(Form("%s_sample%d", nameCumTwo.Data(), iSample));
    listCumTwo->Add(hCumTwo);

    // vn{2}
    TH1D* hFlowTwo = CalcRefFlowTwo(hCumTwo,task);
    if(!hFlowTwo) { Error(Form("vn{2} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
    hFlowTwo->SetName(Form("%s_sample%d", nameFlowTwo.Data(), iSample));
    listFlowTwo->Add(hFlowTwo);

    if(task->fDoFour)
    {
      // cn{4}
      TH1D* hCumFour = CalcRefCumFour(pCorFour, pCorTwo, task, bCorrelated);
      if(!hCumFour) { Error(Form("cn{4} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
      hCumFour->SetName(Form("%s_sample%d", nameCumFour.Data(), iSample));
      listCumFour->Add(hCumFour);

      // vn{4}
      TH1D* hFlowFour = CalcRefFlowFour(hCumFour, task);
      if(!hFlowFour) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessRefs"); return kFALSE; }
      hFlowFour->SetName(Form("%s_sample%d", nameFlowFour.Data(), iSample));
      listFlowFour->Add(hFlowFour);
    }
  } // end-for {iSample}: samples
  Debug("Samples processing done!","ProcessRefs");

  // merging correlation profiles to get central values
  Debug("Merging correlations for central values", "ProcessRefs");
  TProfile* pCorTwoMerged = (TProfile*) MergeListProfiles(listCorTwo);
  if(!pCorTwoMerged) { Error("Merging of 'pCorTwoMerged' failed!","ProcessRefs"); return kFALSE; }
  pCorTwoMerged->SetName(Form("%s_merged", nameCorTwo.Data()));

  TH1D* hCumTwoMerged = CalcRefCumTwo(pCorTwoMerged, task);
  if(!hCumTwoMerged) { Error(Form("cn{2} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
  hCumTwoMerged->SetName(Form("%s_merged", nameCumTwo.Data()));

  TH1D* hFlowTwoMerged = CalcRefFlowTwo(hCumTwoMerged, task);
  if(!hFlowTwoMerged) { Error(Form("vn{2} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
  hFlowTwoMerged->SetName(Form("%s_merged", nameFlowTwo.Data()));

  TProfile* pCorFourMerged = 0x0;
  TH1D* hCumFourMerged = 0x0;
  TH1D* hFlowFourMerged = 0x0;
  if(bDoFour)
  {
    pCorFourMerged = (TProfile*) MergeListProfiles(listCorFour);
    if(!pCorFourMerged) { Error("Merging of 'pCorFourMerged' failed!","ProcessRefs"); return kFALSE; }
    pCorFourMerged->SetName(Form("%s_merged", nameCorFour.Data()));

    hCumFourMerged = CalcRefCumFour(pCorFourMerged, pCorTwoMerged, task, bCorrelated);
    if(!hCumFourMerged) { Error(Form("cn{4} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
    hCumFourMerged->SetName(Form("%s_merged", nameCumFour.Data()));

    hFlowFourMerged = CalcRefFlowFour(hCumFourMerged, task);
    if(!hFlowFourMerged) { Error(Form("vn{4} (merged) not processed correctly!"),"ProcessRefs"); return kFALSE; }
    hFlowFourMerged->SetName(Form("%s_merged", nameFlowFour.Data()));
  }

  // desampling
  Debug("Desampling","ProcessRefs");

  TH1D* hCorTwoDesampled = DesampleList(listCorTwo, pCorTwoMerged->ProjectionX(), task, nameCorTwo, kTRUE); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
  if(!hCorTwoDesampled) { Error("Desampling 'hCorTwoDesampled' failed","ProcessRefs"); return kFALSE; }
  hCorTwoDesampled->SetName(nameCorTwo.Data());

  TH1D* hCumTwoDesampled = DesampleList(listCumTwo, hCumTwoMerged, task, nameCumTwo, kTRUE); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
  if(!hCumTwoDesampled) { Error("Desampling 'hCumTwoDesampled' failed","ProcessRefs"); return kFALSE; }
  hCumTwoDesampled->SetName(nameCumTwo.Data());

  TH1D* hFlowTwoDesampled = DesampleList(listFlowTwo, hFlowTwoMerged, task, nameFlowTwo, kTRUE); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
  if(!hFlowTwoDesampled) { Error("Desampling 'hFlowTwoDesampled' failed","ProcessRefs"); return kFALSE; }
  hFlowTwoDesampled->SetName(nameFlowTwo.Data());

  ffOutputFile->cd();
  hCorTwoDesampled->Write();
  hCumTwoDesampled->Write();
  hFlowTwoDesampled->Write();

  if(bDoFour)
  {
    TH1D* hCorFourDesampled = DesampleList(listCorFour, pCorFourMerged->ProjectionX(), task, nameCorFour, kTRUE); // NOTE skipping desampling (last argument kTRUE) for vn{2} -> nothing to de-correlate
    if(!hCorFourDesampled) { Error("Desampling 'hCorFourDesampled' failed","ProcessRefs"); return kFALSE; }
    hCorFourDesampled->SetName(nameCorFour.Data());

    TH1D* hCumFourDesampled = DesampleList(listCumFour, hCumFourMerged, task, nameCumFour);
    if(!hCumFourDesampled) { Error("Desampling 'hCumFourDesampled' failed","ProcessRefs"); return kFALSE; }
    hCumFourDesampled->SetName(nameCumFour.Data());

    TH1D* hFlowFourDesampled = DesampleList(listFlowFour, hFlowFourMerged, task, nameFlowFour);
    if(!hFlowFourDesampled) { Error("Desampling 'hFlowFourDesampled' failed","ProcessRefs"); return kFALSE; }
    hFlowFourDesampled->SetName(nameFlowFour.Data());

    ffOutputFile->cd();
    hCorFourDesampled->Write();
    hCumFourDesampled->Write();
    hFlowFourDesampled->Write();
  }

  // Comment :: Not sure what this is about =====>
  // if(!task->fRebinning)
  // {
  //   // no rebinning
  //   TH1D* hNoRebin_rebinned = TestRebin(hDesampledFlow,task);
  //   hNoRebin_rebinned->Write(Form("%s",hNoRebin_rebinned->GetName()));
  //
  //   if(task->fSampleMerging)
  //   {
  //     TH1D* hMerged_rebinned = TestRebin(hMerged,task);
  //     hMerged_rebinned->Write(Form("%s",hMerged_rebinned->GetName()));
  //   }
  // }
  // <=====

  Debug("Processing done","ProcessRefs");

  if(listCorTwo) delete listCorTwo;
  if(listCumTwo) delete listCumTwo;
  if(listFlowTwo) delete listFlowTwo;
  if(listCorFour) delete listCorFour;
  if(listCumFour) delete listCumFour;
  if(listFlowFour) delete listFlowFour;

  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefCumTwo(TProfile* hTwoRef, FlowTask* task)
{
  // Calculate reference c_n{2} out of correlations
  // NOTE: it is just a fancier Clone(): for consistency
  // cn{2} = <<2>>

  if(!hTwoRef) { Error("Profile 'hTwoRef' not valid!","CalcRefCumTwo"); return 0x0; }
  if(!task) { Error("FlowTask not found!","CalcRefCumTwo"); return 0x0; }

  TH1D* histCum = (TH1D*) hTwoRef->ProjectionX(Form("hCum2_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: c_{%d}{2%s}",task->GetSpeciesName().Data(), task->fHarmonics, sGap.Data()));

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefCumFour(TProfile* hFourRef, TProfile* hTwoRef, FlowTask* task, Bool_t bCorrel)
{
  // Calculate reference c_n{4} out of correlations
  // cn{4} = <<4>> - 2*<<2>>^2

  if(!task) { Error("FlowTask not found!","CalcRefCumFour"); return 0x0; }
  if(!hFourRef) { Error("Profile 'hFourRef' not valid!","CalcRefCumFour"); return 0x0; }
  if(!hTwoRef) { Error("Profile 'hTwoRef' not valid!","CalcRefCumFour"); return 0x0; }
  if(hFourRef->GetNbinsX() != hTwoRef->GetNbinsX()) { Error("Different number of bins!","CalcRefCumFour"); return 0x0; }

  TH1D* histCum = (TH1D*) hFourRef->ProjectionX(Form("hCum4_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: c_{%d}{4%s}",task->GetSpeciesName().Data(), task->fHarmonics, sGap.Data()));
  histCum->Reset();

  for(Int_t iBin(1); iBin < hFourRef->GetNbinsX()+1; ++iBin)
  {
    Double_t dContInFour = hFourRef->GetBinContent(iBin);
    Double_t dErrInFour = hFourRef->GetBinError(iBin);

    Double_t dContInTwo = hTwoRef->GetBinContent(iBin);
    Double_t dErrInTwo = hTwoRef->GetBinError(iBin);

    Double_t dContOut = dContInFour - 2.0 * dContInTwo * dContInTwo;
    histCum->SetBinContent(iBin, dContOut);

    Double_t dErrOutFour = dErrInFour; // wrt. <4>
    Double_t dErrOutTwo = -4.0 * dContInTwo * dErrInTwo; // wrt. <2>

    Double_t dErrOut = TMath::Power(dErrOutFour, 2.0) + TMath::Power(dErrOutTwo, 2.0);
    if(bCorrel) { dErrOut += 2.0 * dErrOutFour * dErrOutTwo; }
    histCum->SetBinError(iBin, TMath::Sqrt(dErrOut));
    Debug(Form("bin %2d : cn{4} = %g - 2*%g^2 = %g", iBin, dContInFour, dContInTwo, dContOut),"CalcRefCumFour");
  }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumTwo(TProfile* hTwoDif, FlowTask* task)
{
  // Same as with TH1D argument.
  // First cast as TH1D and then call CalcDifCumTwo

  TH1D* histCum = CalcDifCumTwo((TH1D*)hTwoDif->ProjectionX("_temp"),task);
  if(!histCum) { Error("Failed!","CalcDifCumTwo"); return 0x0; }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumTwo(TH1D* hTwoDif, FlowTask* task)
{
  // Calculate reference d_n{2} out of correlations
  // NOTE: it is just a fancier Clone(): for consistency
  // dn{2} = <<2'>>

  if(!hTwoDif) { Error("Input 'hTwoDif' not valid!","CalcDifCumTwo"); return 0x0; }
  if(!task) { Error("FlowTask not found!","CalcDifCumTwo"); return 0x0; }

  TH1D* histCum = (TH1D*) hTwoDif->Clone(Form("hCum2_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: d_{%d}{2%s}",task->GetSpeciesLabel().Data(), task->fHarmonics, sGap.Data()));

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumFour(TProfile* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Same as with TH1D argument.
  // First cast as TH1D and then call CalcDifCumFour

  TH1D* histCum = CalcDifCumFour((TH1D*)hFourDif->ProjectionX("_temp"),hTwoDif, hTwoRef, iRefBin, task, bCorrel);
  if(!histCum) { Error("Failed!","CalcDifCumFour"); return 0x0; }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifCumFour(TH1D* hFourDif, TH1* hTwoDif, TH1* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Calculate reference d_n{4} out of correlations
  // dn{4} = <<4'>> - <<2>><<2'>>

  if(!hFourDif) { Error("Input 'hFourDif' not valid!","CalcDifCumFour"); return 0x0; }
  if(!hTwoDif) { Error("Input 'hTwoDif' not valid!","CalcDifCumFour"); return 0x0; }
  if(!hTwoRef) { Error("Input 'hTwoRef' not valid!","CalcDifCumFour"); return 0x0; }
  if(!task) { Error("FlowTask not found!","CalcDifCumFour"); return 0x0; }
  if(iRefBin < 1) { Error("Bin 'iRefBin; < 1!'","CalcDifCumFour"); return 0x0; }
  if(hFourDif->GetNbinsX() != hTwoDif->GetNbinsX()) { Error("Different number of bins!","CalcDifCumFlow"); return 0x0; }

  TH1D* histCum = (TH1D*) hFourDif->Clone(Form("hCum4_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histCum->SetTitle(Form("%s: d_{%d}{4%s}",task->GetSpeciesLabel().Data(), task->fHarmonics, sGap.Data()));
  histCum->Reset();

  Double_t dContInTwoRef = hTwoRef->GetBinContent(iRefBin);
  Double_t dErrInTwoRef = hTwoRef->GetBinError(iRefBin);

  for(Int_t iBin(1); iBin < histCum->GetNbinsX()+1; ++iBin)
  {
    Double_t dContInFourDif = hFourDif->GetBinContent(iBin);
    Double_t dErrInFourDif = hFourDif->GetBinError(iBin);

    Double_t dContInTwoDif = hTwoDif->GetBinContent(iBin);
    Double_t dErrInTwoDif = hTwoDif->GetBinError(iBin);

    Double_t dContOut = dContInFourDif - 2.0 * dContInTwoDif * dContInTwoRef;
    histCum->SetBinContent(iBin, dContOut);

    Double_t dErrOutFour = dErrInFourDif; // wrt. <4>
    Double_t dErrOutTwoDif =  -2.0 * dContInTwoRef * dErrInTwoDif; // wrt. <2'>
    Double_t dErrOutTwoRef =  -2.0 * dContInTwoDif * dErrInTwoRef; // wrt. <2>

    Double_t dErrOutSq = TMath::Power(dErrOutFour, 2.0) + TMath::Power(dErrOutTwoDif, 2.0) + TMath::Power(dErrOutTwoRef, 2.0);
    if(bCorrel) { dErrOutSq += 2.0 * dErrOutFour * dErrOutTwoDif + 2.0 * dErrOutFour * dErrOutTwoRef + 2.0 * dErrOutTwoDif * dErrOutTwoRef; }
    histCum->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
    Debug(Form("bin %2d : dn{4} = %g - 2*%g*%g = %g", iBin, dContInFourDif, dContInTwoDif, dContInTwoRef, dContOut),"CalcDifCumFour");
  }

  return histCum;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefFlowTwo(TH1D* hTwoRef, FlowTask* task)
{
  // Calculate reference v_n{2} out of c_n{2}
  // vn{2} = cn{2}^(1/2)

  if(!hTwoRef) { Error("Histo 'hTwoRef' not valid!","CalcRefFlowTwo"); return 0x0; }
  if(!task) { Error("FlowTask not found!","CalcRefFlowTwo"); return 0x0; }

  TH1D* histFlow = (TH1D*) hTwoRef->Clone(Form("hFlow2_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histFlow->SetTitle(Form("%s: v_{%d}{2%s}",task->GetSpeciesLabel().Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  for(Int_t iBin(1); iBin < hTwoRef->GetNbinsX()+1; ++iBin)
  {
    Double_t dContIn = hTwoRef->GetBinContent(iBin);
    Double_t dErrIn = hTwoRef->GetBinError(iBin);

    if(dContIn > 0.0 && dErrIn >= 0.0)
    {
      Double_t dContOut = TMath::Sqrt(dContIn);
      histFlow->SetBinContent(iBin, dContOut);
      Double_t dErrOutSq = 0.25 * dErrIn * dErrIn / dContIn;
      histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
    }
    else
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }
  }

  return histFlow;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcRefFlowFour(TH1D* hFourRef, FlowTask* task)
{
  // Calculate reference v_n{4} out of c_n{4}
  // vn{4} = (-cn{4})^(1/4)

  if(!hFourRef) { Error("Histo 'hFourRef' not valid!","CalcRefFlowFour"); return 0x0; }
  if(!task) { Error("FlowTask not found!","CalcRefFlowFour"); return 0x0; }

  TH1D* histFlow = (TH1D*) hFourRef->Clone(Form("hFlow4_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histFlow->SetTitle(Form("%s: v_{%d}{4%s}",task->GetSpeciesLabel().Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  for(Int_t iBin(1); iBin < hFourRef->GetNbinsX()+1; ++iBin)
  {
    Double_t dContIn = hFourRef->GetBinContent(iBin);
    Double_t dErrIn = hFourRef->GetBinError(iBin);

    if(dContIn < 0.0 && dErrIn >= 0.0)
    {
      Double_t dContOut = TMath::Power(-1.0 * dContIn, 0.25);
      histFlow->SetBinContent(iBin, dContOut);
      Double_t dErrOutSq = TMath::Power(0.25 * dErrIn * TMath::Power(-dContIn, -0.75), 2.0);
      histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
      Debug(Form("bin %2d : vn{4} = (-1*%g)^(0,25) = %g", iBin, dContIn, dContOut),"CalcRefFlowFour");
    }
    else
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }
  }

  return histFlow;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifFlowTwo(TH1D* hTwoDif, TH1D* hTwoRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Calculate differential v_n{2} out of d_n{2} & v_n{2} (!)
  // vn'{2} = dn{2} / vn{2}

  if(!hTwoDif) { Error("Histo 'hTwoDif' not valid!","CalcDifFlowTwo"); return 0x0; }
  if(!hTwoRef) { Error("Histo 'hTwoRef' not valid!","CalcDifFlowTwo"); return 0x0; }
  if(!task) { Error("FlowTask not found!","CalcDifFlowTwo"); return 0x0; }
  if(iRefBin < 1) { Error("Bin 'iRefBin; < 1!'","CalcDifFlowTwo"); return 0x0; }

  TH1D* histFlow = (TH1D*) hTwoDif->Clone(Form("hFlow2_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histFlow->SetTitle(Form("%s: v_{%d}{2%s}",task->GetSpeciesLabel().Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  Double_t dContInRef = hTwoRef->GetBinContent(iRefBin);
  Double_t dErrInRef = hTwoRef->GetBinError(iRefBin);

  // flow not real -> putting 'wrong' numbers in
  if(dContInRef < -9.0 || (dContInRef <= 0.0 && dErrInRef > 1000))
  {
    for(Int_t iBin(1); iBin < histFlow->GetNbinsX()+1; ++iBin)
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }

    return histFlow;
  }

  // flow real -> correct analytical calculation
  for(Int_t iBin(1); iBin < histFlow->GetNbinsX()+1; ++iBin)
  {
    Double_t dContInDif = hTwoDif->GetBinContent(iBin);
    Double_t dErrInDif = hTwoDif->GetBinError(iBin);

    Double_t dContOut = dContInDif / dContInRef;
    histFlow->SetBinContent(iBin, dContOut);

    Double_t dErrOutDif = dErrInDif / dContInRef;
    Double_t dErrOutRef = -1.0 * dContInDif * TMath::Power(dContInRef, -2.0) * dErrInRef;

    Double_t dErrOutSq = TMath::Power(dErrOutDif, 2.0) + TMath::Power(dErrOutRef, 2.0);
    if(bCorrel) { dErrOutSq += 2.0 * dErrOutDif * dErrOutRef; }
    histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
  }

  return histFlow;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::CalcDifFlowFour(TH1D* hFourDif, TH1D* hFourRef, Int_t iRefBin, FlowTask* task, Bool_t bCorrel)
{
  // Calculate differential v_n{4} out of d_n{4} & v_n{4} (!)
  // vn'{2} = - dn{4} / vn{4}^3

  if(!hFourDif) { Error("Histo 'hFourDif' not valid!","CalcDifFlowFour"); return 0x0; }
  if(!hFourRef) { Error("Histo 'hFourRef' not valid!","CalcDifFlowFour"); return 0x0; }
  if(!task) { Error("FlowTask not found!","CalcDifFlowFour"); return 0x0; }
  if(iRefBin < 1) { Error("Bin 'iRefBin; < 1!","CalcDifFlowFour"); return 0x0; }

  TH1D* histFlow = (TH1D*) hFourDif->Clone(Form("hFlow4_%s_harm%d_gap%s",task->GetSpeciesName().Data(),task->fHarmonics,task->GetEtaGapString().Data()));

  TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form(",|#Delta#eta| > %g",task->fEtaGap)); }
  histFlow->SetTitle(Form("%s: v_{%d}{4%s}",task->GetSpeciesLabel().Data(), task->fHarmonics, sGap.Data()));
  histFlow->Reset();

  Double_t dContInRef = hFourRef->GetBinContent(iRefBin);
  Double_t dErrInRef = hFourRef->GetBinError(iRefBin);

  // flow not real -> putting 'wrong' numbers in
  if(dContInRef <= 0.0 || (dContInRef < -9.0 && dErrInRef > 1000))
  {
    for(Int_t iBin(1); iBin < histFlow->GetNbinsX()+1; ++iBin)
    {
      histFlow->SetBinContent(iBin, -9.9);
      histFlow->SetBinError(iBin, 99999.9);
    }

    return histFlow;
  }

  // flow real -> correct analytical calculation
  for(Int_t iBin(1); iBin < histFlow->GetNbinsX()+1; ++iBin)
  {
    Double_t dContInDif = hFourDif->GetBinContent(iBin);
    Double_t dErrInDif = hFourDif->GetBinError(iBin);

    Double_t dContOut = -1.0 * dContInDif * TMath::Power(dContInRef, -3.0);
    histFlow->SetBinContent(iBin, dContOut);

    Double_t dErrOutDif = -1.0 * dErrInDif * TMath::Power(dContInRef, -3.0); // wrt. dn
    Double_t dErrOutRef = 3.0 * dContInDif * dErrInRef * TMath::Power(dContInRef, -4.0); // wrt. vn

    Double_t dErrOutSq = TMath::Power(dErrOutDif, 2.0) + TMath::Power(dErrOutRef, 2.0);
    if(bCorrel) { dErrOutSq += 2.0 * dErrOutDif * dErrOutRef; }
    histFlow->SetBinError(iBin, TMath::Sqrt(dErrOutSq));
    Debug(Form("bin %2d : vn{4} = -1 * %g / (%g)^3 = %g / %g = %g", iBin, dContInDif, dContInRef, -1.0*dContInDif, TMath::Power(dContInRef, +3.0),  dContOut), "CalcDifFlowFour");
  }

  return histFlow;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessDirect(FlowTask* task, Short_t iMultBin)
{
  Info(Form("Processing direct task (mult. bin %d)", iMultBin),"ProcesDirect");
  if(!task) { Error("Task not valid!","ProcessDirect"); return kFALSE; }

  Bool_t bDoFour = task->fDoFour; // check if cn{4} should be processed
  Bool_t bCorrelated = task->fConsCorr; // check if correlated uncrt. are considered

  TList* listInput = 0x0;
  switch (task->fSpecies)
  {
    case FlowTask::kCharged:
      listInput = flFlowCharged;
    break;

    case FlowTask::kPion:
      listInput = flFlowPion;
    break;

    case FlowTask::kKaon:
      listInput = flFlowKaon;
    break;

    case FlowTask::kProton:
      listInput = flFlowProton;
    break;

    default:
      Error("Task species not direct!","ProcessDirect");
      return kFALSE;
  }

  if(!listInput) { Error("Input list not loaded!","ProcessDirect"); return kFALSE; }

  // Loading list where reference flow samples are stored
  TList* listRefCorTwo = (TList*) ffDesampleFile->Get(Form("pCor2_Refs_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefCorTwo) { Error("List 'listRefCorTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  TList* listRefTwo = (TList*) ffDesampleFile->Get(Form("hFlow2_Refs_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefTwo) { Error("List 'listRefTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  TList* listRefFour = 0x0;
  if(bDoFour)
  {
    listRefFour = (TList*) ffDesampleFile->Get(Form("hFlow4_Refs_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
    if(!listRefFour) { Error("List 'listRefFour' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }
  }

  // List for desampling
  TList* listCorTwo = new TList(); TString nameCorTwo = Form("pCor2_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listCumTwo = new TList(); TString nameCumTwo = Form("hCum2_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listFlowTwo = new TList(); TString nameFlowTwo = Form("hFlow2_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);

  TList* listCorFour = new TList(); TString nameCorFour = Form("pCor4_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listCumFour = new TList(); TString nameCumFour = Form("hCum4_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
  TList* listFlowFour = new TList(); TString nameFlowFour = Form("hFlow4_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);

  Debug("Processing samples","ProcessDirect");

  // new naming convention for input histos (from FlowTask)
  TString sProfTwoName = Form("<<2>>(%d,-%d)",task->fHarmonics, task->fHarmonics);
  TString sProfFourName = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarmonics, task->fHarmonics, task->fHarmonics, task->fHarmonics);
  if(task->fEtaGap > -1.0) {
    sProfTwoName += Form("_2sub(%.2g)",task->fEtaGap);
    sProfFourName += Form("_2sub(%.2g)",task->fEtaGap);
  }

  for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
  {
    Debug(Form("Processing sample %d",iSample), "ProcessDirect");
    // <<2'>>
    TProfile2D* p2CorTwoDif = 0x0;
    if(task->fMergePosNeg)
    {
      // loading pos & neg if fMergePosNeg is ON
      TProfile2D* prof2pos = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfTwoName.Data(),iSample));
      TProfile2D* prof2neg = (TProfile2D*) listInput->FindObject(Form("%s_Neg_sample%d",sProfTwoName.Data(),iSample));
      if(!prof2pos || !prof2neg) { Error("<2>: Pos & Neg profile merging: 'prof2pos' OR 'prof2neg' not found!.","ProcessDirect"); return kFALSE; }

      // merging pos & neg
      TList* listMerge = new TList();
      listMerge->Add(prof2pos);
      listMerge->Add(prof2neg);
      p2CorTwoDif = (TProfile2D*) MergeListProfiles(listMerge);
      delete listMerge; // first delete, then check (return)
      if(!p2CorTwoDif) { Error("<2>: Pos & Neg profile merging failed!","ProcessDirect"); return kFALSE; }
    }
    else
    {
      // loading single (Pos) profile
      if(task->fInputTag.EqualTo("")) // loading default-ly named profile
      { p2CorTwoDif = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfTwoName.Data(),iSample)); }
      else // loading "non-standardly" named profile
      { p2CorTwoDif = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<2>_harm%d_gap%s_%s_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),task->fInputTag.Data(),iSample)); }
    }
    if(!p2CorTwoDif) { Error(Form("Profile '%s' (sample %d) does not exists.",sProfTwoName.Data(),iSample),"ProcessDirect"); return kFALSE; }

    // <<4'>>
    TProfile2D* p2CorFourDif = 0x0;
    if(bDoFour)
    {
      if(task->fMergePosNeg)
      {
        // loading pos & neg if fMergePosNeg is ON
        TProfile2D* prof2pos = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfFourName.Data(),iSample));
        TProfile2D* prof2neg = (TProfile2D*) listInput->FindObject(Form("%s_Neg_sample%d",sProfFourName.Data(),iSample));
        if(!prof2pos || !prof2neg) { Error("<4>: Pos & Neg profile merging: 'prof2pos' OR 'prof2neg' not found!.","ProcessDirect"); return kFALSE; }

        // merging pos & neg
        TList* listMerge = new TList();
        listMerge->Add(prof2pos);
        listMerge->Add(prof2neg);
        p2CorFourDif = (TProfile2D*) MergeListProfiles(listMerge);
        delete listMerge; // first delete, then check (return)
        if(!p2CorFourDif) { Error("<4>: Pos & Neg profile merging failed!","ProcessDirect"); return kFALSE; }
      }
      else
      {
        // loading single (Pos) profile
        if(task->fInputTag.EqualTo("")) // loading default-ly named profile
        { p2CorFourDif = (TProfile2D*) listInput->FindObject(Form("%s_Pos_sample%d",sProfFourName.Data(),iSample)); }
        else // loading "non-standardly" named profile
        { p2CorFourDif = (TProfile2D*) listInput->FindObject(Form("fp2%s_%s<4>_harm%d_gap%s_%s_sample%d",task->GetSpeciesName().Data(),fsGlobalProfNameLabel.Data(),task->fHarmonics,task->GetEtaGapString().Data(),task->fInputTag.Data(),iSample)); }
      }
      if(!p2CorFourDif) { Error(Form("Profile '%s' (sample %d) does not exists.",sProfFourName.Data(),iSample),"ProcessDirect"); listInput->ls(); return kFALSE; }
    }

    Debug("Rebinning profiles","ProcessDirect");
    // rebinning according in mult bin
    Short_t binMultLow = p2CorTwoDif->GetXaxis()->FindFixBin(fdMultBins[iMultBin]);
    Short_t binMultHigh = p2CorTwoDif->GetXaxis()->FindFixBin(fdMultBins[iMultBin+1]) - 1;

    TProfile* pCorTwoDif = p2CorTwoDif->ProfileY(nameCorTwo.Data(),binMultLow,binMultHigh);
    TProfile* pCorFourDif = 0x0; if(bDoFour) { pCorFourDif = p2CorFourDif->ProfileY(nameCorFour.Data(),binMultLow,binMultHigh); }

    // rebinning according to pt bins
    if(task->fNumPtBins > 0)
    {
      pCorTwoDif = (TProfile*) pCorTwoDif->Rebin(task->fNumPtBins,Form("%s_sample%d_rebin", nameCorTwo.Data(), iSample), task->fPtBinsEdges.data());
      if(bDoFour) { pCorFourDif = (TProfile*) pCorFourDif->Rebin(task->fNumPtBins,Form("%s_sample%d_rebin", nameCorFour.Data(), iSample), task->fPtBinsEdges.data()); }
    }
    else
    {
      pCorTwoDif = (TProfile*) pCorTwoDif->Clone(Form("%s_sample%d_rebin", nameCorTwo.Data(), iSample));
      if(bDoFour) { pCorFourDif = (TProfile*) pCorFourDif->Clone(Form("%s_sample%d_rebin", nameCorFour.Data(), iSample)); }
    }

    // renaming
    TString sGap = TString(); if(task->fEtaGap > -1.0) { sGap.Append(Form("{|#Delta#eta| > %g}",task->fEtaGap)); }
    pCorTwoDif->SetTitle(Form("%s: <<2'>>_{%d} %s; #it{p}_{T} (GeV/#it{c});",task->GetSpeciesName().Data(), task->fHarmonics, sGap.Data()));
    if(bDoFour) { pCorFourDif->SetTitle(Form("%s: <<4'>>_{%d} %s; #it{p}_{T} (GeV/#it{c});",task->GetSpeciesName().Data(), task->fHarmonics, sGap.Data())); }

    // NOTE: Here the <X'> is ready & rebinned
    listCorTwo->Add(pCorTwoDif);
    if(bDoFour) { listCorFour->Add(pCorFourDif); }

    Debug("Calculating flow","ProcessDirect");
    // loading reference vn{2}
    TH1D* hFlowRefTwo = (TH1D*) listRefTwo->FindObject(Form("hFlow2_Refs_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
    if(!hFlowRefTwo) { Error(Form("Histo 'hFlowRefTwo' (sample %d) does not exists",iSample),"ProcessDirect"); listRefTwo->ls(); ffDesampleFile->ls(); return kFALSE; }

    // dn{2}
    TH1D* hCumTwoDif = CalcDifCumTwo(pCorTwoDif, task);
    if(!hCumTwoDif) { Error(Form("dn{2} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
    hCumTwoDif->SetName(Form("%s_sample%d", nameCumTwo.Data(), iSample));
    listCumTwo->Add(hCumTwoDif);

    // v'n{2}
    TH1D* hFlowTwoDif = CalcDifFlowTwo(hCumTwoDif, hFlowRefTwo, iMultBin+1, task, bCorrelated);
    if(!hFlowTwoDif) { Error(Form("vn{2} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
    hFlowTwoDif->SetName(Form("%s_sample%d", nameFlowTwo.Data(), iSample));
    listFlowTwo->Add(hFlowTwoDif);

    if(bDoFour)
    {
      // loading reference <<2>>
      TProfile* pCorTwoRef = (TProfile*) listRefCorTwo->FindObject(Form("pCor2_Refs_harm%d_gap%s_sample%d_rebin",task->fHarmonics, task->GetEtaGapString().Data(), iSample));
      if(!pCorTwoRef) { Error(Form("Profile 'pCorTwoRef' (sample %d) does not exists",iSample),"ProcessDirect"); listRefCorTwo->ls(); return kFALSE; }

      // loading reference vn{4}
      TH1D* hFlowRefFour = (TH1D*) listRefFour->FindObject(Form("hFlow4_Refs_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
      if(!hFlowRefFour) { Error(Form("Histo 'hFlowRefFour' (sample %d) does not exists",iSample),"ProcessDirect"); listRefFour->ls(); return kFALSE; }

      // dn{4}
      TH1D* hCumFourDif = CalcDifCumFour(pCorFourDif, pCorTwoDif, pCorTwoRef, iMultBin+1, task, bCorrelated);
      if(!hCumFourDif) { Error(Form("dn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hCumFourDif->SetName(Form("%s_sample%d", nameCumFour.Data(), iSample));
      listCumFour->Add(hCumFourDif);

      // v'n{4}
      TH1D* hFlowFourDif = CalcDifFlowFour(hCumFourDif, hFlowRefFour, iMultBin+1, task, bCorrelated);
      if(!hFlowFourDif) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hFlowFourDif->SetName(Form("%s_sample%d", nameFlowFour.Data(), iSample));
      listFlowFour->Add(hFlowFourDif);
    }
  } // end-for {iSample} : loop over samples

  Debug("Merging correlations for central values", "ProcessDirect");

  // loading reference vn{2} merged for vn{2} dif. merged
  TH1D* hFlowTwoRefMerged = (TH1D*) ffOutputFile->Get(Form("hFlow2_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!hFlowTwoRefMerged) { Error(Form("Reference vn{2} (merged) not loaded!"),"ProcessDirect"); return kFALSE; }

  // <<2>>
  TProfile* pCorTwoMerged = (TProfile*) MergeListProfiles(listCorTwo);
  if(!pCorTwoMerged) { Error("Merging of 'pCorTwoMerged' failed!","ProcessDirect"); return kFALSE; }
  pCorTwoMerged->SetName(Form("%s_merged", nameCorTwo.Data()));

  // dn{2}
  TH1D* hCumTwoMerged = CalcDifCumTwo(pCorTwoMerged, task);
  if(!hCumTwoMerged) { Error(Form("dn{2} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
  hCumTwoMerged->SetName(Form("%s_merged", nameCumTwo.Data()));

  // vn{2}
  TH1D* hFlowTwoMerged = CalcDifFlowTwo(hCumTwoMerged, hFlowTwoRefMerged, iMultBin+1, task);
  if(!hFlowTwoMerged) { Error(Form("vn{2} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
  hFlowTwoMerged->SetName(Form("%s_merged", nameFlowTwo.Data()));

  TProfile* pCorFourMerged = 0x0;
  TH1D* hCumFourMerged = 0x0;
  TH1D* hFlowFourMerged = 0x0;
  if(bDoFour)
  {
    // loading reference <<2>> merged
    TProfile* pCorTwoRefMerged = (TProfile*) ffOutputFile->Get(Form("pCor2_Refs_harm%d_gap%s", task->fHarmonics, task->GetEtaGapString().Data() ));
    if(!pCorTwoRefMerged) { Error(Form("Reference <<2>> (merged) not loaded!"),"ProcessDirect"); ffOutputFile->ls(); return kFALSE; }

    // loading reference vn{4} merged
    TH1D* hFlowFourRefMerged = (TH1D*) ffOutputFile->Get(Form("hFlow4_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
    if(!hFlowFourRefMerged) { Error(Form("Reference vn{4} (merged) not loaded!"),"ProcessDirect"); ffOutputFile->ls(); return kFALSE; }

    // <<4>>
    pCorFourMerged = (TProfile*) MergeListProfiles(listCorFour);
    if(!pCorFourMerged) { Error("Merging of 'pCorFourMerged' failed!","ProcessDirect"); return kFALSE; }
    pCorFourMerged->SetName(Form("%s_merged", nameCorFour.Data()));

    // dn{4}
    hCumFourMerged = CalcDifCumFour(pCorFourMerged, pCorTwoMerged, pCorTwoRefMerged, iMultBin+1, task, bCorrelated);
    if(!hCumFourMerged) { Error(Form("cn{4} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
    hCumFourMerged->SetName(Form("%s_merged", nameCumFour.Data()));

    // vn{4}
    hFlowFourMerged = CalcDifFlowFour(hCumFourMerged, hFlowFourRefMerged, iMultBin+1, task, bCorrelated);
    if(!hFlowFourMerged) { Error(Form("vn{4} (merged) not processed correctly!"),"ProcessDirect"); return kFALSE; }
    hFlowFourMerged->SetName(Form("%s_merged", nameFlowFour.Data()));
  }

  Debug("Desampling","ProcessDirect");

  Debug(Form("<<2'>>: Number of samples in list pre-merging: %d",listCorTwo->GetEntries()),"ProcessDirect");
  TH1D* hDesampledTwo_Cor = DesampleList(listCorTwo, pCorTwoMerged->ProjectionX(), task, nameCorTwo, kTRUE); // skipping desampling for output structure
  if(!hDesampledTwo_Cor) { Error("Desampling <<2'>> failed","ProcessDirect"); return kFALSE; }
  hDesampledTwo_Cor->SetName(nameCorTwo.Data());

  Debug(Form("dn{2}: Number of samples in list pre-merging: %d",listCumTwo->GetEntries()),"ProcessDirect");
  TH1D* hDesampledTwo_Cum = DesampleList(listCumTwo, hCumTwoMerged, task, nameCumTwo);
  if(!hDesampledTwo_Cum) { Error("Desampling dn{2} failed","ProcessDirect"); return kFALSE; }
  hDesampledTwo_Cum->SetName(nameCumTwo.Data());

  Debug(Form("vn{2}: Number of samples in list pre-merging: %d",listFlowTwo->GetEntries()),"ProcessDirect");
  TH1D* hDesampledTwo = DesampleList(listFlowTwo, hFlowTwoMerged, task, nameFlowTwo);
  if(!hDesampledTwo) { Error("Desampling vn{2} failed","ProcessDirect"); return kFALSE; }
  hDesampledTwo->SetName(nameFlowTwo.Data());

  // saving to output file & cleaning
  ffOutputFile->cd();
  if(fSaveInterSteps) {
    hDesampledTwo_Cor->Write();
    hDesampledTwo_Cum->Write();
  }
  hDesampledTwo->Write();

  delete hDesampledTwo_Cum;
  delete hDesampledTwo;

  if(bDoFour)
  {
    Debug(Form("<<4'>>: Number of samples in list pre-merging: %d",listCorFour->GetEntries()),"ProcessDirect");
    TH1D* hDesampledFour_Cor = DesampleList(listCorFour, pCorFourMerged->ProjectionX(), task, nameCorFour, kTRUE); // skipping desampling for structure
    if(!hDesampledFour_Cor) { Error("Desampling dn{4} failed","ProcessDirect"); return kFALSE; }
    hDesampledFour_Cor->SetName(nameCorFour.Data());

    Debug(Form("dn{4}: Number of samples in list pre-merging: %d",listCumFour->GetEntries()),"ProcessDirect");
    TH1D* hDesampledFour_Cum = DesampleList(listCumFour, hCumFourMerged, task, nameCumFour);
    if(!hDesampledFour_Cum) { Error("Desampling dn{4} failed","ProcessDirect"); return kFALSE; }
    hDesampledFour_Cum->SetName(nameCumFour.Data());

    Debug(Form("vn{4}: Number of samples in list pre-merging: %d",listFlowFour->GetEntries()),"ProcessDirect");
    TH1D* hDesampledFour = DesampleList(listFlowFour, hFlowFourMerged, task, nameFlowFour);
    if(!hDesampledFour) { Error("Desampling vn{4} failed","ProcessDirect"); return kFALSE; }
    hDesampledFour->SetName(nameFlowFour.Data());

    // saving to output file & cleaning
    ffOutputFile->cd();
    if(fSaveInterSteps) {
      hDesampledFour_Cor->Write();
      hDesampledFour_Cum->Write();
    }
    hDesampledFour->Write();
    delete hDesampledFour_Cum;
    delete hDesampledFour;
  };

  delete listCorTwo;
  delete listCumTwo;
  delete listFlowTwo;

  delete listCorFour;
  delete listCumFour;
  delete listFlowFour;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ProcessReconstructed(FlowTask* task,Short_t iMultBin)
{
  Info("Processing task","ProcessReconstructed");
  if(!task) { Error("Task not valid!","ProcessReconstructed"); return kFALSE; }
  if(task->fNumPtBins < 1) { Error("Num of pt bins too low!","ProcessReconstructed"); return kFALSE; }

  if(fFlowFitCumulants) { fFlowFitCumulants = kFALSE; Warning("Fitting cumulants currently not available! WIP! switching flag off"); }

  // new naming convention for input histos (from FlowTask)
  TString sProfTwoName = Form("<<2>>(%d,-%d)",task->fHarmonics, task->fHarmonics);
  TString sProfFourName = Form("<<4>>(%d,%d,-%d,-%d)",task->fHarmonics, task->fHarmonics, task->fHarmonics, task->fHarmonics);
  if(task->fEtaGap > -1.0) {
    sProfTwoName += Form("_2sub(%.2g)",task->fEtaGap);
    sProfFourName += Form("_2sub(%.2g)",task->fEtaGap);
  }

  TString sSpeciesName = task->GetSpeciesName();
  TString sSpeciesLabel = task->GetSpeciesLabel();

  // ### Preparing slices of pt
  // if(!PrepareSlices(iMultBin,task,profFlow,histEntries,histEntriesBg,profFlowFour)) { return kFALSE; }

  if(!PrepareSlicesNew(task,sProfTwoName,kTRUE)) { Error(Form("PrepareSlicesNew '%s' failed!",sProfTwoName.Data()),"ProcessReconstructed"); return kFALSE; }
  if(task->fDoFour && !PrepareSlicesNew(task,sProfFourName,kFALSE)) { Error(Form("PrepareSlicesNew for '%s' failed!",sProfFourName.Data()),"ProcessReconstructed"); return kFALSE; }


  // things for later

  Bool_t bCorrelated = 0;

  // Loading list where reference flow samples are stored
  TList* listRefCorTwo = (TList*) ffDesampleFile->Get(Form("pCor2_Refs_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefCorTwo) { Error("List 'listRefCorTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  TList* listRefTwo = (TList*) ffDesampleFile->Get(Form("hFlow2_Refs_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
  if(!listRefTwo) { Error("List 'listRefTwo' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }

  TList* listRefFour = 0x0;
  if(task->fDoFour)
  {
    listRefFour = (TList*) ffDesampleFile->Get(Form("hFlow4_Refs_harm%d_gap%s_list",task->fHarmonics,task->GetEtaGapString().Data()));
    if(!listRefFour) { Error("List 'listRefFour' not found!","ProcessDirect"); ffDesampleFile->ls(); return kFALSE; }
  }


  Int_t iSample = 0;

  for(Int_t binMult(0); binMult < fiNumMultBins; ++binMult)
  {
    iMultBin = binMult;


    // List for desampling : later
    TString nameCorTwo = Form("pCor2_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameCumTwo = Form("hCum2_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameFlowTwo = Form("hFlow2_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);

    TString nameCorFour = Form("pCor4_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameCumFour = Form("hCum4_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    TString nameFlowFour = Form("hFlow4_%s_harm%d_gap%s_cent%d", task->GetSpeciesName().Data(), task->fHarmonics, task->GetEtaGapString().Data(),iMultBin);
    //

    TH1D* hCumTwoDif = 0x0;
    TH1D* hFlowTwoDif = 0x0;
    TH1D* hCumFourDif = 0x0;
    TH1D* hFlowFourDif = 0x0;

    Double_t dFlow = 0.0, dFlowError = 0.0; // containers for flow extraction results
    TH1D* hFlowMass = 0x0;
    TH1D* pCorTwoDif = 0x0;
    if(!fFlowFitCumulants) { pCorTwoDif = new TH1D(Form("pCor2_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin),Form("%s: <<2>>_{%d}{|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); <<2>>_{%d}{|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }
    else { pCorTwoDif = new TH1D(Form("hCum2_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin),Form("%s: d_{%d}{2,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); d_{%d}{2,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }


    Double_t dFlowFour = 0.0, dFlowFourError = 0.0; // containers for flow extraction results
    TH1D* hFlowMassFour = 0x0;
    TH1D* pCorFourDif = 0x0;
    if(!fFlowFitCumulants) { pCorFourDif = new TH1D(Form("pCor4_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin),Form("%s: <<4>>_{%d}{|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); <<4>>_{%d}{|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }
    else { pCorFourDif = new TH1D(Form("hCum4_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin),Form("%s: d_{%d}{4,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); d_{%d}{4,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap), task->fNumPtBins,task->fPtBinsEdges.data()); }


    for(Short_t binPt(0); binPt < task->fNumPtBins; binPt++)
    {
      TH1D* hInvMass = (TH1D*) task->fListHistos->FindObject(Form("hInvMass_mult%d_pt%d",binMult,binPt));
      if(!hInvMass) { Error("hInvMass histo not found among slices!","ProcessReconstructed"); task->fListHistos->ls(); return kFALSE; }
      hInvMass->SetTitle(Form("%s: InvMass dist (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
      hInvMass->SetMarkerStyle(kFullCircle);

      TH1D* hFlowMass = (TH1D*) task->fListProfiles->FindObject(Form("%s_Pos_sample0_mult%d_pt%d",sProfTwoName.Data(),binMult,binPt));
      if(!hFlowMass) { Error(Form("hFlowMass histo '%s' not found among slices!",sProfTwoName.Data()),"ProcessReconstructed"); task->fListProfiles->ls(); return kFALSE; }
      hFlowMass->SetTitle(Form("%s: FlowMass (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
      hFlowMass->SetMarkerStyle(kFullCircle);

      TList* listFits = new TList();
      listFits->SetOwner(1);

      TList* listFitsFour = 0x0;
      TF1 fitCorSigFour;
      TF1 fitCorBgFour;

      if(task->fDoFour)
      {
        hFlowMassFour = (TH1D*) task->fListProfiles->FindObject(Form("%s_Pos_sample0_mult%d_pt%d",sProfFourName.Data(),binMult,binPt));
        if(!hFlowMassFour) { Error(Form("hFlowMassFour histo '%s' not found among slices!",sProfFourName.Data()),"ProcessReconstructed"); task->fListProfiles->ls(); return kFALSE; }
        hFlowMassFour->SetTitle(Form("%s: FlowMassFour (|#Delta#eta| > %02.2g, cent %d, pt %d)",sSpeciesLabel.Data(),task->fEtaGap,iMultBin,binPt));
        hFlowMassFour->SetMarkerStyle(kFullCircle);
        listFitsFour = new TList();
        listFitsFour->SetOwner(1);
      }

      // ### Fitting the correlations

      TF1 fitOutSig;
      TF1 fitOutBg;

      Bool_t fitMass = FitInvMass(hInvMass, task, fitOutSig, fitOutBg);
      if(!fitMass) { Error("FitMass failed!","ProcessReconstructed"); ffFitsFile->cd(); fitOutSig.Write("fitMassSig"); fitOutBg.Write("fitMassBg"); return kFALSE; }

      TF1 fitCorSig;
      TF1 fitCorBg;

      Bool_t fitCor = FitCorrelations(hFlowMass, task, fitCorSig, fitCorBg, fitOutSig, fitOutBg);
      if(!fitCor) { Error("FitCor failed!","ProcessReconstructed"); ffFitsFile->cd(); fitCorSig.Write("fitCorSig"); fitCorBg.Write("fitCorBg"); return kFALSE; }

      // storing fits
      listFits->Add(hInvMass);
      listFits->Add(&fitOutSig);
      listFits->Add(&fitOutBg);
      listFits->Add(hFlowMass);
      listFits->Add(&fitCorSig);
      listFits->Add(&fitCorBg);

      Int_t iParFlow = fitCorSig.GetNpar() - 1;
      dFlow = fitCorSig.GetParameter(iParFlow);
      dFlowError = fitCorSig.GetParError(iParFlow);

      if(task->fDoFour)
      {
        Bool_t fitCorFour = FitCorrelations(hFlowMassFour, task, fitCorSigFour, fitCorBgFour, fitOutSig, fitOutBg);
        if(!fitCorFour) { Error("FitCorFour failed!","ProcessReconstructed"); return kFALSE; }

        // storing fits
        listFitsFour->Add(hInvMass);
        listFitsFour->Add(&fitOutSig);
        listFitsFour->Add(&fitOutBg);
        listFitsFour->Add(hFlowMassFour);
        listFitsFour->Add(&fitCorSigFour);
        listFitsFour->Add(&fitCorBgFour);

        Int_t iParFlow = fitCorSigFour.GetNpar() - 1;
        dFlowFour = fitCorSigFour.GetParameter(iParFlow);
        dFlowFourError = fitCorSigFour.GetParError(iParFlow);
      }

      // setting the flow
      pCorTwoDif->SetBinContent(binPt+1,dFlow);
      pCorTwoDif->SetBinError(binPt+1,dFlowError);

      Double_t dFlowRel = -999.9; if(TMath::Abs(dFlow) > 0.0) { dFlowRel = dFlowError / dFlow; }
      Info(Form("Final vn{2}: (mult %d | pt %d) %g +- %g (rel. %.3f)",iMultBin,binPt,dFlow,dFlowError,dFlowRel), "ProcessReconstructed");

      if(task->fDoFour)
      {
        pCorFourDif->SetBinContent(binPt+1,dFlowFour);
        pCorFourDif->SetBinError(binPt+1,dFlowFourError);

        Double_t dFlowRel = -999.9; if(TMath::Abs(dFlowFour) > 0.0) { dFlowRel = dFlowFourError / dFlowFour; }
        Info(Form("Final vn{4}: (mult %d | pt %d) %g +- %g (rel. %.3f)",iMultBin,binPt,dFlowFour,dFlowFourError,dFlowRel), "ProcessReconstructed");
      }

      // processing / plotting fits
      ffFitsFile->cd();
      listFits->Write(Form("fits_%s_cent%d_pt%d",sSpeciesName.Data(),iMultBin,binPt),TObject::kSingleKey);
      if(task->fDoFour) { listFitsFour->Write(Form("fitsFour_%s_cent%d_pt%d",sSpeciesName.Data(),iMultBin,binPt),TObject::kSingleKey); }

      // HERE THE CORRELATIONS (NOT FLOW) are ready!!!
      // TODO: Rename that stuff!!!

    } // endfor {binPt}

    pCorTwoDif->SetName(Form("pCor2_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));

    Debug("Calculating flow","ProcessDirect");
    // loading reference vn{2}
    TH1D* hFlowRefTwo = (TH1D*) listRefTwo->FindObject(Form("hFlow2_Refs_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
    if(!hFlowRefTwo) { Error(Form("Histo 'hFlowRefTwo' (sample %d) does not exists",iSample),"ProcessReconstructed"); listRefTwo->ls(); ffDesampleFile->ls(); return kFALSE; }

    // dn{2}
    hCumTwoDif = CalcDifCumTwo(pCorTwoDif, task);
    if(!hCumTwoDif) { Error(Form("dn{2} (sample %d) not processed correctly!",iSample),"ProcessReconstructed"); return kFALSE; }
    hCumTwoDif->SetName(Form("%s", nameCumTwo.Data()));

    // v'n{2}
    hFlowTwoDif = CalcDifFlowTwo(hCumTwoDif, hFlowRefTwo, iMultBin+1, task, bCorrelated);
    if(!hFlowTwoDif) { Error(Form("vn{2} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
    hFlowTwoDif->SetName(Form("%s", nameFlowTwo.Data()));

    if(task->fDoFour)
    {
      pCorFourDif->SetName(Form("pCor4_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
      // loading reference <<2>>
      TProfile* pCorTwoRef = (TProfile*) listRefCorTwo->FindObject(Form("pCor2_Refs_harm%d_gap%s_sample%d_rebin",task->fHarmonics, task->GetEtaGapString().Data(), iSample));
      if(!pCorTwoRef) { Error(Form("Profile 'pCorTwoRef' (sample %d) does not exists",iSample),"ProcessDirect"); listRefCorTwo->ls(); return kFALSE; }

      // loading reference vn{4}
      TH1D* hFlowRefFour = (TH1D*) listRefFour->FindObject(Form("hFlow4_Refs_harm%d_gap%s_sample%d",task->fHarmonics,task->GetEtaGapString().Data(),iSample));
      if(!hFlowRefFour) { Error(Form("Histo 'hFlowRefFour' (sample %d) does not exists",iSample),"ProcessDirect"); listRefFour->ls(); return kFALSE; }

      // dn{4}
      hCumFourDif = CalcDifCumFour(pCorFourDif, pCorTwoDif, pCorTwoRef, iMultBin+1, task, bCorrelated);
      if(!hCumFourDif) { Error(Form("dn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hCumFourDif->SetName(Form("%s", nameCumFour.Data()));

      // v'n{4}
      hFlowFourDif = CalcDifFlowFour(hCumFourDif, hFlowRefFour, iMultBin+1, task, bCorrelated);
      if(!hFlowFourDif) { Error(Form("vn{4} (sample %d) not processed correctly!",iSample),"ProcessDirect"); return kFALSE; }
      hFlowFourDif->SetName(Form("%s", nameFlowFour.Data()));
    }


    ffOutputFile->cd();

    if(fSaveInterSteps) {
      pCorTwoDif->Write();
      hCumTwoDif->Write();
    }

    hFlowTwoDif->Write();

    if(task->fDoFour) {
      if(fSaveInterSteps) {
        pCorFourDif->Write();
        hCumFourDif->Write();
      }
      hFlowFourDif->Write();
    }

    // if(fFlowFitCumulants)
    // {
    //   TH1D* hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow2_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
    //   if(!hRefFlow) { Error("Something went wrong when running automatic refs flow task:","ProcessReconstructed"); return kFALSE; }
    //
    //   TH1D* hFlow_vn = CalcDifFlowTwo(hFlow, hRefFlow, iMultBin+1 ,task, task->fConsCorr);
    //   hFlow_vn->SetName(Form("hFlow2_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
    //   hFlow_vn->SetTitle(Form("%s: v_{%d}{2,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); v_{%d}{2,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap));
    //   hFlow_vn->Write();
    //
    //   if(task->fDoFour)
    //   {
    //     TH1D* hRefFlow = (TH1D*) ffOutputFile->Get(Form("hFlow4_Refs_harm%d_gap%s",task->fHarmonics,task->GetEtaGapString().Data()));
    //     if(!hRefFlow) { Error("Something went wrong when running automatic refs flow task:","ProcessReconstructed"); return kFALSE; }
    //     TH1D* hFlow_vn = CalcDifFlowFour(hFlowFour, hRefFlow, iMultBin+1, task, task->fConsCorr);
    //     hFlow_vn->SetName(Form("hFlow4_%s_harm%d_gap%s_cent%d",sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin));
    //     hFlow_vn->SetTitle(Form("%s: v_{%d}{4,|#Delta#eta|>%g} (%g - %g); #it{p}_{T} (GeV/#it{c}); v_{%d}{4,|#Delta#eta|>%g}",sSpeciesLabel.Data(),task->fHarmonics,task->fEtaGap,fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fHarmonics,task->fEtaGap));
    //     hFlow_vn->Write();
    //   }
    // }

    TCanvas* cFlow = new TCanvas("cFlow","cFlow");
    cFlow->cd();
    hFlowTwoDif->SetStats(0);
    hFlowTwoDif->Draw();
    cFlow->SaveAs(Form("%s/Flow_%s_n%d2_gap%s_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());

    if(task->fDoFour)
    {
      TCanvas* cFlow = new TCanvas("cFlow","cFlow");
      cFlow->cd();
      hFlowFourDif->SetStats(0);
      hFlowFourDif->Draw();
      cFlow->SaveAs(Form("%s/FlowFour_%s_n%d2_gap%s_cent%d.%s",fsOutputFilePath.Data(),sSpeciesName.Data(),task->fHarmonics,task->GetEtaGapString().Data(),iMultBin,fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
    }
  }


  return kTRUE;
}
//_____________________________________________________________________________
TH1* ProcessUniFlow::MergeListProfiles(TList* list)
{
  // merge list of TProfiles into single TProfile and return it
  if(!list || list->IsEmpty()) { Error("List not valid or empty","MergeListProfiles"); return 0x0; }

  TH1* merged = (TH1*) list->At(0)->Clone();
  // merged->SetName(Form("%s_merged",merged->GetName()));

  if(list->GetEntries() < 2) // only 1 entry
  {
    Warning("Only one entry for merging; returning it directly instead!","MergeListProfiles");
    return merged;
  }

  merged->Reset();
  Double_t mergeStatus = merged->Merge(list);
  if(mergeStatus == -1) { Error("Merging failed!","MergeListProfiles"); return 0x0; }

  return merged;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::DesampleList(TList* list, TH1D* merged, FlowTask* task, TString name, Bool_t bSkipDesampling)
{
  // bSkipDesampling [kFALSE] is used to keep the output structure cosistent (for check) even if Desampling can be skipped
  // e.g. vn{2} ~ cn{2} = <<2>>

  if(!merged) { Error("Merged histogram not valid","DesampleList"); return 0x0; }
  if(!list) { Error("List does not valid","DesampleList"); return 0x0; }
  if(list->GetEntries() < 1) { Error("List is empty","DesampleList"); return 0x0; }
  if(list->GetEntries() != task->fNumSamples) { Warning("Number of list entries is different from task number of samples","DesampleList"); }
  if(!task) { Error("FlowTask does not exists","DesampleList"); return 0x0; }

  Debug(Form("Number of samples in list pre-desampling: %d",list->GetEntries()),"DesampleList");

  TH1D* hDesampled = (TH1D*) merged->Clone(Form("%s_desampled",name.Data()));
  if(!hDesampled) { Error("Histo 'hDesampled' cloning failed","DesampleList"); return 0x0; }

  // saving objects to separate desampling file
  ffDesampleFile->cd();
  list->SetName(Form("%s_list",name.Data()));
  list->Write(0,TObject::kSingleKey);
  merged->SetName(Form("%s_merged",name.Data()));
  merged->Write();

  if(bSkipDesampling || task->fNumSamples < 2 || list->GetEntries() < 2)  // only one sample -> no sampling needed
  {
    // Warning("Only 1 sample for desampling; returning merged instead!","DesampleList");
    ffDesampleFile->cd();
    hDesampled->Write();
    return hDesampled;
  }

  for(Int_t iBin(0); iBin < hDesampled->GetNbinsX()+2; ++iBin)
  {
    const Double_t dDesMean = hDesampled->GetBinContent(iBin);

    Double_t dSum = 0.0;
    Int_t iCount = 0;

    for(Short_t iSample(0); iSample < list->GetEntries(); ++iSample)
    {
      TH1D* hTemp = (TH1D*) list->At(iSample);
      if(!hTemp) { Error(Form("Histo 'hTemp' (bin %d, sample %d) not found in list",iBin,iSample),"DesampleList"); return 0x0; }

      Double_t dContent = hTemp->GetBinContent(iBin);
      Double_t dError = hTemp->GetBinError(iBin);

      // TODO Check content & error

      dSum += TMath::Power(dDesMean - dContent, 2.0);
      iCount++;
    } // end-for {iSample} : loop over samples in list

    Double_t dError = TMath::Sqrt(dSum / (iCount*iCount));
    hDesampled->SetBinError(iBin,dError);

  } // end-for {iBin} : loop over bins in histo

  PlotDesamplingQA(list, hDesampled, task);

  hDesampled->Write();
  return hDesampled;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PlotDesamplingQA(TList* list, TH1D* hDesampled, FlowTask* task)
{
  if(!list) { Error("Input samples list not found!","PlotDesamplingQA"); return kFALSE; }
  if(!hDesampled) { Error("Desampled result not found!","PlotDesamplingQA"); return kFALSE; }
  if(!task) { Error("FlowTask not found!","PlotDesamplingQA"); return kFALSE; }

  Warning("Not implemented yet!","PlotDesamplingQA");
    //
    //
    // // getting copy which does not affect histo which is returned
    // TH1D* hDesampledClone = (TH1D*) hDesampled->Clone(Form("%sClone",hDesampled->GetName()));
    //
    // TList* listOutput = new TList(); // list for collecting all QA histos
    //
    // // doing QA plots with spread, etc.
    // TCanvas* canDesample = new TCanvas(Form("canDesample_%s",task->fName.Data()),Form("canDesample_%s",task->fName.Data()),1200,400);
    // canDesample->Divide(3,1);
    //
    // TH1D* hTempRatio = 0x0;
    // TH1D* hTempError = 0x0;
    //
    // TLine* lineUnity = new TLine();
    // lineUnity->SetLineColor(kRed);
    // lineUnity->SetLineWidth(3);
    //
    // canDesample->cd(1);
    // hDesampledClone->SetStats(kFALSE);
    // hDesampledClone->SetFillColor(kBlue);
    // hDesampledClone->SetStats(kFALSE);
    // hDesampledClone->SetMarkerStyle(20);
    // hDesampledClone->SetMarkerSize(0.5);
    // hDesampledClone->SetMarkerColor(kRed);
    // hDesampledClone->DrawCopy("E2");
    //
    // for(Short_t iSample(0); iSample < task->fNumSamples; iSample++)
    // {
    //   hTempSample = (TH1D*) list->At(iSample);
    //   if(!hTempSample) { Warning(Form("Sample %d not found during plotting QA! Skipping!",iSample),"DesampleList"); continue; }
    //
    //   canDesample->cd(1);
    //   hTempSample->SetStats(kFALSE);
    //   hTempSample->SetLineColor(30+2*iSample);
    //   hTempSample->SetMarkerColor(30+2*iSample);
    //   hTempSample->SetMarkerStyle(24);
    //   hTempSample->SetMarkerSize(0.5);
    //   hTempSample->DrawCopy("hist p same");
    //
    //   hTempRatio = (TH1D*) hTempSample->Clone(Form("%s_ratio",hTempSample->GetName()));
    //   hTempRatio->Divide(hDesampled);
    //   hTempRatio->SetYTitle("Value: final / sample");
    //   hTempRatio->SetTitleOffset(1.2,"Y");
    //
    //   canDesample->cd(2);
    //   hTempRatio->SetMinimum(0.6);
    //   hTempRatio->SetMaximum(1.4);
    //   hTempRatio->Draw("hist p same");
    //
    //   hTempError = (TH1D*) hTempSample->Clone(Form("%s_error",hTempSample->GetName()));
    //   for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hTempSample->GetBinError(bin)); }
    //
    //   canDesample->cd(3);
    //   hTempError->SetMinimum(0.);
    //   hTempError->SetMaximum(1.5*hTempError->GetMaximum());
    //   hTempError->SetYTitle("Uncertainty");
    //   hTempError->SetTitleOffset(1.2,"Y");
    //
    //   hTempError->Draw("hist p same");
    //
    //   listOutput->Add(hTempSample);
    //   listOutput->Add(hTempRatio);
    //   listOutput->Add(hTempError);
    // }
    //
    // canDesample->cd(1);
    // hDesampledClone->DrawCopy("hist p same");
    //
    // canDesample->cd(2);
    // lineUnity->DrawLine(hTempRatio->GetXaxis()->GetXmin(),1,hTempRatio->GetXaxis()->GetXmax(),1);
    //
    // hTempError = (TH1D*) hDesampledClone->Clone(Form("%s_error",hDesampled->GetName()));
    // for(Short_t bin(1); bin < hTempSample->GetNbinsX()+1; bin++) { hTempError->SetBinContent(bin,hDesampledClone->GetBinError(bin)); }
    // listOutput->Add(hTempError);
    //
    // canDesample->cd(3);
    // hTempError->Draw("hist p same");
    //
    // // saving QA plots
    // canDesample->SaveAs(Form("%s/Desampling_%s_harm%d_gap%g_cent%d_%s.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),task->fHarmonics,10*task->fEtaGap,iMultBin,task->fName.Data(),fsOutputFileFormat.Data()));
    //
    // Info("Saving desampling QA into output file","DesampleList");
    // ffDesampleFile->cd();
    // listOutput->Add(canDesample);
    // listOutput->Write(Form("Desampling_%s_cent%d_%s",task->GetSpeciesName().Data(),iMultBin,task->fName.Data()),TObject::kSingleKey);
    //
    // // deleting created stuff
    // delete listOutput;
    // // delete canDesample;
    // delete lineUnity;
    // // if(hTempSample) delete hTempSample;
    // // delete hTempRatio;
    // // delete hTempError;
    // // delete hDesampledClone;
    //
    //

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PrepareSlices(const Short_t multBin, FlowTask* task, TProfile3D* p3Cor, TH3D* h3Entries, TH3D* h3EntriesBG, TProfile3D* p3CorFour)
{

  if(!task) { Error("Input task not found!","PrepareSlices"); return kFALSE; }
  if(!h3Entries) { Error("Input hist with entries not found!","PrepareSlices"); return kFALSE; }
  if(!p3Cor) { Error("Input profile with correlations not found!","PrepareSlices"); return kFALSE; }
  if(task->fDoFour && !p3CorFour) { Error("Input profile with <<4'>> not found!","PrepareSlices"); return kFALSE; }
  if(multBin < 0 || multBin > fiNumMultBins) { Error("Wrong multiplicity bin index (not in range)!","PrepareSlices"); return kFALSE; }

  // cleaning the vectros with flow-mass and inv. mass plots
  if(task->fVecHistInvMass->size() > 0) { task->fVecHistInvMass->clear(); }
  if(task->fVecHistInvMassBG->size() > 0) { task->fVecHistInvMassBG->clear(); }
  if(task->fVecHistFlowMass->size() > 0) { task->fVecHistFlowMass->clear(); }
  if(task->fDoFour && task->fVecHistFlowMassFour->size() > 0) { task->fVecHistFlowMassFour->clear(); }

  const Short_t binMultLow = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin]);
  const Short_t binMultHigh = h3Entries->GetXaxis()->FindFixBin(fdMultBins[multBin+1]) - 1;
  // printf("Mult: %g(%d) -  %g(%d)\n",fdMultBins[multBin],binMultLow,fdMultBins[multBin+1],binMultHigh);

  TH1D* hRefFlow = 0x0;
  if(!fFlowFitCumulants)
  {
    // loading reference flow, if not found, it will be prepared
    hRefFlow = (TH1D*) ffDesampleFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g_desampled",task->fHarmonics,10*task->fEtaGap));
    if(!hRefFlow)
    {
      Warning("Relevant Reference flow not found within output file.","PrepareSlices");
      ffOutputFile->ls();

      Info("Creating relevant reference flow task.","PrepareSlices");
      FlowTask* taskRef = new FlowTask(FlowTask::kRefs,"Ref");
      taskRef->SetHarmonics(task->fHarmonics);
      taskRef->SetEtaGap(task->fEtaGap);
      taskRef->SetNumSamples(task->fNumSamples);
      taskRef->SetInputTag(task->fInputTag);
      taskRef->SetDoFourCorrelations(task->fDoFour);

      if(ProcessRefs(taskRef))
      {
        hRefFlow = (TH1D*) ffDesampleFile->Get(Form("hFlow2_Refs_harm%d_gap%02.2g",task->fHarmonics,10*task->fEtaGap));
        if(!hRefFlow) {  Error("Automated Refs task completed, but RefFlow not found!","PrepareSlices"); ffDesampleFile->ls(); return kFALSE; }
      }
      else { Error("Something went wrong when running automatic refs flow task:","PrepareSlices"); taskRef->PrintTask(); return kFALSE; }
    }
  }

  TH1D* hRefCor2 = 0x0;
  TH1D* hRefFlow4 = 0x0;
  if(task->fDoFour)
  {
    // loading <<2>> refs for dn{4}
    const char* name = Form("pCor2_Refs_harm%d_gap%02.2g_desampled",task->fHarmonics,10*task->fEtaGap);
    hRefCor2 = (TH1D*) ffDesampleFile->Get(name);
    if(!hRefCor2) { Error(Form("Input '%s' (hRefCor2) not found!",name),"PrepareSlices"); ffDesampleFile->ls(); return kFALSE; }

    const char* nameFour = Form("hFlow4_Refs_harm%d_gap%02.2g_desampled",task->fHarmonics,10*task->fEtaGap);
    hRefFlow4 = (TH1D*) ffDesampleFile->Get(nameFour);
    if(!hRefFlow4) { Error(Form("Input '%s' (hRefFlow4) not found!",nameFour),"PrepareSlices"); ffDesampleFile->ls(); return kFALSE; }
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
  // prof2FlowMass_temp = (TProfile2D*) prof3Flow_temp->Project3DProfile("yz"); // NOTE: standard ROOT way - working properly in ROOTv6-12 onwards
  prof2FlowMass_temp = Project3DProfile(prof3Flow_temp);

  TProfile3D* prof3FlowFour_temp = 0x0;
  TProfile2D* prof2FlowMassFour_temp = 0x0;
  TProfile* profFlowMassFour_temp = 0x0;
  TH1D* hFlowMassFour_temp = 0x0;
  if(task->fDoFour)
  {
    prof3FlowFour_temp = (TProfile3D*) p3CorFour->Clone(Form("prof3FlowFour_temp_cent%d",multBin));
    prof3FlowFour_temp->GetXaxis()->SetRange(binMultLow,binMultHigh);
    // prof2FlowMassFour_temp = (TProfile2D*) prof3FlowFour_temp->Project3DProfile("yz"); // NOTE: standard ROOT way - working properly in ROOTv6-12 onwards
    prof2FlowMassFour_temp = Project3DProfile(prof3FlowFour_temp);
  }

  Short_t iNumPtBins = task->fNumPtBins;

  TCanvas* canInvMass = new TCanvas("canInvMass","InvMass",1400,600);
  TCanvas* canInvMassBG = new TCanvas("canInvMassBG","InvMassBG",1400,600);
  TCanvas* canFlowMass = new TCanvas("canFlowMass","FlowMass",1400,600);
  TCanvas* canFlowMassFour = new TCanvas("canFlowMassFour","FlowMassFour",1400,600);
  canInvMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canInvMassBG->Divide(5,std::ceil(iNumPtBins/5)+1);
  canFlowMass->Divide(5,std::ceil(iNumPtBins/5)+1);
  canFlowMassFour->Divide(5,std::ceil(iNumPtBins/5)+1);

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

    task->fVecHistInvMass->push_back(hInvMass_temp);
    if(h3EntriesBG) { task->fVecHistInvMassBG->push_back(hInvMassBG_temp); }
    // done with InvMass

    // ### projection of flow-mass profile
    profFlowMass_temp = (TProfile*) prof2FlowMass_temp->ProfileX(Form("profFlowMass_cent%d_pt%d",multBin,binPt),binPtLow,binPtHigh);
    if(task->fDoFour) { profFlowMassFour_temp = (TProfile*) prof2FlowMassFour_temp->ProfileX(Form("profFlowMassFour_cent%d_pt%d",multBin,binPt),binPtLow,binPtHigh); }

    // checking for rebinning the flow-mass profile
    if(task->fRebinFlowMass > 1) {
      profFlowMass_temp->Rebin(task->fRebinFlowMass);
      if(task->fDoFour) { profFlowMassFour_temp->Rebin(task->fRebinFlowMass); }
    }

    // prepare slices for <<2>>
    hFlowMass_temp = CalcDifCumTwo(profFlowMass_temp,task);
    if(!hFlowMass_temp) { Error("<<2'>> not ready!","PrepareSlices"); return kFALSE; }
    hFlowMass_temp->SetName(Form("hFlowMass_cent%d_pt%d",multBin,binPt));
    if(!fFlowFitCumulants) { hFlowMass_temp = CalcDifFlowTwo(hFlowMass_temp, hRefFlow, multBin+1, task, task->fConsCorr); }
    if(!hFlowMass_temp) { Error("hFlowMass_temp not ready! Something went wrong!","PrepareSlices"); return kFALSE; }
    task->fVecHistFlowMass->push_back(hFlowMass_temp);

    // prepare slices for <<4>>
    if(task->fDoFour) {
      hFlowMassFour_temp = CalcDifCumFour(profFlowMassFour_temp, profFlowMass_temp, hRefCor2, multBin+1, task, task->fConsCorr);
      if(!hFlowMassFour_temp) { Error("<<4'>> not ready!","PrepareSlices"); return kFALSE; }
      hFlowMassFour_temp->SetName(Form("hFlowMassFour_cent%d_pt%d",multBin,binPt));
      if(!fFlowFitCumulants) { hFlowMassFour_temp = CalcDifFlowFour(hFlowMassFour_temp, hRefFlow4, multBin+1, task, task->fConsCorr); }
      if(!hFlowMassFour_temp) { Error("FlowMassFour_temp not ready! Something went wrong!","PrepareSlices"); return kFALSE; }
      task->fVecHistFlowMassFour->push_back(hFlowMassFour_temp);
    }

    canInvMass->cd(binPt+1);
    hInvMass_temp->Draw();

    if(h3EntriesBG)
    {
      canInvMassBG->cd(binPt+1);
      hInvMassBG_temp->Draw();
    }

    canFlowMass->cd(binPt+1);
    hFlowMass_temp->Draw();

    if(task->fDoFour)
    {
      canFlowMassFour->cd(binPt+1);
      hFlowMassFour_temp->Draw();
    }

  } // endfor {binPt}: over Pt bins

  printf(" # of slices: InvMass: %lu | InvMassBG %lu | FlowMass %lu\n",task->fVecHistInvMass->size(),task->fVecHistInvMassBG->size(),task->fVecHistFlowMass->size());

  gSystem->mkdir(Form("%s/slices/",fsOutputFilePath.Data()));
  canInvMass->SaveAs(Form("%s/slices/Slices_InvMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  if(h3EntriesBG) canInvMassBG->SaveAs(Form("%s/slices/Slices_InvMassBG_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  canFlowMass->SaveAs(Form("%s/slices/Slices_FlowMass_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data()));
  if(task->fDoFour) { canFlowMassFour->SaveAs(Form("%s/slices/Slices_FlowMassFour_%s_gap%g_cent%d.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),10*task->fEtaGap,multBin,fsOutputFileFormat.Data())); }

  if(task->fVecHistInvMass->size() < 1 || task->fVecHistFlowMass->size() < 1 || task->fVecHistFlowMass->size() != task->fVecHistInvMass->size()) { Error("Output vector empty. Something went wrong","PrepareSlices"); return kFALSE; }
  if(h3EntriesBG && (task->fVecHistInvMassBG->size() < 1 || task->fVecHistInvMassBG->size() != task->fVecHistInvMass->size()) ) { Error("Output vector empty. Something went wrong with BG histograms","PrepareSlices"); return kFALSE; }
  if(task->fDoFour && (task->fVecHistFlowMassFour->size() < 1 || task->fVecHistFlowMass->size() != task->fVecHistInvMass->size()) ) { Error("Output vector for <<4>> empty. Something went wrong","PrepareSlices"); return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::PrepareSlicesNew(FlowTask* task, TString histName, Bool_t bDoCand)
{
  // wrapper for making/preparing per-task slices
  if(!task) { Error("FlowTask does not exists!","PrepareSlicesNew"); return kFALSE; }

  FlowTask::PartSpecies species = task->fSpecies;
  Bool_t bReco = kFALSE;
  if(species == FlowTask::kK0s || species == FlowTask::kLambda || species == FlowTask::kPhi) { bReco = kTRUE; }

  TList* inputList = 0x0;
  switch (species) {
    case FlowTask::kCharged : inputList = flFlowCharged; break;
    case FlowTask::kPion : inputList = flFlowPion; break;
    case FlowTask::kKaon : inputList = flFlowKaon; break;
    case FlowTask::kProton : inputList = flFlowProton; break;
    case FlowTask::kK0s : inputList = flFlowK0s; break;
    case FlowTask::kLambda : inputList = flFlowLambda; break;
    case FlowTask::kPhi : inputList = flFlowPhi; break;
    default: Error("Invalid species!","PrepareSlicesNew"); return kFALSE;
  }

  // preparing flow slices
  TH1* prof = 0x0;
  if(task->fMergePosNeg)
  {
    TH1* profPos = (TH1*) inputList->FindObject(Form("%s_Pos_sample0",histName.Data()));
    if(!profPos) { Error("Positive profile 'profNeg' not found!","PrepareSlicesNew"); inputList->ls(); return kFALSE; }
    TH1* profNeg = (TH1*) inputList->FindObject(Form("%s_Neg_sample0",histName.Data()));
    if(!profNeg) { Error("Negative profile 'profNeg' not found!","PrepareSlicesNew"); inputList->ls(); return kFALSE; }

    TList* listMerge = new TList();
    listMerge->Add(profPos);
    listMerge->Add(profNeg);
    prof = (TH1*) MergeListProfiles(listMerge);
    delete listMerge;
  } else {
    prof = (TH1*) inputList->FindObject(Form("%s_Pos_sample0",histName.Data()));
  }

  if(!prof) { Error(Form("Profile '%s' not found!",histName.Data()),"PrepareSlicesNew"); inputList->ls(); return kFALSE; }
  if(!MakeProfileSlices(task,prof,task->fListProfiles)) { Error("Profile Slices failed!","PrepareSlicesNew"); return kFALSE; };

  // preparing inv. mass slices (NB: merging pos/neg done in MakeSparseSlices() )
  if(bReco && bDoCand) {

    TString sNameCand = Form("fhsCand%s",task->GetSpeciesName().Data());

    THnSparseD* sparse = (THnSparseD*) inputList->FindObject(sNameCand.Data());
    if(!MakeSparseSlices(task,sparse,task->fListHistos)) { Error("Histo Slices failed!","PrepareSlicesNew"); return kFALSE; };

    if(species == FlowTask::kPhi) {
      THnSparseD* sparseBg = (THnSparseD*) inputList->FindObject(Form("%sBg",sNameCand.Data()));
      if(!MakeSparseSlices(task,sparseBg,task->fListHistos,"hInvMassBg")) { Error("Histo Slices for Phi BG failed!","PrepareSlicesNew"); return kFALSE; };
    }
  }

  if(fbSaveSlices && bReco) {
    // due to naming convention
    if(!task->fProcessMixed) { Warning("SaveSlices() implemented for Mixed harmonics only so far! Skipping!","PrepareSlicesNew"); return kTRUE; }
    if(!SaveSlices(task)) { Error("Save Slices failed!","PrepareSlicesNew"); return kFALSE; }
  }


  return kTRUE;
}

//_____________________________________________________________________________
Bool_t ProcessUniFlow::MakeProfileSlices(FlowTask* task, TH1* inputProf, TList* outList)
{
  // prepare slices out of inputHist
  if(!task) { Error("FlowTask does not exists!","MakeProfileSlices"); return kFALSE; }
  if(!inputProf) { Error("Input profile does not exists!","MakeProfileSlices"); return kFALSE; }
  if(!outList) { Error("Output TList does not exists!","MakeProfileSlices"); return kFALSE; }

  FlowTask::PartSpecies spec = task->fSpecies;
  if(spec == FlowTask::kRefs) { Error("Species is 'kRefs': no slicing required!","MakeProfileSlices"); return kFALSE; }

  Bool_t bReco = kFALSE;
  if(spec == FlowTask::kK0s || spec == FlowTask::kLambda || spec == FlowTask::kPhi) { bReco = kTRUE; }

  Int_t iNumBinsMult = fiNumMultBins;
  Int_t iNumBinsPt = task->fNumPtBins;

  TAxis* axisMult = inputProf->GetXaxis();
  TAxis* axisPt = inputProf->GetYaxis();

  TList trashCol;
  trashCol.SetOwner(kTRUE);

  for(Int_t iBinMult(0); iBinMult < iNumBinsMult; ++iBinMult) {
    const Int_t iBinMultLow = axisMult->FindFixBin(fdMultBins[iBinMult]);
    const Int_t iBinMultHigh = axisMult->FindFixBin(fdMultBins[iBinMult+1]) - 1;

    if(!bReco) {
      // direct species
      TProfile* prof1D_preRebin = ((TProfile2D*)inputProf)->ProfileY("",iBinMultLow,iBinMultHigh);
      if(!prof1D_preRebin) { Error("Profile 'prof1D_preRebin' failed!","MakeProfileSlices"); return kFALSE; }
      trashCol.Add(prof1D_preRebin); // to ensure to be deleted

      TProfile* prof1D = 0x0;
      if(iNumBinsPt > 0) { prof1D = (TProfile*) prof1D_preRebin->Rebin(iNumBinsPt,Form("%s_rebin", inputProf->GetName()), task->fPtBinsEdges.data()); }
      else { prof1D = (TProfile*) prof1D_preRebin->Clone(); }
      if(!prof1D) { Error("Profile 'prof1D' does not exists!","MakeProfileSlices"); return kFALSE; }

      prof1D->SetName(Form("%s_mult%d",inputProf->GetName(),iBinMult));
      prof1D->GetXaxis()->SetTitle(inputProf->GetYaxis()->GetTitle());
      outList->Add(prof1D);
    } else {
      // reconstructed species
      axisMult->SetRange(iBinMultLow,iBinMultHigh);
      TProfile2D* prof2D = Project3DProfile((TProfile3D*) inputProf);
      if(!prof2D) { Error("Mult projection failed!","MakeProfileSlices"); return kFALSE; }
      trashCol.Add(prof2D); // NB: to ensure that it will be deleted

      for(Int_t iBinPt(0); iBinPt < iNumBinsPt; ++iBinPt) {
        const Double_t dEdgePtLow = task->fPtBinsEdges[iBinPt];
        const Double_t dEdgePtHigh = task->fPtBinsEdges[iBinPt+1];

        const Int_t iBinPtLow = axisPt->FindFixBin(dEdgePtLow);
        const Int_t iBinPtHigh = axisPt->FindFixBin(dEdgePtHigh) - 1;

        TProfile* prof1D = prof2D->ProfileX("",iBinPtLow,iBinPtHigh);
        if(!prof1D) { Error("Profile 'prof1D' failed!","MakeProfileSlices"); return kFALSE; }

        if(task->fRebinFlowMass > 1) {
          prof1D->Rebin(task->fRebinFlowMass);
        }

        prof1D->SetName(Form("%s_mult%d_pt%d",inputProf->GetName(),iBinMult,iBinPt));
        prof1D->GetXaxis()->SetTitle(inputProf->GetZaxis()->GetTitle());
        outList->Add(prof1D);
      } // end-for {binPt}
    } // end-else {!bReco}
  } // end-for {binMult}

  Info("Successfull!","MakeProfileSlices");
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::MakeSparseSlices(FlowTask* task, THnSparse* inputSparse, TList* outList, const char* outName)
{
  // prepare slices out of inputHist
  if(!task) { Error("FlowTask does not exists!","MakeSparseSlices"); return kFALSE; }
  if(!inputSparse) { Error("Input THnSparse does not exists!","MakeSparseSlices"); return kFALSE; }
  if(!outList) { Error("Output TList does not exists!","MakeSparseSlices"); return kFALSE; }
  // if(outList->GetEntries() > 0) { Error("Output TList is not empty!","MakeSparseSlices"); return kFALSE; }

  FlowTask::PartSpecies species = task->fSpecies;
  if(species != FlowTask::kK0s && species != FlowTask::kLambda && species != FlowTask::kPhi) {
    Error("Not a reconstructed species!","MakeSparseSlices");
    return kFALSE;
  }

  TList trashCol;
  trashCol.SetOwner(kTRUE);

  Double_t dEtaGap = task->fEtaGap;
  Bool_t bHasGap = dEtaGap < 0.0 ? kFALSE : kTRUE;

  TH3D* histEntries = 0x0;
  TH3D* histEntriesPos = 0x0;
  TH3D* histEntriesNeg = 0x0;

  if(!bHasGap) {
    histEntries = (TH3D*) inputSparse->Projection(1,2,0);
    trashCol.Add(histEntries);
  } // end-if {!bHasGap}
  else {
    TAxis* axisEta = inputSparse->GetAxis(3);
    if(!axisEta) { Error("TAxis 'axisEta' does not exists!","MakeSparseSlices"); return kFALSE; }

    // positive POIs
    axisEta->SetRangeUser(dEtaGap/2.0,axisEta->GetXmax());
    TH3D* histEntriesPos = (TH3D*) inputSparse->Projection(1,2,0);
    if(!histEntriesPos) { Error("Projection 'histEntriesPos' failed!","MakeSparseSlices"); return kFALSE; }
    trashCol.Add(histEntriesPos);

    // negative POIs
    axisEta->SetRangeUser(axisEta->GetXmin(),-dEtaGap/2.0);
    TH3D* histEntriesNeg = (TH3D*) inputSparse->Projection(1,2,0);
    if(!histEntriesNeg) { Error("Projection 'histEntriesNeg' failed!","MakeSparseSlices"); return kFALSE; }
    trashCol.Add(histEntriesNeg);

    if(task->fMergePosNeg) {
      TList* listMerge = new TList();
      // No need to add ownership since it is collected by trashCol TList
      listMerge->Add(histEntriesPos);
      listMerge->Add(histEntriesNeg);
      histEntries = (TH3D*) MergeListProfiles(listMerge);
      delete listMerge;
      trashCol.Add(histEntries);
    } // end-if {fMergePosNeg}
    else {
      // loading single histo (positive by default)
      if(task->fInputTag.EqualTo("")) {
        histEntries = histEntriesPos;
      }
      else if (task->fInputTag.EqualTo("Neg")) {
        histEntries = histEntriesNeg;
      }
      else {
        Error(Form("Invalid InputTag '%s'!",task->fInputTag.Data()),"ProcessReconstructed");
        return kFALSE;
      }
    } // end-else {fMergePosNeg}
  } // end-else {!bHasGap}

  if(!histEntries) { Error("Histo 'histEntries' failed!","MakeSparseSlices"); return kFALSE; }
  // histEntries (TH3D*) ready for slicing in mult & pt bins (NB should be put in trashCol)

  Int_t iNumBinsMult = fiNumMultBins;
  Int_t iNumBinsPt = task->fNumPtBins;

  TAxis* axisMult = histEntries->GetXaxis();
  TAxis* axisPt = histEntries->GetYaxis();

  for(Int_t iBinMult(0); iBinMult < iNumBinsMult; ++iBinMult) {
    const Int_t iBinMultLow = axisMult->FindFixBin(fdMultBins[iBinMult]);
    const Int_t iBinMultHigh = axisMult->FindFixBin(fdMultBins[iBinMult+1]) - 1;

    for(Int_t iBinPt(0); iBinPt < iNumBinsPt; ++iBinPt) {
      const Double_t dEdgePtLow = task->fPtBinsEdges[iBinPt];
      const Double_t dEdgePtHigh = task->fPtBinsEdges[iBinPt+1];
      const Int_t iBinPtLow = axisPt->FindFixBin(dEdgePtLow);
      const Int_t iBinPtHigh = axisPt->FindFixBin(dEdgePtHigh) - 1;

      TH1D* histInvMass = (TH1D*) histEntries->ProjectionZ(Form("%s_mult%d_pt%d",outName,iBinMult,iBinPt),iBinMultLow,iBinMultHigh,iBinPtLow,iBinPtHigh,"e");
      if(!histInvMass) { Error("Projection 'histInvMass' failed!","MakeSparseSlices"); return kFALSE; }

      if(task->fRebinInvMass > 1) {
        histInvMass->Rebin(task->fRebinInvMass);
      }

      outList->Add(histInvMass);
    } // end-for {binPt}
  } // end-for {binMult}

  Debug("Successfull!","MakeSparseSlices");
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::SaveSlices(FlowTask* task)
{
  if(!task) { Error("FlowTask does not exist!","SaveSlices"); return kFALSE; }

  TList* listProfiles = task->fListProfiles;
  TList* listHistos = task->fListHistos;

  gSystem->mkdir(Form("%s/slices/",fsOutputFilePath.Data()),1);

  ffSlicesFile->cd();
  listProfiles->Write(Form("%s_listProfiles",task->GetSpeciesName().Data()),TObject::kSingleKey);
  listHistos->Write(Form("%s_listHistos",task->GetSpeciesName().Data()),TObject::kSingleKey);

  TLatex latex;
  latex.SetNDC();

  TCanvas* canSlices = new TCanvas("canSlices", "canSlices", 1600,800);

  for(Int_t iMultBin(0); iMultBin < fiNumMultBins; ++iMultBin) {
    for(Int_t iPtBin(0); iPtBin < task->fNumPtBins; ++iPtBin) {

      TString sHistoName = Form("hInvMass_mult%d_pt%d",iMultBin,iPtBin);
      TString sProfName = Form("%s_Pos_sample0_mult%d_pt%d",task->fMixedDiff.Data(),iMultBin,iPtBin);

      TH1D* hInvMass = (TH1D*) listHistos->FindObject(sHistoName.Data());
      if(!hInvMass) { Error(Form("Inv.Mass slice '%s' not found!",sHistoName.Data()),"SaveSlices"); listHistos->ls(); delete canSlices; return kFALSE; }

      TProfile* hFlowMass = (TProfile*) listProfiles->FindObject(sProfName.Data());
      if(!hFlowMass) { Error(Form("Flow.Mass slice '%s' not found!",sProfName.Data()),"SaveSlices"); listProfiles->ls(); delete canSlices; return kFALSE; }

      canSlices->Clear();
      canSlices->Divide(2,1);

      canSlices->cd(1);
      hInvMass->GetXaxis()->SetTitle("M_{inv} (GeV/c^{2})");
      hInvMass->SetMarkerStyle(20);
      hInvMass->SetStats(0);
      hInvMass->SetMinimum(0);
      hInvMass->DrawCopy();
      latex.DrawLatex(0.17,0.80,Form("#color[9]{Enries = %.3g}",hInvMass->GetEntries()));
      latex.DrawLatex(0.17,0.75,Form("#color[9]{mult %.0f-%.0f | pt %0.1f - %0.1f}",fdMultBins[iMultBin],fdMultBins[iMultBin+1],task->fPtBinsEdges[iPtBin],task->fPtBinsEdges[iPtBin+1]));

      canSlices->cd(2);
      hFlowMass->GetXaxis()->SetTitle("M_{inv} (GeV/c^{2})");
      hFlowMass->SetMarkerStyle(20);
      hFlowMass->SetStats(0);
      hFlowMass->DrawCopy();

      canSlices->SaveAs(Form("%s/slices/slices_%s_%s.%s",fsOutputFilePath.Data(),task->GetSpeciesName().Data(),sProfName.Data(),fsOutputFileFormat.Data()),fsOutputFileFormat.Data());
    }
  }

  delete canSlices;
  return kTRUE;
}
//_____________________________________________________________________________
TH1D* ProcessUniFlow::TestRebin(TH1D* hOrig, FlowTask* task)
{
  if(!hOrig) { Error("Original histogram not found!","TestRebin"); return 0x0; }
  if(!task) { Error("Task not found!","TestRebin"); return 0x0; }
  TH1D* hRebin = (TH1D*) hOrig->Rebin(fiNumMultBins,Form("%s_testRebin",hOrig->GetName()),fdMultBins.data());

  // rebinning
  Short_t numBins = fiNumMultBins;
  Double_t* multBins = fdMultBins.data();
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
  if(!prof3dorig) { Error("Input profile does not exists!","Project3DProfile"); return 0x0; }

  TProfile3D* prof3d = (TProfile3D*) prof3dorig->Clone();
  if(!prof3d) { Error("Cloning failed!","Project3DProfile"); return 0x0; }

  Int_t iBinFirst = prof3d->GetXaxis()->GetFirst();
  Int_t iBinLast = prof3d->GetXaxis()->GetLast();
  Int_t iNumBins = prof3d->GetNbinsX();
  Int_t iNumBinsAxis = prof3d->GetXaxis()->GetNbins();
  // printf("Bins:  %d - %d (%d | %d) \n", iBinFirst,iBinLast,iNumBins,iNumBinsAxis);

  // // making 3d hist from 3d profile
  // TH3D* hist3d = prof3d->ProjectionXYZ();   //NOTE do not care about range !!!
  // TH3D* hist3d_entry = prof3d->ProjectionXYZ("hist3d_entry","B");   //NOTE do not care about range !!!
  // TH3D* hist3d_weight = prof3d->ProjectionXYZ("hist3d_weight","W");   //NOTE do not care about range !!!

  TProfile2D* prof2d_test = DoProjectProfile2D(prof3d,Form("%s_px",prof3d->GetName()),prof3d->GetTitle(),prof3d->GetYaxis(),prof3d->GetZaxis(),1,0,0);
  if(!prof2d_test) { Error("DoProjectProfile2D failed!","Project3DProfile"); delete prof3d; return 0x0; }
  prof2d_test->GetXaxis()->SetTitle(prof3d->GetZaxis()->GetTitle());
  prof2d_test->GetYaxis()->SetTitle(prof3d->GetYaxis()->GetTitle());

  // prof2d_test->Draw("colz");
  delete prof3d;
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
 bool useWeights = (h3->GetBinSumw2()->fN != 0); //array elements
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

     Double_t stats[TH1::kNstat] = {0};
     h3->GetStats(stats);
     double dfTsumw = stats[0];

     if (h3->IsA() == TH3F::Class() ) eps = 1.E-6;
     if (dfTsumw != 0 && TMath::Abs( dfTsumw - totcont) <  TMath::Abs(dfTsumw) * eps) resetStats = false;

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
        // since this  is the only way to CalcCum them
        Double_t entries =  h2->GetEffectiveEntries();
        if (!computeErrors) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
        h2->SetEntries( entries );
     }
     else {
        h2->SetEntries( h3->GetEntries() );
     }


     return h2;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::FitInvMass(TH1* hist, FlowTask* task, TF1& fitOutSig, TF1& fitOutBg)
{
  if(!hist) { Error("Input histo not found!","FitInvMass"); return kFALSE; }
  if(!task) { Error("FlowTask not found!","FitInvMass"); return kFALSE; }

  FlowTask::PartSpecies species = task->fSpecies;
  if(species != FlowTask::kPhi && species != FlowTask::kK0s && species != FlowTask::kLambda) { Error("Invalid species!","FitInvMass"); return kFALSE; }

  Double_t dMassRangeLow = hist->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hist->GetXaxis()->GetXmax();
  Double_t dMaximum = hist->GetMaximum();

  Int_t iNpx = 10000;
  TString sFitOptMass = "RNL";

  TString sMassBG = TString(); Int_t iNumParsMassBG = 0; // function for inv. mass dist. (BG component)
  TString sMassSig = TString();  Int_t iNumParsMassSig = 0; // function for inv. mass dist. (sig component)

  Int_t iParMass = 0;
  Int_t iParWidth = 0;
  Int_t iParWidth_2 = 0;

  std::vector<Double_t> dParDef;
  std::vector<Double_t> dParLimLow;
  std::vector<Double_t> dParLimHigh;

  if(species == FlowTask::kPhi) {
    Debug("Setting paramters for Phi","FitInvMass");

    dMassRangeLow = 0.994;
    // dMassRangeHigh = 1.134;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::BreitWigner(x,[5],[6])"; iNumParsMassSig = 3;

    iParMass = 5;
    iParWidth = 6;

    dParDef =     {dMaximum,1000,1.0,0.0,   dMaximum,1.019445,0.0046};
    dParLimLow =  {-1,-1,-1,-1,    0.0,1.018,0.001};
    dParLimHigh = {-1,-1,-1,-1,  1.2*dMaximum,1.022,0.008};

    // assignment to external arrays
  }

  if(species == FlowTask::kK0s)
  {
    Debug("Setting paramters for K0s","FitInvMass");

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;

    dParDef =       {1.0,1.0,1.0,1.0,   dMaximum,0.4976,0.003,dMaximum,0.01};
    dParLimLow =    {-1,-1,-1,-1,    0.0,0.48,0.003,0.0,0.003};
    dParLimHigh =   {-1,-1,-1,-1,  1.2*dMaximum,0.52,0.006,2.0*dMaximum,0.015};
  }

  if(species == FlowTask::kLambda)
  {
    Debug("Setting paramters for Lambda","FitInvMass");

    dMassRangeLow = 1.096;
    dMassRangeHigh = 1.150;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;

    dParDef = {1.0,1.0,1.0,1.0,   dMaximum,1.115, 0.001,dMaximum,0.01};
    dParLimLow = {-1,-1,-1,-1,    0.0,1.10,0.001,0.0,0.001};
    dParLimHigh = {-1,-1,-1,-1,  1.2*dMaximum,1.13,0.008,2.0*dMaximum,0.01};
  }

  // check if parametrisation is setup manually
  if(task->fFlowFitRangeLow > 0.0) { dMassRangeLow = task->fFlowFitRangeLow; }
  if(task->fFlowFitRangeHigh > 0.0) { dMassRangeHigh = task->fFlowFitRangeHigh; }
  if(task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow >= task->fFlowFitRangeHigh) { Error("Wrong fitting ranges set!","FitInvMass"); return kFALSE; }

  Bool_t bUserPars = kFALSE;
  if(task->fNumParMassSig > 0) { bUserPars = kTRUE; sMassSig = task->fFlowFitMassSig; iNumParsMassSig = task->fNumParMassSig; Debug(" Task massSig set","FitInvMass"); }
  if(task->fNumParMassBG > 0) { bUserPars = kTRUE; sMassBG = task->fFlowFitMassBG; iNumParsMassBG = task->fNumParMassBG; Debug(" Task massBG set","FitInvMass"); }
  if(bUserPars && (task->fNumParMassSig == 0 || task->fNumParMassBG == 0)) { Error("Only a subset of functions has been changed. Provide all, or non.","FitInvMass"); return kFALSE; }


  if(bUserPars) {
    Info("Setting UserParameters","FitInvMass");
    dParDef.clear();
    dParLimLow.clear();
    dParLimHigh.clear();

    Int_t iNumParTot = task->fNumParMassSig + task->fNumParMassBG;
    for(Int_t par(0); par < iNumParTot; ++par)
    {
      dParDef.push_back(task->fFitParDefaults[par]);
      dParLimLow.push_back(task->fFitParLimLow[par]);
      dParLimHigh.push_back(task->fFitParLimHigh[par]);
    }
  }

  Int_t iNumParTot = iNumParsMassSig+iNumParsMassBG;
  Int_t iNumParDefs = dParDef.size();
  Int_t iNumParLimLow = dParLimLow.size();
  Int_t iNumParLimHigh = dParLimHigh.size();

  if(iNumParDefs != iNumParTot) { Error(Form("Length of dParDef array does not match number of parameters (%d != %d)",iNumParDefs,iNumParTot),"FitInvMass"); return kFALSE; }
  if(iNumParDefs != iNumParLimLow) { Error(Form("Different length of arrays with parameter defauls and low limit values (%d != %d).",iNumParDefs,iNumParLimLow),"FitInvMass"); return kFALSE; }
  if(iNumParDefs != iNumParLimHigh) { Error(Form("Different length of arrays with parameter defauls and high limit values (%d != %d).",iNumParDefs,iNumParLimHigh),"FitInvMass"); return kFALSE; }

  // check the output of the vector assigment
  if(fbDebug) {
    Debug("Fittin setting done","FitInvMass");
    printf("Form: (%s) + (%s)\n",sMassBG.Data(), sMassSig.Data());
    for(Int_t par(0); par < iNumParDefs; ++par) { printf("  par %d: %g (%g<%g)\n",par, dParDef.at(par), dParLimLow.at(par), dParLimHigh.at(par)); }
  }

  // === Initialision ===
  Int_t iNumParMass = iNumParTot;

  // master formula used in the fitting procedure
  if(!fbDebug) { sFitOptMass += "Q"; } // quite fitting option if NOT in debug

  TString sFuncMass = Form("%s + %s",sMassBG.Data(),sMassSig.Data());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh), "FitInvMass");
  Debug(Form("Fit :\n    %s",sFuncMass.Data()), "FitInvMass");

  // changes the axis
  hist->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // === Fitting procedure ===

  //  fitting invariant mass distribution
  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetNpx(iNpx);

  for(Int_t par(0); par < iNumParMass; ++par) {
    fitMass->SetParameter(par, dParDef.at(par));

    Double_t dLimLow = dParLimLow.at(par);
    Double_t dLimHigh = dParLimHigh.at(par);

    if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
    else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"FitInvMass"); return kFALSE; }
  }

  // fitting
  Int_t nfitsA = 1;
  Bool_t bFitOK = kFALSE;

  TVirtualFitter::SetMaxIterations(10000);

  while(!bFitOK && (nfitsA < 30))
  {
    if(nfitsA > 1)
    {
      fitMass->SetParameter(0, fitMass->GetParameter(0)/nfitsA);

      for(Int_t par(1); par < iNumParMass; ++par)
      {
        fitMass->SetParameter(par, fitMass->GetParameter(par));

        Double_t dLimLow = dParLimLow.at(par);
        Double_t dLimHigh = dParLimHigh.at(par);

        if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
        else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"FitInvMass"); return kFALSE; }
      }
    }

    hist->Fit(fitMass, sFitOptMass.Data());

    TString statusA = gMinuit->fCstatu.Data();

    if(statusA.Contains("CONVERGED")) { bFitOK = kTRUE; }
    nfitsA++;
  }

  // === Extracting fitting components to separated TF1's ===

  TF1 fitBg = TF1("fitMassBG",sMassBG.Data(),dMassRangeLow,dMassRangeHigh);
  fitBg.SetLineColor(kBlue);
  fitBg.SetLineStyle(2);

  TF1 fitSig = TF1("fitMassSig", sMassSig.Data(), dMassRangeLow,dMassRangeHigh);
  fitSig.SetLineColor(kGreen+2);
  fitSig.SetLineStyle(2);

  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar)
  {
    fitBg.SetParameter(iPar, fitMass->GetParameter(iPar));
    fitBg.SetParError(iPar, fitMass->GetParError(iPar));

    fitSig.SetParameter(iPar, 0.0);
    fitSig.SetParError(iPar, 0.0);
  }

  for(Int_t iPar(iNumParsMassBG); iPar < iNumParTot; ++iPar)
  {
    fitSig.SetParameter(iPar, fitMass->GetParameter(iPar));
    fitSig.SetParError(iPar, fitMass->GetParError(iPar));
  }

  fitOutSig = fitSig;
  fitOutBg = fitBg;

  if(!bFitOK) { Error(Form("Inv.mass fit does not converged (%d iterations)",nfitsA)); delete fitMass; return kFALSE; }
  Info(Form("Inv.mass distribution fit: SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g; %d iterations)",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF(),fitMass->GetProb(),nfitsA), "FitInvMass");

  delete fitMass;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::FitCorrelations(TH1* hist, FlowTask* task, TF1& fitOutSig, TF1& fitOutBg, TF1& fitInSig, TF1& fitInBg)
{
  if(!hist) { Error("Input histo not found!","FitCorrelations"); return kFALSE; }
  if(!task) { Error("FlowTask not found!","FitCorrelations"); return kFALSE; }

  FlowTask::PartSpecies species = task->fSpecies;
  if(species != FlowTask::kPhi && species != FlowTask::kK0s && species != FlowTask::kLambda) { Error("Invalid species!","FitCorrelations"); return kFALSE; }

  // === Fitting parametrisation (species dependent default) ===

  Double_t dMassRangeLow = hist->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hist->GetXaxis()->GetXmax();
  Double_t dMaximum = hist->GetMaximum();

  Int_t iNpx = 10000;
  TString sFitOptFlow = "RN";

  TString sMassBG = TString(); Int_t iNumParsMassBG = 0; // function for inv. mass dist. (BG component)
  TString sMassSig = TString();  Int_t iNumParsMassSig = 0; // function for inv. mass dist. (sig component)
  TString sFlowBG = TString();  Int_t iNumParsFlowBG = 0; // function for flow-mass (BG component)

  Int_t iParMass = 0;
  Int_t iParWidth = 0;
  Int_t iParWidth_2 = 0;

  if(species == FlowTask::kPhi)
  {
    Debug("Setting paramters for Phi","FitCorrelations");

    dMassRangeLow = 0.994;
    // dMassRangeHigh = 1.134;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::BreitWigner(x,[5],[6])"; iNumParsMassSig = 3;

    iParMass = 5;
    iParWidth = 6;
  }

  if(species == FlowTask::kK0s)
  {
    Debug("Setting paramters for K0s","FitCorrelations");

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;
  }

  if(species == FlowTask::kLambda)
  {
    Debug("Setting paramters for Lambda","FitCorrelations");

    dMassRangeLow = 1.096;
    dMassRangeHigh = 1.150;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;
  }

  // NB: Because fitInSig includes by construction BG terms
  if(fitInSig.GetNpar() != (iNumParsMassSig + iNumParsMassBG)) { Error(Form("Wrong number of mass signal parameters: %d (%d expected)", iNumParsMassSig+iNumParsMassBG,fitInSig.GetNpar()),"FitCorrelations"); return kFALSE; }

  std::vector<Double_t> dParDef;
  std::vector<Double_t> dParLimLow;
  std::vector<Double_t> dParLimHigh;

  // Species (independent) flow shape
  sFlowBG = Form("[%d]*x+[%d]", iNumParsMassSig+iNumParsMassBG,iNumParsMassSig+iNumParsMassBG+1); iNumParsFlowBG = 2;

  dParDef = {0.0,1.0};
  dParLimLow = {-1,-1};
  dParLimHigh = {-1,-1};

  // check if parametrisation is setup manually
  // TODO: can be extracted from fitIn ???
  if(task->fFlowFitRangeLow > 0.0) { dMassRangeLow = task->fFlowFitRangeLow; }
  if(task->fFlowFitRangeHigh > 0.0) { dMassRangeHigh = task->fFlowFitRangeHigh; }
  if(task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow >= task->fFlowFitRangeHigh) { Error("Wrong fitting ranges set!","FitCorrelations"); return kFALSE; }

  Bool_t bUserPars = kFALSE;
  if(task->fNumParFlowBG > 0) { bUserPars = kTRUE; sFlowBG = task->fFlowFitFlowBG; iNumParsFlowBG = task->fNumParFlowBG; Debug(" Task flowBG set","FitCorrelations"); }

  if(bUserPars)
  {
    Warning("Setting User parameters for correlations not implemented yet","FitCorrelations");
    // TODO: Due to way how default parameters are passed -> need reimplementation
    // using std::vector<> a = {} syntax: Setparames({2,1,2}) ...

    // dParDef.clear();
    // dParLimLow.clear();
    // dParLimHigh.clear();
    //
    // Int_t iNumParTot = task->fNumParFlowBG;
    // for(Int_t par(0); par < iNumParTot; ++par)
    // {
    //   dParDef.push_back(task->fFitParDefaults[par]);
    //   dParLimLow.push_back(task->fFitParLimLow[par]);
    //   dParLimHigh.push_back(task->fFitParLimHigh[par]);
    // }
  }

  Int_t iNumParDefs = dParDef.size();
  Int_t iNumParLimLow = dParLimLow.size();
  Int_t iNumParLimHigh = dParLimHigh.size();

  // check the output of the vector assigment
  if(fbDebug) {
    Debug("Fittin setting done","FitCorrelations");
    printf("Form: %s\n",sFlowBG.Data());
    for(Int_t par(0); par < iNumParDefs; ++par) { printf("  par %d: %g (%g<%g)\n",par, dParDef.at(par), dParLimLow.at(par), dParLimHigh.at(par)); }
  }

  // === Initialision ===
  Int_t iNumParMass = iNumParsMassSig+iNumParsMassBG;
  Int_t iParFlow = iNumParsMassBG+iNumParsMassSig+iNumParsFlowBG; // index of Flow (vn/dn) parameter

  // master formula used in the fitting procedure
  if(!fbDebug) { sFitOptFlow += "Q"; } // quite fitting option if NOT in debug

  TString sFuncMass = Form("(%s) + (%s)",sMassBG.Data(),sMassSig.Data());
  TString sFuncVn = Form("[%d]*(%s)/(%s + %s) + (%s)*(%s)/(%s + %s)", iParFlow, sMassSig.Data(), sMassSig.Data(), sMassBG.Data(), sFlowBG.Data(),sMassBG.Data(),sMassSig.Data(),sMassBG.Data());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh), "FitCorrelations");
  Debug(Form("Fit Dist :\n    %s",sFuncMass.Data()), "FitCorrelations");
  Debug(Form("Fit Flow :\n    %s\n",sFuncVn.Data()), "FitCorrelations");

  // changes the axis
  // hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  // hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // === Fitting procedure ===

  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);

  // fixing fraction from input fits

  for(Int_t par(0); par < iNumParsMassBG; ++par) { fitVn->FixParameter(par, fitInBg.GetParameter(par)); }
  for(Int_t par(iNumParsMassBG); par < iNumParMass; ++par) { fitVn->FixParameter(par, fitInSig.GetParameter(par)); }
  for(Int_t par(iNumParMass); par < iParFlow; ++par)
  {
    // Here par-iNumParMass is to account for a fact that dParDef takes only flow part (vector index != parameter index)
    fitVn->SetParameter(par, dParDef.at(par-iNumParMass));
    Debug(Form("Parameter %d : %f",par,dParDef.at(par-iNumParMass) ),"FitCorrelations");
    Double_t dLimLow = dParLimLow.at(par-iNumParMass);
    Double_t dLimHigh = dParLimHigh.at(par-iNumParMass);

    if(dLimLow > -1.0 && dLimHigh > -1.0) { fitVn->SetParLimits(par, dLimLow,dLimHigh); }
    else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Flow-mass: Only one of the parameter limits is set (par %d). Fix this!",par),"FitCorrelations"); return kFALSE; }
  }

  fitVn->SetParameter(iParFlow, 0.5);
  fitVn->SetParLimits(iParFlow, 0.0,1.0);

  // NB: Currently only one iteration
  // // fitting
  // Int_t nfitsA = 1;
  // Bool_t bFitOK = kFALSE;
  //
  // while(!bFitOK && (nfitsA < 15))
  // {
  //   if(nfitsA > 1)
  //   {
  //     fitMass->SetParameter(0, fitMass->GetParameter()/nfitsA);
  //
  //     for(Int_t par(1); par < iNumParMass; ++par)
  //     {
  //       fitMass->SetParameter(par, fitMass->GetParameter(par));
  //
  //       Double_t dLimLow = dParLimLow.at(par);
  //       Double_t dLimHigh = dParLimHigh.at(par);
  //
  //       if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
  //       else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"FitInvMass"); return kFALSE; }
  //     }
  //   }

  //   hist->Fit(fitVn, sFitOptFlow.Data());
  //
  //   TString statusA = gMinuit->fCstatu.Data();
  //   if(statusA.Contains("CONVERGED")) { bFitOK = kTRUE; }
  //   nfitsA++;
  // }

  hist->Fit(fitVn, sFitOptFlow.Data());

  if(!gMinuit->fCstatu.Contains("CONVERGED") ) { Error(Form("Flow-mass fit does not converged within iterations limit (1)!"), "FitCorrelations"); delete fitVn; return kFALSE; }
  Info(Form("Flow-mass fit: SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g)",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF(),fitVn->GetProb()), "FitCorrelations");

  // saving flow to output
  Double_t dFlow = fitVn->GetParameter(iParFlow);
  Double_t dFlowError = fitVn->GetParError(iParFlow);
  Double_t dFlowRel = 0.0; if(TMath::Abs(dFlow) > 0.0) { dFlowRel = dFlowError / dFlow; }
  Debug(Form("Final flow: %g +- %g (rel. %.3f)", dFlow,dFlowError,dFlowRel), "FitCorrelations");

  // === Extracting fitting components to separated TF1's ===

  TF1 fitFlowBg = TF1("fitFlowBg", Form("(%s)*(%s)/(%s + %s)", sFlowBG.Data(), sMassBG.Data(), sMassSig.Data(), sMassBG.Data()), dMassRangeLow,dMassRangeHigh);
  fitFlowBg.SetLineColor(kBlue);
  fitFlowBg.SetLineStyle(2);
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar)
  {
    fitFlowBg.SetParameter(iPar, fitInBg.GetParameter(iPar));
    fitFlowBg.SetParError(iPar, fitInBg.GetParError(iPar));
  }
  for(Int_t iPar(iNumParsMassBG); iPar < iNumParMass; ++iPar)
  {
    fitFlowBg.SetParameter(iPar, fitInSig.GetParameter(iPar));
    fitFlowBg.SetParError(iPar, fitInSig.GetParError(iPar));
  }
  for(Int_t iPar(iNumParMass); iPar < iParFlow; ++iPar)
  {
    fitFlowBg.SetParameter(iPar, fitVn->GetParameter(iPar));
    fitFlowBg.SetParError(iPar, fitVn->GetParError(iPar));
  }

  TF1 fitFlowSig = TF1("fitFlowSig", Form("[%d]*(%s)/(%s + %s)", iParFlow, sMassSig.Data(), sMassSig.Data(), sMassBG.Data()), dMassRangeLow,dMassRangeHigh);
  fitFlowSig.SetLineColor(kGreen+2);
  fitFlowSig.SetLineStyle(2);
  for(Int_t iPar(0); iPar < iNumParMass; ++iPar)
  {
    fitFlowSig.SetParameter(iPar, fitFlowBg.GetParameter(iPar));
    fitFlowSig.SetParError(iPar, fitFlowBg.GetParError(iPar));
  }
  for(Int_t iPar(iNumParMass); iPar < iParFlow; ++iPar)
  {
    fitFlowSig.SetParameter(iPar, 0.0);
    fitFlowSig.SetParError(iPar, 0.0);
  }
  fitFlowSig.SetParameter(iParFlow, fitVn->GetParameter(iParFlow));
  fitFlowSig.SetParError(iParFlow, fitVn->GetParError(iParFlow));

  fitOutSig = fitFlowSig;
  fitOutBg = fitFlowBg;

  delete fitVn;

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowOneGo"); return kFALSE; }
  if(listFits->GetEntries() > 0) { Error("TList not empty!","ExtractFlowOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowOneGo"); return kFALSE; }
  if(task->fSpecies == FlowTask::kPhi && task->fFlowFitPhiSubtLS && !hInvMassBG) { Error("Inv. Mass (BG) histogram does not exists!","ExtractFlowOneGo"); return kFALSE; }
  if(task->fFlowFitPhiSubtLS) { Error("Phi like-sign subtraction not implemented ATM. Please turn the switch off.","ExtractFlowOneGo"); return kFALSE; }

  // === Fitting parametrisation (species dependent default) ===

  Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  Double_t dMaximum = hInvMass->GetMaximum();

  Int_t iNpx = 10000;
  TString sFitOptMass = "RNL";
  TString sFitOptFlow = "RN";

  TString sMassBG = TString(); Int_t iNumParsMassBG = 0; // function for inv. mass dist. (BG component)
  TString sMassSig = TString();  Int_t iNumParsMassSig = 0; // function for inv. mass dist. (sig component)
  TString sFlowBG = TString();  Int_t iNumParsFlowBG = 0; // function for flow-mass (BG component)

  Int_t iParMass = 0;
  Int_t iParWidth = 0;
  Int_t iParWidth_2 = 0;

  std::vector<Double_t> dParDef;
  std::vector<Double_t> dParLimLow;
  std::vector<Double_t> dParLimHigh;

  if(task->fSpecies == FlowTask::kPhi)
  {
    dMassRangeLow = 0.994;
    // dMassRangeHigh = 1.134;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::BreitWigner(x,[5],[6])"; iNumParsMassSig = 3;
    sFlowBG = "[7]*x+[8]"; iNumParsFlowBG = 2;

    iParMass = 5;
    iParWidth = 6;

    Double_t dDef[] =      {1.0,1.0,1.0,1.0,   dMaximum,1.019445,0.0046, 1.0,1.0};
    Double_t dLimLow[] =   {-1,-1,-1,-1,    0.0,1.018,0.001, -1,-1};
    Double_t dLimHigh[] =  {-1,-1,-1,-1,  2.0*dMaximum,1.022,0.006,  -1,-1};

    // assignment to external arrays
    for(Int_t par(0); par < (Int_t) (sizeof(dDef)/sizeof(dDef[0])); ++par) { dParDef.push_back(dDef[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimLow)/sizeof(dLimLow[0])); ++par) { dParLimLow.push_back(dLimLow[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimHigh)/sizeof(dLimHigh[0])); ++par) { dParLimHigh.push_back(dLimHigh[par]); }
  }
  else if(task->fSpecies == FlowTask::kK0s)
  {
    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;
    sFlowBG = "[9]*x+[10]"; iNumParsFlowBG = 2;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;

    Double_t dDef[] =      {1.0,1.0,1.0,1.0,   dMaximum,0.4976,0.003,dMaximum,0.01, 1.0,1.0};
    Double_t dLimLow[] =   {-1,-1,-1,-1,    0.0,0.48,0.003,0.0,0.003, -1,-1};
    Double_t dLimHigh[] =  {-1,-1,-1,-1,  2.0*dMaximum,0.52,0.006,2.0*dMaximum,0.01,  -1,-1};

    // assignment to external arrays
    for(Int_t par(0); par < (Int_t) (sizeof(dDef)/sizeof(dDef[0])); ++par) { dParDef.push_back(dDef[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimLow)/sizeof(dLimLow[0])); ++par) { dParLimLow.push_back(dLimLow[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimHigh)/sizeof(dLimHigh[0])); ++par) { dParLimHigh.push_back(dLimHigh[par]); }
  }
  else if(task->fSpecies == FlowTask::kLambda)
  {
    dMassRangeLow = 1.096;
    dMassRangeHigh = 1.150;

    sMassBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; iNumParsMassBG = 4;
    sMassSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; iNumParsMassSig = 5;
    sFlowBG = "[9]*x+[10]"; iNumParsFlowBG = 2;

    iParMass = 5;
    iParWidth = 6;
    iParWidth_2 = 8;

    Double_t dDef[] =      {1.0,1.0,1.0,1.0,   dMaximum,1.115, 0.001,dMaximum,0.01, 1.0,1.0};
    Double_t dLimLow[] =   {-1,-1,-1,-1,    0.0,1.10,0.001,0.0,0.001, -1,-1};
    Double_t dLimHigh[] =  {-1,-1,-1,-1,  2.0*dMaximum,1.13,0.008,2.0*dMaximum,0.01,  -1,-1};

    // assignment to external arrays
    for(Int_t par(0); par < (Int_t) (sizeof(dDef)/sizeof(dDef[0])); ++par) { dParDef.push_back(dDef[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimLow)/sizeof(dLimLow[0])); ++par) { dParLimLow.push_back(dLimLow[par]); }
    for(Int_t par(0); par < (Int_t) (sizeof(dLimHigh)/sizeof(dLimHigh[0])); ++par) { dParLimHigh.push_back(dLimHigh[par]); }
  }
  else { Error("Invalid species","ExtractFlowOneGo"); return kFALSE; }

  // check if parametrisation is setup manually
  if(task->fFlowFitRangeLow > 0.0) { dMassRangeLow = task->fFlowFitRangeLow; }
  if(task->fFlowFitRangeHigh > 0.0) { dMassRangeHigh = task->fFlowFitRangeHigh; }
  if(task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow > 0.0 && task->fFlowFitRangeLow >= task->fFlowFitRangeHigh) { Error("Wrong fitting ranges set!","ExtractFlowOneGo"); return kFALSE; }

  Bool_t bUserPars = kFALSE;
  if(task->fNumParMassSig > 0) { bUserPars = kTRUE; sMassSig = task->fFlowFitMassSig; iNumParsMassSig = task->fNumParMassSig; Debug(" Task massSig set","ExtractFlowOneGo"); }
  if(task->fNumParMassBG > 0) { bUserPars = kTRUE; sMassBG = task->fFlowFitMassBG; iNumParsMassBG = task->fNumParMassBG; Debug(" Task massBG set","ExtractFlowOneGo"); }
  if(task->fNumParFlowBG > 0) { bUserPars = kTRUE; sFlowBG = task->fFlowFitFlowBG; iNumParsFlowBG = task->fNumParFlowBG; Debug(" Task flowBG set","ExtractFlowOneGo"); }

  if(bUserPars && (task->fNumParMassSig == 0 || task->fNumParMassBG == 0|| task->fNumParFlowBG == 0)) { Error("Only a subset of functions has been changed. Provide all, or non.","ExtractFlowOneGo"); return kFALSE; }

  if(bUserPars)
  {
    dParDef.clear();
    dParLimLow.clear();
    dParLimHigh.clear();

    Int_t iNumParTot = task->fNumParMassSig + task->fNumParMassBG + task->fNumParFlowBG;
    for(Int_t par(0); par < iNumParTot; ++par)
    {
      dParDef.push_back(task->fFitParDefaults[par]);
      dParLimLow.push_back(task->fFitParLimLow[par]);
      dParLimHigh.push_back(task->fFitParLimHigh[par]);
    }
  }

  Int_t iNumParDefs = dParDef.size();
  Int_t iNumParLimLow = dParLimLow.size();
  Int_t iNumParLimHigh = dParLimHigh.size();

  // check the output of the vector assigment
  Debug("Post ifs","ExtractFlowOneGo");
  for(Int_t par(0); par < iNumParDefs; ++par) { Debug(Form("  par %d: %g (%g<%g)",par, dParDef.at(par), dParLimLow.at(par), dParLimHigh.at(par)),"ExtractFlowOneGo"); }

  if(iNumParDefs != iNumParsMassBG+iNumParsMassSig+iNumParsFlowBG) { Error(Form("Length of dParDef array does not match number of parameters (%d != %d)",iNumParDefs,iNumParsMassBG+iNumParsMassSig+iNumParsFlowBG),"ExtractFlowOneGo"); return kFALSE; }
  if(iNumParDefs != iNumParLimLow) { Error(Form("Different length of arrays with parameter defauls and low limit values (%d != %d).",iNumParDefs,iNumParLimLow),"ExtractFlowOneGo"); return kFALSE; }
  if(iNumParDefs != iNumParLimHigh) { Error(Form("Different length of arrays with parameter defauls and high limit values (%d != %d).",iNumParDefs,iNumParLimHigh),"ExtractFlowOneGo"); return kFALSE; }

  // === Initialision ===
  Int_t iNumParMass = iNumParsMassSig+iNumParsMassBG;
  Int_t iParFlow = iNumParsMassBG+iNumParsMassSig+iNumParsFlowBG; // index of Flow (vn/dn) parameter

  // master formula used in the fitting procedure
  if(!fbDebug) { sFitOptFlow += "Q"; sFitOptMass += "Q"; } // quite fitting option if NOT in debug

  TString sFuncMass = Form("%s + %s",sMassBG.Data(),sMassSig.Data());
  TString sFuncVn = Form("[%d]*(%s)/(%s + %s) + (%s)*(%s)/(%s + %s)", iParFlow, sMassSig.Data(), sMassSig.Data(), sMassBG.Data(), sFlowBG.Data(),sMassBG.Data(),sMassSig.Data(),sMassBG.Data());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh), "ExtractFlowOneGo");
  Debug(Form("Fit Dist :\n    %s",sFuncMass.Data()), "ExtractFlowOneGo");
  Debug(Form("Fit Flow :\n    %s\n",sFuncVn.Data()), "ExtractFlowOneGo");

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // === Fitting procedure ===

  // fitting invariant mass distribution
  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetNpx(iNpx);

  for(Int_t par(0); par < iNumParMass; ++par)
  {
    fitMass->SetParameter(par, dParDef.at(par));

    Double_t dLimLow = dParLimLow.at(par);
    Double_t dLimHigh = dParLimHigh.at(par);

    if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
    else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"ExtractFlowOneGo"); return kFALSE; }
  }

  hInvMass->Fit(fitMass, sFitOptMass.Data());

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 15))
  {
    fitMass->SetParameter(0, fitMass->GetParameter(0)/nfitsA);
    for(Int_t par(0); par < iNumParMass; ++par)
    {
      fitMass->SetParameter(par, fitMass->GetParameter(par));

      Double_t dLimLow = dParLimLow.at(par);
      Double_t dLimHigh = dParLimHigh.at(par);

      if(dLimLow > -1.0 && dLimHigh > -1.0) { fitMass->SetParLimits(par, dLimLow, dLimHigh); }
      else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Inv.mass (def): Only one of the parameter limits is set (par %d : %g :%g < %g). Fix this!",par,dParDef[par], dLimLow, dLimHigh),"ExtractFlowOneGo"); return kFALSE; }
    }

    hInvMass->Fit(fitMass, sFitOptMass.Data());

    statusA = gMinuit->fCstatu.Data();
    nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv.mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("Inv.mass distribution fit: SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g; %d iterations)",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF(),fitMass->GetProb(),nfitsA), "ExtractFlowOneGo");

  // fitting invariant mass distribution
  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  // fixing Nsig & Nbg terms extracted in previous step
  for(Int_t par(0); par < iNumParMass; ++par) { fitVn->FixParameter(par, fitMass->GetParameter(par)); }
  for(Int_t par(iNumParMass); par < iParFlow; ++par)
  {
    fitVn->SetParameter(par, dParDef.at(par));

    Double_t dLimLow = dParLimLow.at(par);
    Double_t dLimHigh = dParLimHigh.at(par);

    if(dLimLow > -1.0 && dLimHigh > -1.0) { fitVn->SetParLimits(par, dLimLow,dLimHigh); }
    else if(dLimLow > -1.0 || dLimHigh > -1.0) { Error(Form("Flow-mass: Only one of the parameter limits is set (par %d). Fix this!",par),"ExtractFlowOneGo"); return kFALSE; }
  }
  hFlowMass->Fit(fitVn, sFitOptFlow.Data());

  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Flow-mass fit does not converged within iterations limit!"), "ExtractFlowOneGo"); return kFALSE; }
  Info(Form("Flow-mass fit: SUCCESSFULL (chi2/ndf = %.3g/%d = %.3g; prob = %0.2g)",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF(),fitVn->GetProb()), "ExtractFlowOneGo");

  // saving flow to output
  dFlow = fitVn->GetParameter(iParFlow);
  dFlowError = fitVn->GetParError(iParFlow);
  Double_t dFlowRel = 0.0; if(TMath::Abs(dFlow) > 0.0) { dFlowRel = dFlowError / dFlow; }
  Debug(Form("Final flow: %g +- %g (rel. %.3f)\n==================================================================", dFlow,dFlowError,dFlowRel), "ExtractFlowOneGo");

  // === Extracting fitting components to separated TF1's ===

  TF1* fitBg = new TF1("fitMassBG",sMassBG.Data(),dMassRangeLow,dMassRangeHigh);
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar)
  {
    fitBg->SetParameter(iPar, fitMass->GetParameter(iPar));
    fitBg->SetParError(iPar, fitMass->GetParError(iPar));
  }

  TF1* fitSig = new TF1("fitMassSig", sMassSig.Data(), dMassRangeLow,dMassRangeHigh);
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);
  for(Int_t iPar(0); iPar < iNumParsMassBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsMassBG); iPar < iNumParsMassBG+iNumParsMassSig; ++iPar)
  {
    fitSig->SetParameter(iPar, fitMass->GetParameter(iPar));
    fitSig->SetParError(iPar, fitMass->GetParError(iPar));
  }

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("(%s)*(%s)/(%s + %s)", sFlowBG.Data(), sMassBG.Data(), sMassSig.Data(), sMassBG.Data()), dMassRangeLow,dMassRangeHigh);
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  for(Int_t iPar(0); iPar < iParFlow; ++iPar)
  {
    fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar));
    fitFlowBg->SetParError(iPar, fitVn->GetParError(iPar));
  }

  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[%d]*(%s)/(%s + %s)", iParFlow, sMassSig.Data(), sMassSig.Data(), sMassBG.Data()), dMassRangeLow,dMassRangeHigh);
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);
  fitFlowSig->SetParameter(iParFlow, fitVn->GetParameter(iParFlow));
  for(Int_t iPar(0); iPar < iNumParsMassBG+iNumParsMassSig; ++iPar)
  {
    fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar));
    fitFlowSig->SetParError(iPar, fitVn->GetParError(iPar));
  }
  for(Int_t iPar(iNumParsMassBG+iNumParsMassSig); iPar < iParFlow; ++iPar) { fitFlowSig->SetParameter(iPar, 0.0); }

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);

  // === Drawing stuff to canvas ===

  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  latex->SetNDC();

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hInvMass->SetMarkerStyle(20);
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitBg->DrawCopy("same");
  fitSig->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(iParMass),fitMass->GetParError(iParMass)));

  if(task->fSpecies == FlowTask::kPhi)
  {
    latex->DrawLatex(0.17,0.70,Form("#color[9]{#Gamma = %.6f #pm %.6f}",fitMass->GetParameter(iParWidth),fitMass->GetParError(iParWidth)));
  }
  else if(task->fSpecies == FlowTask::kK0s || task->fSpecies == FlowTask::kLambda)
  {
    latex->DrawLatex(0.17,0.70,Form("#color[9]{#sigma_{1} = %.6f #pm %.6f}",fitMass->GetParameter(iParWidth),fitMass->GetParError(iParWidth)));
    latex->DrawLatex(0.17,0.65,Form("#color[9]{#sigma_{2} = %.6f #pm %.6f}",fitMass->GetParameter(iParWidth_2),fitMass->GetParError(iParWidth_2)));
  }

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hFlowMass->SetMarkerStyle(20);
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  // fitFlowSig->DrawCopy("same");
  // fitFlowBg->DrawCopy("same");

  TString sResult = "v_{2}"; if(fFlowFitCumulants) { sResult = "d_{2}"; }
  latex->DrawLatex(0.17,0.80,Form("#color[9]{%s = %.4f #pm %.4f}",sResult.Data(),dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{p = %.3g}",fitVn->GetProb()));

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowPhiOneGo(FlowTask* task, TH1* hInvMass, TH1* hInvMassBG, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowPhiOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowPhiOneGo"); return kFALSE; }

  // subtraction of LS?

  // fitting parametrisation
  TString sFuncBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; Int_t iNumParsFuncBG = 4;
  TString sFuncSig = "[4]*TMath::BreitWigner(x,[5],[6])";  Int_t iNumParsFuncSig = 3;

  TString sFuncVn = Form("[9]*(%s)/(%s + %s) + ([7]*x+[8])*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data(), sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data());
  TString sFuncMass = Form("%s + %s",sFuncBG.Data(),sFuncSig.Data());

  // Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  Double_t dMassRangeLow = 0.994;
  // Double_t dMassRangeHigh =1.134;
  Double_t dMaximum = hInvMass->GetMaximum();

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh),"ExtractFlowPhiOneGo");
  Debug(Form("Fit func invmass:\n%s",sFuncVn.Data()));
  Debug(Form("Fit func vn:\n%s",sFuncVn.Data()));


  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetParameters(dMaximum/10.0, 1.0, 1.0, 1.0, dMaximum, 1.019445, 0.0046);
  fitMass->SetNpx(5000);
  fitMass->SetParLimits(5,1.018,1.022);
  fitMass->SetParLimits(6,0.001,0.006);
  hInvMass->Fit(fitMass, "RNL");

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 10))
  {
    fitMass->SetParameters(fitMass->GetParameter(0)/nfitsA, fitMass->GetParameter(1), fitMass->GetParameter(2), fitMass->GetParameter(3), fitMass->GetParameter(4), 1.019445, fitMass->GetParameter(5)*nfitsA);
    fitMass->SetParLimits(5,1.018,1.022);
    fitMass->SetParLimits(6,0.001,0.006);
    hInvMass->Fit(fitMass, "RNL");

    statusA = gMinuit->fCstatu.Data();
    nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("Number of iterations: %d\n",nfitsA));

  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  fitVn->SetParameter(7, 1.0);
  fitVn->SetParameter(8, 1.0);
  fitVn->SetParameter(9, 0.1);

  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitVn->FixParameter(iPar, fitMass->GetParameter(iPar)); }
  hFlowMass->Fit(fitVn, "RN");
  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged!")); return kFALSE; }

  // saving flow to output
  dFlow = fitVn->GetParameter(9);
  dFlowError = fitVn->GetParError(9);
  Info(Form("=================================\n Final flow: %g +- %g\n =================================\n", dFlow,dFlowError));

  TF1* fitBg = new TF1("fitMassBG",sFuncBG.Data(),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitBg->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);

  TF1* fitSig = new TF1("fitMassSig", sFuncSig.Data(), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsFuncBG); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitSig->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("([7]*x+[8])*(%s)/(%s + %s)",sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data()),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowBg->SetParameter(7, fitVn->GetParameter(7));
  fitFlowBg->SetParameter(8, fitVn->GetParameter(8));
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  //
  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[9]*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data()), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowSig->SetParameter(7, 0.0);
  fitFlowSig->SetParameter(8, 0.0);
  fitFlowSig->SetParameter(9, fitVn->GetParameter(9));
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);

  // Drawing stuff
  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hInvMass->SetMarkerStyle(20);
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitBg->DrawCopy("same");
  fitSig->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(5),fitMass->GetParError(5)));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{#Gamma = %.6f #pm %.6f}",fitMass->GetParameter(6),fitMass->GetParError(6)));

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{#phi} (GeV/c^{2})");
  hFlowMass->SetMarkerStyle(20);
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  // fitFlowSig->DrawCopy("same");
  // fitFlowBg->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{v_{2} = %.4f #pm %.4f}",dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowK0sOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowK0sOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowK0sOneGo"); return kFALSE; }


  Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  // Double_t dMassRangeLow = 0.41;
  // Double_t dMassRangeHigh = 0.59;
  Double_t dMaximum = hInvMass->GetMaximum();

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);

  // fitting parametrisation
  TString sFuncBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; Int_t iNumParsFuncBG = 4;
  TString sFuncSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; Int_t iNumParsFuncSig = 5;

  TString sFuncVn = Form("[11]*(%s)/(%s + %s) + ([9]*x+[10])*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data(), sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data());
  TString sFuncMass = Form("%s + %s",sFuncBG.Data(),sFuncSig.Data());

  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh),"ExtractFlowPhiOneGo");
  Debug(Form("Fit func invmass:\n%s",sFuncVn.Data()));
  Debug(Form("Fit func vn:\n%s",sFuncVn.Data()));

  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetNpx(10000);
  fitMass->SetParameters(dMaximum/10.0, 1.0, 1.0, 1.0, dMaximum, 0.4976, 0.001,dMaximum,0.001);
  fitMass->SetParLimits(4, 0, dMaximum*2.0);
  fitMass->SetParLimits(5, 0.48, 0.52);
  fitMass->SetParLimits(6, 0.003, 0.006);
  fitMass->SetParLimits(7, 0, dMaximum*2.0);
  fitMass->SetParLimits(8, 0.003, 0.01);
  hInvMass->Fit(fitMass, "RNL");

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 15)){

      fitMass->SetParameters(fitMass->GetParameter(0)/nfitsA, fitMass->GetParameter(1), fitMass->GetParameter(2), fitMass->GetParameter(3), fitMass->GetParameter(4), 0.4976, fitMass->GetParameter(6)*nfitsA,fitMass->GetParameter(7), fitMass->GetParameter(8)*nfitsA);
      // fitMass->SetParLimits(4, 0, dMaximum*2.0);
      fitMass->SetParLimits(5, 0.48, 0.52);
      fitMass->SetParLimits(6, 0.003, 0.006);
      // fitMass->SetParLimits(7, 0, dMaximum*2.0);
      fitMass->SetParLimits(8, 0.003, 0.01);
      hInvMass->Fit(fitMass, "RNL");

      statusA = gMinuit->fCstatu.Data();
      nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("Number of iterations: %d\n",nfitsA));

  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  fitVn->SetParameter(9, 1.0);
  fitVn->SetParameter(10, 1.0);
  fitVn->SetParameter(11, 0.1);

  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitVn->FixParameter(iPar, fitMass->GetParameter(iPar)); }
  hFlowMass->Fit(fitVn, "RN");
  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged!")); return kFALSE; }

  // saving flow to output
  dFlow = fitVn->GetParameter(11);
  dFlowError = fitVn->GetParError(11);
  Info(Form("=================================\n Final flow: %g +- %g\n =================================\n", dFlow,dFlowError));

  // Drawing stuff
  TF1* fitBg = new TF1("fitBG",sFuncBG.Data(),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitBg->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);

  TF1* fitSig = new TF1("fitSig", sFuncSig.Data(), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsFuncBG); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitSig->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("([7]*x+[8])*(%s)/(%s + %s)",sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data()),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowBg->SetParameter(9, fitVn->GetParameter(9));
  fitFlowBg->SetParameter(10, fitVn->GetParameter(10));
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  //
  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[9]*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data()), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowSig->SetParameter(9, 0.0);
  fitFlowSig->SetParameter(10, 0.0);
  fitFlowSig->SetParameter(11, fitVn->GetParameter(11));
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);


  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitSig->DrawCopy("same");
  fitBg->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(5),fitMass->GetParError(5)));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{#sigma_{1} = %.6f #pm %.6f}",fitMass->GetParameter(6),fitMass->GetParError(6)));
  latex->DrawLatex(0.17,0.65,Form("#color[9]{#sigma_{2} = %.6f #pm %.6f}",fitMass->GetParameter(8),fitMass->GetParError(8)));

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{K^{0}} (GeV/c^{2})");
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{v_{2} = %.4f #pm %.4f}",dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::ExtractFlowLambdaOneGo(FlowTask* task, TH1* hInvMass, TH1* hFlowMass, Double_t &dFlow, Double_t &dFlowError, TCanvas* canFitInvMass, TList* listFits)
{
  if(!task) { Error("Coresponding FlowTask not found!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!hInvMass) { Error("Inv. Mass histogram does not exists!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!hFlowMass) { Error("Flow Mass histogram does not exists!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!canFitInvMass) { Error("Canvas not found!","ExtractFlowLambdaOneGo"); return kFALSE; }
  if(!listFits) { Error("TList for fits not found!","ExtractFlowLambdaOneGo"); return kFALSE; }


  TString sFuncBG = "[0] + [1]*x + [2]*x*x + [3]*x*x*x"; Int_t iNumParsFuncBG = 4;
  TString sFuncSig = "[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[5],[8])"; Int_t iNumParsFuncSig = 5;
  // TString sFuncSig = "gaus(4)+gaus(7)"; Int_t iNumParsFuncSig = 6;

  TString sFuncMass = Form("%s + %s",sFuncBG.Data(),sFuncSig.Data());
  TString sFuncVn = Form("[11]*(%s)/(%s + %s) + ([9]*x+[10])*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data(), sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data());

  // Double_t dMassRangeLow = hInvMass->GetXaxis()->GetXmin();
  // Double_t dMassRangeHigh = hInvMass->GetXaxis()->GetXmax();
  Double_t dMassRangeLow = 1.096;
  Double_t dMassRangeHigh =1.150;
  Double_t dMaximum = hInvMass->GetMaximum();

  // changes the axis
  hInvMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);
  hFlowMass->GetXaxis()->SetRangeUser(dMassRangeLow,dMassRangeHigh);


  Debug(Form("Mass range %g-%g",dMassRangeLow,dMassRangeHigh),"ExtractFlowPhiOneGo");
  Debug(Form("Fit func invmass:\n%s",sFuncVn.Data()));
  Debug(Form("Fit func vn:\n%s",sFuncVn.Data()));

  TF1* fitMass = new TF1(Form("fitMass"), sFuncMass.Data(), dMassRangeLow,dMassRangeHigh);
  fitMass->SetParameters(dMaximum/10.0, 1.0, 1.0, 1.0, dMaximum, 1.115, 0.001,dMaximum, 0.001);
  fitMass->SetNpx(5000);
  fitMass->SetParLimits(4, 0, dMaximum*2.0);
  fitMass->SetParLimits(5, 1.10, 1.13);
  fitMass->SetParLimits(6, 0.001,0.006);
  fitMass->SetParLimits(7, 0, dMaximum*2.0);
  fitMass->SetParLimits(8, 0.001,0.01);
  hInvMass->Fit(fitMass, "RNL");

  // checking the status of convergence
  Int_t nfitsA = 1;
  TString statusA = gMinuit->fCstatu.Data();

  while ((!statusA.Contains("CONVERGED")) && (nfitsA < 15))
  {
    fitMass->SetParameters(fitMass->GetParameter(0)/nfitsA, fitMass->GetParameter(1), fitMass->GetParameter(2), fitMass->GetParameter(3), fitMass->GetParameter(4), 1.115, fitMass->GetParameter(6)*nfitsA, fitMass->GetParameter(7), fitMass->GetParameter(9)*nfitsA);
    fitMass->SetParLimits(4, 0, dMaximum*2.0);
    fitMass->SetParLimits(5, 1.10, 1.13);
    fitMass->SetParLimits(6, 0.001,0.008);
    fitMass->SetParLimits(7, 0, dMaximum*2.0);
    fitMass->SetParLimits(8, 0.001,0.01);
    hInvMass->Fit(fitMass, "RNL");

    statusA = gMinuit->fCstatu.Data();
    nfitsA++;
  }

  if(!statusA.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged (%d iterations)",nfitsA)); return kFALSE; }
  Info(Form("Number of iterations: %d\n",nfitsA));

  TF1* fitVn = new TF1(Form("fitVn"), sFuncVn.Data(), dMassRangeLow,dMassRangeHigh);
  fitVn->SetParameter(9, 1.0);
  fitVn->SetParameter(10, 1.0);
  fitVn->SetParameter(11, 0.1);

  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitVn->FixParameter(iPar, fitMass->GetParameter(iPar)); }
  hFlowMass->Fit(fitVn, "RN");
  if(!gMinuit->fCstatu.Contains("CONVERGED")) { Error(Form("Inv. mass fit does not converged!")); return kFALSE; }

  // saving flow to output
  dFlow = fitVn->GetParameter(11);
  dFlowError = fitVn->GetParError(11);
  Info(Form("=================================\n Final flow: %g +- %g\n =================================\n", dFlow,dFlowError));

  // Drawing stuff
  TF1* fitBg = new TF1("fitBG",sFuncBG.Data(),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitBg->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitBg->SetLineColor(kBlue);
  fitBg->SetLineStyle(2);

  TF1* fitSig = new TF1("fitSig", sFuncSig.Data(), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG; ++iPar) { fitSig->SetParameter(iPar, 0.0); }
  for(Int_t iPar(iNumParsFuncBG); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitSig->SetParameter(iPar, fitMass->GetParameter(iPar)); }
  fitSig->SetLineColor(kGreen+2);
  fitSig->SetLineStyle(2);

  TF1* fitFlowBg = new TF1("fitFlowBG", Form("([9]*x+[10])*(%s)/(%s + %s)",sFuncBG.Data(),sFuncSig.Data(),sFuncBG.Data()),dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowBg->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  fitFlowBg->SetParameter(9, fitVn->GetParameter(9));
  fitFlowBg->SetParameter(10, fitVn->GetParameter(10));
  fitFlowBg->SetLineColor(kBlue);
  fitFlowBg->SetLineStyle(2);
  //
  TF1* fitFlowSig = new TF1("fitFlowSig", Form("[12]*(%s)/(%s + %s)",sFuncSig.Data(), sFuncSig.Data(), sFuncBG.Data()), dMassRangeLow,dMassRangeHigh);
  for(Int_t iPar(0); iPar < iNumParsFuncBG+iNumParsFuncSig; ++iPar) { fitFlowSig->SetParameter(iPar, fitVn->GetParameter(iPar)); }
  // fitFlowSig->SetParameter(10, 0.0);
  // fitFlowSig->SetParameter(11, 0.0);
  fitFlowSig->SetParameter(11, fitVn->GetParameter(11));
  fitFlowSig->SetLineColor(kGreen+2);
  fitFlowSig->SetLineStyle(2);

  // saving fitting related stuff to TList listFits
  listFits->Add(hInvMass);
  listFits->Add(fitMass);
  listFits->Add(fitBg);
  listFits->Add(fitSig);
  listFits->Add(hFlowMass);
  listFits->Add(fitVn);
  listFits->Add(fitFlowBg);
  listFits->Add(fitFlowSig);

  // Reseting the canvas (removing drawn things)
  canFitInvMass->Clear();
  canFitInvMass->Divide(2,1);

  TLatex* latex = new TLatex();
  // latex->SetLineColor(kRed);
  latex->SetNDC();

  canFitInvMass->cd(1);
  // gPad->SetLogy();
  hInvMass->SetStats(0);
  hInvMass->SetMinimum(0);
  hInvMass->GetXaxis()->SetTitle("M_{#Lambda} (GeV/c^{2})");
  hInvMass->SetMarkerStyle(20);
  hInvMass->DrawCopy();
  fitMass->DrawCopy("same");
  fitSig->DrawCopy("same");
  fitBg->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitMass->GetChisquare(), fitMass->GetNDF(),fitMass->GetChisquare()/fitMass->GetNDF()));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#mu = %.6f #pm %.6f}",fitMass->GetParameter(5),fitMass->GetParError(5)));
  latex->DrawLatex(0.17,0.70,Form("#color[9]{#sigma_{1} = %.6f #pm %.6f}",fitMass->GetParameter(6),fitMass->GetParError(6)));
  latex->DrawLatex(0.17,0.65,Form("#color[9]{#sigma_{2} = %.6f #pm %.6f}",fitMass->GetParameter(8),fitMass->GetParError(8)));

  canFitInvMass->cd(2);
  hFlowMass->GetXaxis()->SetTitle("M_{K^{0}} (GeV/c^{2})");
  hFlowMass->SetStats(0);
  hFlowMass->DrawCopy();
  fitVn->DrawCopy("same");
  latex->DrawLatex(0.17,0.80,Form("#color[9]{v_{2} = %.4f #pm %.4f}",dFlow,dFlowError));
  latex->DrawLatex(0.17,0.75,Form("#color[9]{#chi^{2}/ndf = %.3g/%d = %.3g}",fitVn->GetChisquare(), fitVn->GetNDF(),fitVn->GetChisquare()/fitVn->GetNDF()));

  return kTRUE;
}
//_____________________________________________________________________________
void ProcessUniFlow::Fatal(TString sMsg, TString sMethod)
{
	printf("\033[91mFatal::%s  %s. Terminating!\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Error(TString sMsg, TString sMethod)
{
	printf("\033[91mError::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Info(TString sMsg, TString sMethod)
{
	printf("\033[96mInfo::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Warning(TString sMsg, TString sMethod)
{
	printf("\033[93mWarning::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Debug(TString sMsg, TString sMethod)
{
	if(fbDebug) printf("\033[95mDebug::%s  %s\033[0m\n", sMethod.Data(), sMsg.Data());
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
