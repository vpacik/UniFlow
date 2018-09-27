/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUniFlow_H
#define AliAnalysisTaskUniFlow_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPicoTrack.h"

struct FlowTask
{
  Bool_t                fbDoRefs; // which particles are procesed (RFPs / POIs / both )
  Bool_t                fbDoPOIs; // which particles are procesed (RFPs / POIs / both )
  Int_t                 fiNumHarm; // correlation order <M>
  Int_t                 fiNumGaps; // number of subevents
  std::vector<Int_t>    fiHarm; // harmonics n1,n2,...,nM
  std::vector<Double_t> fdGaps; // gaps between subevents (standard GF notation)
  TString               fsName; // automatically generated name: see Init() for format
  TString               fsLabel; // automatically generated label see Init() for format


              FlowTask() : fbDoRefs(kTRUE), fbDoPOIs(kTRUE), fiNumHarm(0), fiNumGaps(0), fsName("DummyName"), fsLabel("DummyLabel") {}; // default ctor
              FlowTask(Bool_t refs, Bool_t pois, std::vector<Int_t> harm, std::vector<Double_t> gaps = std::vector<Double_t>()) { fbDoRefs = refs; fbDoPOIs = pois; fiHarm = harm; fdGaps = gaps; fiNumHarm = harm.size(); fiNumGaps = gaps.size(); Init(); } // actual ctor
              ~FlowTask() { fiHarm.clear(); fdGaps.clear(); }
  void        Init(); // initialization
  void        Print(); // print FlowTask properties
  Bool_t      HasGap() { return (Bool_t) fiNumGaps; }; // check if Gap
};
//_____________________________________________________________________________
void FlowTask::Init()
{
  // Initilization of FlowTask

  if(fiNumHarm < 2)
  {
    fsName = "NA";
    fsLabel = "NA";
    return;
  }

  // generating name
  TString sName = Form("<<%d>>(%d",fiNumHarm,fiHarm[0]);
  for(Int_t i(1); i < fiNumHarm; ++i) { sName += Form(",%d",fiHarm[i]); }
  sName += ")";

  if(fiNumGaps > 0)
  {
    sName += Form("_%dsub(%.2g",fiNumGaps+1,fdGaps[0]);
    for(Int_t i(1); i < fiNumGaps; ++i) { sName += Form(",%.2g",fdGaps[i]); }
    sName += ")";
  }

  // generating label
  TString sLabel = Form("<<%d>>_{%d",fiNumHarm,fiHarm[0]);
  for(Int_t i(1); i < fiNumHarm; ++i) { sLabel += Form(",%d",fiHarm[i]); }
  sLabel += "}";

  if(fiNumGaps > 0)
  {
    sLabel += Form(" %dsub(|#Delta#eta| > %.2g",fiNumGaps+1,fdGaps[0]);
    for(Int_t i(1); i < fiNumGaps; ++i) { sLabel += Form(", |#Delta#eta| > %.2g",fdGaps[i]); }
    sLabel += ")";
  }

  fsName = sName;
  fsLabel = sLabel;
}
//_____________________________________________________________________________
void FlowTask::Print()
{
  printf("FlowTask::Print() : '%s' (%s) | fbDoRefs %d | fbDoPOIs %d | fiHarm[%d] = { ",fsName.Data(), fsLabel.Data(), fbDoRefs, fbDoPOIs, fiNumHarm);
  for(Int_t i(0); i < fiNumHarm; ++i) { printf("%d ",fiHarm[i]); }
  printf("} | fgGaps[%d] = { ",fiNumGaps);
  for(Int_t i(0); i < fiNumGaps; ++i) { printf("%0.2f ",fdGaps[i]); }
  printf("}\n");
}
//_____________________________________________________________________________


class AliAnalysisTaskUniFlow : public AliAnalysisTaskSE
{
    public:
      enum    RunMode {kFull = 0, kTest, kSkipFlow}; // task running mode (NOT GRID MODE)
      enum    ColSystem {kPP = 0, kPPb, kPbPb}; // tag for collisional system
      enum    AnalType {kAOD = 0, kESD}; // tag for analysis type
      enum    MultiEst {kRFP = 0, kV0A, kV0C, kV0M, kCL0, kCL1, kZNA, kZNC}; // multiplicity estimator as AliMultSelection
      enum    PartSpecies {kRefs = 0, kCharged, kPion, kKaon, kProton, kK0s, kLambda, kPhi, kUnknown}; // list of all particle species of interest; NB: kUknown last as counter
      enum    SparseCand {kInvMass = 0, kCent, kPt, kEta, kDim}; // reconstructed candidates dist. dimensions

                              AliAnalysisTaskUniFlow(); // constructor
                              AliAnalysisTaskUniFlow(const char *name); // named (primary) constructor
      virtual                 ~AliAnalysisTaskUniFlow(); // destructor

      virtual void            UserCreateOutputObjects(); //
      virtual void            UserExec(Option_t* option); // main methond - called for each event
      virtual void            Terminate(Option_t* option); // called after all events are processed
      // analysis setters
      void                    SetRunMode(RunMode mode = kFull) { fRunMode = mode; }
      void                    SetNumEventsAnalyse(Short_t num) { fNumEventsAnalyse = num; }
      void					          SetAnalysisType(AnalType type = kAOD) { fAnalType = type; }
      void                    SetMC(Bool_t mc = kTRUE) { fMC = mc; }
      void                    SetSampling(Bool_t sample = kTRUE) { fSampling = sample; }
      void                    SetFillQAhistos(Bool_t fill = kTRUE) { fFillQA = fill; }
      //void                    SetNumberOfSamples(Short_t numSamples = 10) { fNumSamples = numSamples; } // not implemented yet
      void                    SetProcessPID(Bool_t use = kTRUE) { fProcessSpec[kPion] = use; fProcessSpec[kKaon] = use; fProcessSpec[kProton] = use; }
      void                    SetProcessV0s(Bool_t use = kTRUE) { fProcessSpec[kK0s] = use; fProcessSpec[kLambda] = use; }
      void                    SetProcessPhi(Bool_t use = kTRUE) { fProcessSpec[kPhi] = use; }
      // flow related setters
      void                    AddTwo(Int_t n1, Int_t n2, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecFlowTask.push_back(new FlowTask(refs, pois, {n1,n2})); }
      void                    AddTwoGap(Int_t n1, Int_t n2, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecFlowTask.push_back(new FlowTask(refs, pois, {n1,n2}, {gap})); }
      void                    AddThree(Int_t n1, Int_t n2, Int_t n3, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecFlowTask.push_back(new FlowTask(refs, pois, {n1,n2,n3})); }
      void                    AddThreeGap(Int_t n1, Int_t n2, Int_t n3, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecFlowTask.push_back(new FlowTask(refs, pois, {n1,n2,n3} ,{gap})); }
      void                    AddFour(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecFlowTask.push_back(new FlowTask(refs, pois, {n1,n2,n3,n4})); }
      void                    AddFourGap(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Double_t gap, Bool_t refs = kTRUE, Bool_t pois = kTRUE) { fVecFlowTask.push_back(new FlowTask(refs, pois, {n1,n2,n3,n4}, {gap})); }

      void                    SetUseFixedMultBins(Bool_t fixed = kTRUE) { fUseFixedMultBins = fixed; }
      void                    SetFlowRFPsPtMin(Double_t pt) { fCutFlowRFPsPtMin = pt; }
      void                    SetFlowRFPsPtMax(Double_t pt) { fCutFlowRFPsPtMax = pt; }
      void                    SetFlowFillWeights(Bool_t weights = kTRUE) { fFlowFillWeights = weights; }
      void                    SetUseWeigthsFile(const char* file, Bool_t bRunByRun) { fFlowWeightsPath = file; fFlowRunByRunWeights = bRunByRun; fFlowUseWeights = kTRUE; } //! NOTE file has to include "alien:///" if the file is on grid
      void                    SetUseWeights3D(Bool_t use = kTRUE) { fFlowUse3Dweights = use; }
      // events setters
      void                    SetCollisionSystem(ColSystem colSystem = kPP) { fColSystem = colSystem; }
      void                    SetMultEstimator(MultiEst est) { fMultEstimator = est; }
      void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
      void					          SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
      void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
      // track setters
      void                    SetChargedEtaMax(Double_t eta) { fCutChargedEtaMax = eta; }
      void                    SetChargedDCAzMax(Double_t dcaz) {  fCutChargedDCAzMax = dcaz; }
      void                    SetChargedDCAxyMax(Double_t dcaxy) {  fCutChargedDCAxyMax = dcaxy; }
      void                    SetChargedNumTPCclsMin(UShort_t tpcCls) { fCutChargedNumTPCclsMin = tpcCls; }
      void                    SetChargedTrackFilterBit(UInt_t filter) { fCutChargedTrackFilterBit = filter; }
      // PID (pi,K,p) setters
      void                    SetPIDUseAntiProtonOnly(Bool_t use = kTRUE) { fCutPIDUseAntiProtonOnly = use; }
      void                    SetPIDNumSigmasPionMax(Float_t numSigmas) { fCutPIDnSigmaPionMax = numSigmas; }
      void                    SetPIDNumSigmasKaonMax(Float_t numSigmas) { fCutPIDnSigmaKaonMax = numSigmas; }
      void                    SetPIDNumSigmasProtonMax(Float_t numSigmas) { fCutPIDnSigmaProtonMax = numSigmas; }
      void                    SetPIDNumSigmasTPCRejectElectron(Float_t numSigmas) { fCutPIDnSigmaTPCRejectElectron = numSigmas; }
      void                    SetPIDNumSigmasCombinedNoTOFrejection(Bool_t reject = kTRUE) { fCutPIDnSigmaCombinedTOFrejection = reject; }
      void                    SetUseBayesPID(Bool_t bayes = kTRUE) { fCutUseBayesPID = bayes; }
      void                    SetPIDBayesProbPionMin(Double_t probPi) { fCutPIDBayesPionMin = probPi; }
      void                    SetPIDBayesProbKaonMin(Double_t probK) { fCutPIDBayesKaonMin = probK; }
      void                    SetPIDBayesProbProtonMin(Double_t probP) { fCutPIDBayesProtonMin = probP; }
      void                    SetPIDBayesRejectElectron(Double_t prob) { fCutPIDBayesRejectElectron = prob; }
      void                    SetPIDBayesRejectMuon(Double_t prob) { fCutPIDBayesRejectMuon = prob; }
      // V0s setters
      void					          SetV0sOnFly(Bool_t onFly) { fCutV0sOnFly = onFly; }
      void					          SetV0sTPCRefit(Bool_t refit) { fCutV0srefitTPC = refit; }
      void					          SetV0sRejectKinks(Bool_t reject) { fCutV0srejectKinks = reject; }
      void                    SetV0sDaughterNumTPCClsMin(UShort_t cls) { fCutV0sDaughterNumTPCClsMin = cls; }
      void                    SetV0sDaughterNumTPCrossMin(UShort_t cls) { fCutV0sDaughterNumTPCCrossMin = cls; }
      void                    SetV0sDaughterNumTPCFindMin(UShort_t cls) { fCutV0sDaughterNumTPCFindMin = cls; }
      void                    SetV0sDaughterNumTPCClsPIDMin(UShort_t cls) { fCutV0sDaughterNumTPCClsPIDMin = cls; }
      void                    SetV0sDaughterRatioCrossFindMin(Double_t ratio) { fCutV0sDaughterRatioCrossFindMin = ratio; }
      void					          SetV0sUseCrossMassRejection(Bool_t reject) { fCutV0sCrossMassRejection = reject; }
      void					          SetV0sCrossMassCutK0s(Double_t mass) { fCutV0sCrossMassCutK0s = mass; }
      void					          SetV0sCrossMassCutLambda(Double_t mass) { fCutV0sCrossMassCutLambda = mass; }
      void					          SetV0sDCAPVMin(Double_t dca) { fCutV0sDCAtoPVMin = dca; }
      void					          SetV0sDCAPVMax(Double_t dca) { fCutV0sDCAtoPVMax = dca; }
      void					          SetV0sDCAPVzMax(Double_t dca) { fCutV0sDCAtoPVzMax = dca; }
      void                    SetV0sDaughtersFilterBit(UInt_t filter) { fCutV0sDaughterFilterBit = filter; }
      void					          SetV0sDCADaughtersMin(Double_t dca) { fCutV0sDCADaughtersMin = dca; }
      void					          SetV0sDCADaughtersMax(Double_t dca) { fCutV0sDCADaughtersMax = dca; }
      void					          SetV0sDecayRadiusMin(Double_t radius) { fCutV0sDecayRadiusMin = radius; }
      void					          SetV0sDecayRadiusMax(Double_t radius) { fCutV0sDecayRadiusMax = radius; }
      void					          SetV0sDaughterEtaMax(Double_t eta) { fCutV0sDaughterEtaMax = eta; }
      void					          SetV0sDaughterPtMin(Double_t pt) { fCutV0sDaughterPtMin = pt; }
      void					          SetV0sDaughterPtMax(Double_t pt) { fCutV0sDaughterPtMax = pt; }
      void					          SetV0sMotherEtaMax(Double_t eta) { fCutV0sMotherEtaMax = eta; }
      void                    SetV0sMotherRapMax(Double_t rap) { fCutV0sMotherRapMax = rap; }
      void					          SetV0sK0sInvMassMin(Double_t mass) { fCutV0sInvMassK0sMin = mass; }
      void					          SetV0sK0sInvMassMax(Double_t mass) { fCutV0sInvMassK0sMax = mass; }
      void					          SetV0sLambdaInvMassMin(Double_t mass) { fCutV0sInvMassLambdaMin = mass; }
      void					          SetV0sLambdaInvMassMax(Double_t mass) { fCutV0sInvMassLambdaMax = mass; }
      void					          SetV0sK0sCPAMin(Double_t cpa) { fCutV0sCPAK0sMin = cpa; }
      void					          SetV0sLambdaCPAMin(Double_t cpa) { fCutV0sCPALambdaMin = cpa; }
      void					          SetV0sK0sNumTauMax(Double_t nTau) { fCutV0sNumTauK0sMax = nTau; }
      void					          SetV0sLambdaNumTauMax(Double_t nTau) { fCutV0sNumTauLambdaMax = nTau; }
      void					          SetV0sK0sArmenterosAlphaMin(Double_t alpha) { fCutV0sArmenterosAlphaK0sMin = alpha; }
      void					          SetV0sLambdaArmenterosAlphaMax(Double_t alpha) { fCutV0sArmenterosAlphaLambdaMax = alpha; }
      void                    SetV0sK0sPionNumTPCSigmaMax(Float_t nSigma) { fCutV0sK0sPionNumTPCSigmaMax = nSigma; }
      void                    SetV0sLambdaPionNumTPCSigmaMax(Float_t nSigma) { fCutV0sLambdaPionNumTPCSigmaMax = nSigma; }
      void                    SetV0sLambdaProtonNumTPCSigmaMax(Float_t nSigma) { fCutV0sLambdaProtonNumTPCSigmaMax = nSigma; }
      // phi setters
      void					          SetPhiMotherEtaMax(Double_t eta) { fCutPhiMotherEtaMax = eta; }
      void					          SetPhiInvMassMin(Double_t mass) { fCutPhiInvMassMin = mass; }
      void					          SetPhiInvMassMax(Double_t mass) { fCutPhiInvMassMax = mass; }


      AliEventCuts fEventCuts; //

    private:
      std::vector<FlowTask*>  fVecFlowTask; //
      // array lenghts & constants
      const Double_t          fPDGMassPion; // [DPGMass] DPG mass of charged pion
      const Double_t          fPDGMassKaon; // [DPGMass] DPG mass of charged kaon
      const Double_t          fPDGMassProton; // [DPGMass] DPG mass of proton
      const Double_t          fPDGMassPhi; // [DPGMass] DPG mass of phi (333) meson
      const Double_t          fPDGMassK0s; // [DPGMass] DPG mass of K0s
      const Double_t          fPDGMassLambda; // [DPGMass] DPG mass of (Anti)Lambda
      static const Short_t    fFlowNumHarmonicsMax = 7; // maximum harmonics length of flow vector array
      static const Short_t    fFlowNumWeightPowersMax = 5; // maximum weight power length of flow vector array
      const Double_t          fFlowPOIsPtMin; // [0] (GeV/c) min pT treshold for POIs for differential flow
      const Double_t          fFlowPOIsPtMax; // [15] (GeV/c) max pT treshold for POIs for differential flow
      Int_t                   fFlowCentMin; // [set in InitializeTask()] min range for centrality/multiplicity histos
      Int_t                   fFlowCentMax; // [set in InitializeTask()] max range for centrality/multiplicity histos
      static const Short_t    fV0sNumBinsMass = 60; // number of InvMass bins for V0s distribution
      static const Short_t    fPhiNumBinsMass = 60; // number of InvMass bins for phi distribution
      static const Short_t    fiNumIndexQA = 2; // QA indexes: 0: before cuts // 1: after cuts

      const static Short_t    fNumSamples = 10; // overall number of samples (from random sampling) used
      const static Int_t      fNumMultBins = 6; // number of multiplicity bins
      static Double_t         fMultBins[fNumMultBins+1]; // multiplicity bins

      const char*             GetSpeciesName(PartSpecies species);
      const char*             GetSpeciesLabel(PartSpecies species);
      const char*             GetEtaGapName(Double_t dEtaGap) { return Form("%02.2g",10.0*dEtaGap); }

      Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
      void                    ListParameters(); // list all task parameters
      void                    ClearVectors(); // properly clear all particle vectors

      Bool_t                  IsEventSelected(); // event selection for Run 2 using AliEventCuts
      Bool_t                  IsEventSelected_oldsmall2016(); // (old/manual) event selection for LHC2016 pp & pPb data
      Bool_t                  IsEventSelected_Naghmeh(); //
      Bool_t                  IsEventRejectedAddPileUp(); // additional pile-up rejection for Run2 Pb-Pb
      Bool_t                  LoadWeights(Bool_t init = kFALSE); // load weights histograms
      void                    FillEventsQA(const Short_t iQAindex); // filling QA plots related to event selection
      Short_t                 GetSamplingIndex(); // returns sampling index based on sampling selection (number of samples)
      Short_t                 GetCentralityIndex(); // returns centrality index based centrality estimator or number of selected tracks
      const char*             GetMultiEstimatorLabel(MultiEst est); // returns mult/cent estimator string with label or 'n/a' if not available

      void                    CalculateCorrelations(FlowTask* task, PartSpecies species, Double_t dPt = -1.0, Double_t dMass = -1.0); // wrapper for correlations methods
      Bool_t                  ProcessFlowTask(FlowTask* task); // procesisng of FlowTask
      Bool_t                  CalculateFlow(); // main (envelope) method for flow calculations in selected events

      void                    FilterCharged(); // charged tracks filtering
      void                    FilterPID(); // pi,K,p filtering
      void                    FilterV0s(); // K0s, Lambda, ALambda filtering
      void                    FilterPhi(); // reconstruction and filtering of Phi meson candidates
      AliAODMCParticle*       GetMCParticle(Int_t label); // find corresponding MC particle from fArrayMC depending of AOD track label
      Double_t                GetRapidity(Double_t mass, Double_t Pt, Double_t Eta); // calculate particle / track rapidity
      Bool_t                  HasTrackPIDTPC(const AliAODTrack* track); // is TPC PID OK for this track ?
      Bool_t                  HasTrackPIDTOF(const AliAODTrack* track); // is TOF PID OK for this track ?
      Bool_t                  IsWithinRefs(const AliAODTrack* track); // check if track fulfill requirements for Refs (used for refs selection & autocorelations)
      Bool_t                  IsChargedSelected(const AliAODTrack* track = 0x0); // charged track selection
      PartSpecies             IsPIDSelected(const AliAODTrack* track); // PID tracks selections
      Bool_t                  IsV0Selected(const AliAODv0* v0 = 0x0); // general (common) V0 selection
      Bool_t                  IsV0aK0s(const AliAODv0* v0 = 0x0); // V0 selection: K0s specific
      Short_t                 IsV0aLambda(const AliAODv0* v0 = 0x0); // V0 selection: (A)Lambda specific
      AliPicoTrack*           MakeMother(const AliAODTrack* part1, const AliAODTrack* part2); // Combine two prongs into a mother particle stored in AliPicoTrack object
      void                    FillSparseCand(THnSparse* sparse, AliVTrack* track); // Fill sparse histogram for inv. mass distribution of candidates (V0s,Phi)
      void                    FillQARefs(const Short_t iQAindex, const AliAODTrack* track = 0x0); // filling QA plots for RFPs selection
      void                    FillQACharged(const Short_t iQAindex, const AliAODTrack* track = 0x0); // filling QA plots for charged track selection
      void                    FillQAPID(const Short_t iQAindex, const AliAODTrack* track = 0x0, const PartSpecies species = kUnknown); // filling pi,K,p QA histograms
      void                    FillQAV0s(const Short_t iQAindex, const AliAODv0* v0 = 0x0, const Bool_t bIsK0s = kTRUE, const Short_t bIsLambda = 2); // filling QA plots for V0s candidates
      void                    FillQAPhi(const Short_t iQAindex, const AliPicoTrack* part = 0x0); // filling QA plots for V0s candidates

      // Flow related methods
      void                    FillRefsVectors(const Double_t dGap); // fill flow vector Q with RFPs for reference flow
      void                    FillPOIsVectors(const Double_t dEtaGap, const PartSpecies species, const Double_t dPtLow, const Double_t dPtHigh, const Double_t dMassLow = 0.0, const Double_t dMassHigh = 0.0); // fill flow vectors p,q and s with POIs (for given species) for differential flow calculations
      void                    ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]); // set values to TComplex(0,0,0) for given array
      void                    ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]); // printf all values of given Flow vector array

      TComplex                Q(const Short_t n, const Short_t p);
      TComplex                QGapPos(const Short_t n, const Short_t p);
      TComplex                QGapNeg(const Short_t n, const Short_t p);
      TComplex                QGapMid(const Short_t n, const Short_t p);
      TComplex                P(const Short_t n, const Short_t p);
      TComplex                PGapPos(const Short_t n, const Short_t p);
      TComplex                PGapNeg(const Short_t n, const Short_t p);
      TComplex                S(const Short_t n, const Short_t p);

      TComplex                Two(const Short_t n1, const Short_t n2); // Two particle reference correlation calculations (no eta gap)
      TComplex                TwoGap(const Short_t n1, const Short_t n2); // Two particle reference correlation calculations (with eta gap)
      TComplex                Three(const Short_t n1, const Short_t n2, const Short_t n3); // Three particle reference correlation calculations (no eta gap)
      TComplex                Four(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (no eta gap)
      TComplex                FourGap(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (no eta gap)
      TComplex                Four3sub(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (with 3 sub-events)

      TComplex                TwoDiff(const Short_t n1, const Short_t n2); // Two particle diff. correlation calculations (no eta gap)
      TComplex                TwoDiffGapPos(const Short_t n1, const Short_t n2); // Two particle diff. correlation calculations (with eta gap)
      TComplex                TwoDiffGapNeg(const Short_t n1, const Short_t n2); // Two particle diff. correlation calculations (with eta gap)
      TComplex                ThreeDiff(const Short_t n1, const Short_t n2, const Short_t n3); // Three particle diff. correlation calculations (no eta gap)
      TComplex                ThreeDiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3); // Three particle diff. correlation calculations (with eta gap)
      TComplex                ThreeDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3); // Three particle diff. correlation calculations (with eta gap)
      TComplex                FourDiff(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (no eta gap)
      TComplex                FourDiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (with eta gap)
      TComplex                FourDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (with eta gap)

      // properties
      AliAODEvent*            fEventAOD; //! AOD event countainer
      Double_t                fPVz; // PV z-coordinate used for weights
      AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
      AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
      TFile*                  fFlowWeightsFile; //! source file containing weights
      TClonesArray*           fArrayMC; //! input list of MC particles
      Bool_t                  fMC; // is running on mc?
      Bool_t                  fInit; // initialization check
      Short_t                 fIndexSampling; // sampling index (randomly generated)
      Short_t                 fIndexCentrality; // centrality bin index (based on centrality est. or number of selected tracks)
      Short_t                 fEventCounter; // event counter (used for local test runmode purpose)
      Short_t                 fNumEventsAnalyse; // [50] number of events to be analysed / after passing selection (only in test mode)
      Int_t                   fRunNumber; // [-1] run number of previous event (not the current one)

      TComplex                fFlowVecQpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecQneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecQmid[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecPpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecPneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecS[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation

      // selected POIs containers
      std::vector<AliVTrack*>*  fVector[kUnknown]; //! container for selected Refs charged particles

      //cuts & selection: analysis
      RunMode                 fRunMode; // running mode (not grid related)
      AnalType                fAnalType; // analysis type: AOD / ESD
      Bool_t                  fSampling;      // Do random sampling ? (estimation of vn stat. uncertanity)
      Bool_t                  fFillQA; //[kTRUE] flag for filling the QA plots
      Bool_t                  fProcessSpec[kUnknown];  // flag for processing species
      // cuts & selection: flow related
      Bool_t                  fUseFixedMultBins; // [kFALSE] setting fixed multiplicity bins
      Double_t                fCutFlowRFPsPtMin; // [0] (GeV/c) min pT treshold for RFPs particle for reference flow
      Double_t                fCutFlowRFPsPtMax; // [0] (GeV/c) max pT treshold for RFPs particle for reference flow
      Bool_t                  fFlowFillWeights; //[kFALSE] flag for filling weights
      Bool_t                  fFlowUseWeights; //[kFALSE] flag for using the previously filled weights (NOTE: this is turned on only when path to file is applied via fFlowWeightsPath)
      Bool_t                  fFlowUse3Dweights; // [kFALSE] flag for using 3D GF weights, if kFALSE, 2D weights are expected
      Bool_t                  fFlowRunByRunWeights; // [kTRUE] flag for using rub-by-run weigths from weigths file; if false, only one set of histrograms is provided
      TString                 fFlowWeightsPath; //[] path to source root file with weigthts (if empty unit weights are applied) e.g. "alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_LHC16kl.root"

      //cuts & selection: events
      ColSystem               fColSystem; // collisional system
      AliVEvent::EOfflineTriggerTypes    fTrigger; // physics selection trigger
      MultiEst                fMultEstimator; // multiplicity/centrality estimator as in AliMultSelection
      Double_t                fPVtxCutZ; // (cm) PV z cut
      Bool_t                  fEventRejectAddPileUp; // additional pile-up rejection for Pb-Pb collisions in Run2 (17n, 15o)
      //cuts & selection: tracks
      UInt_t                  fCutChargedTrackFilterBit; // (-) tracks filter bit
      UShort_t                fCutChargedNumTPCclsMin;  // (-) Minimal number of TPC clusters used for track reconstruction
      Double_t                fCutChargedEtaMax; // (-) Maximum pseudorapidity range
      Double_t                fCutChargedDCAzMax; // (cm) Maximal DCA-z cuts for tracks (pile-up rejection suggested for LHC16)
      Double_t                fCutChargedDCAxyMax; // (cm) Maximal DCA-xy cuts for tracks (pile-up rejection suggested for LHC16)
      // cuts & selection: PID selection
      Bool_t                  fCutPIDUseAntiProtonOnly; // [kFALSE] check proton PID charge to select AntiProtons only
      Bool_t                  fCutPIDnSigmaCombinedTOFrejection; // [kTRUE] flag for rejection candidates in TPC+TOF pt region if TOF is not available (if true and no TOF track is skipped, otherwise only TPC is used)
      Float_t                 fCutPIDnSigmaPionMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for pion candidates
      Float_t                 fCutPIDnSigmaKaonMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for kaon candidates
      Float_t                 fCutPIDnSigmaProtonMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for proton candidates
      Float_t                 fCutPIDnSigmaTPCRejectElectron; // [3] number of TPC nSigma for electron rejection
      Bool_t                  fCutUseBayesPID; // [kFALSE] flag for using Bayes PID for pi,K,p instead nsigma cut
      Double_t                fCutPIDBayesPionMin; // [0.9] minimal value of Bayes PID probability for pion
      Double_t                fCutPIDBayesKaonMin; // [0.9] minimal value of Bayes PID probability for Kaon
      Double_t                fCutPIDBayesProtonMin; // [0.9] minimal value of Bayes PID probability for proton
      Double_t                fCutPIDBayesRejectElectron; // [0.5] maximal value of Bayes PID probability for electron rejection
      Double_t                fCutPIDBayesRejectMuon; // [0.5] maximal value of Bayes PID probability for muon rejection
      //cuts & selection: V0 reconstruction
	    Bool_t 					        fCutV0sOnFly;		// V0 reconstruction method: is On-the-fly? (or offline)
  		Bool_t					        fCutV0srefitTPC; // Check TPC refit of V0 daughters ?
  		Bool_t					        fCutV0srejectKinks; // Reject Kink V0 daughter tracks ?
      UShort_t                fCutV0sDaughterNumTPCClsMin; // min number of TPC clusters
      UShort_t                fCutV0sDaughterNumTPCCrossMin; // min number of crossed TPC rows
      UShort_t                fCutV0sDaughterNumTPCFindMin; // min number of findable TPC clusters
      UShort_t                fCutV0sDaughterNumTPCClsPIDMin; // min number of TPC clusters used for PID
      Double_t                fCutV0sDaughterRatioCrossFindMin; // min ratio of crossed / findable TPC clusters
  		Bool_t					        fCutV0sCrossMassRejection; // competing V0 rejection based on InvMass
      Double_t                fCutV0sCrossMassCutK0s; // [0.005] (GeV/c2) restricted vicinity of Lambda/ALambda inv. mass peak for K0s candidates
      Double_t                fCutV0sCrossMassCutLambda; // [0.020] (GeV/c2) restricted vicinity of K0s inv. mass peak for Lambda/ALambda candidates
  		Double_t                fCutV0sDCAtoPVMin;   // (cm) min DCA of V0 daughter to PV
      Double_t				        fCutV0sDCAtoPVMax;	// (cm) max DCA of V0 daughter to PV
      Double_t                fCutV0sDCAtoPVzMax; // (cm) max DCA-z coordinate of V0 daughters to PV
		  Double_t				        fCutV0sDCADaughtersMin;	// (cm) min DCA of V0 daughters among themselves
		  Double_t				        fCutV0sDCADaughtersMax;	// (cm) max DCA of V0 daughters among themselves
      Double_t                fCutV0sDecayRadiusMin; // (cm) min distance of secondary vertex from z-axis in transverse plane
		  Double_t				        fCutV0sDecayRadiusMax; // (cm) max distance of secondary vertex from z-axis in transverse plane
      UInt_t                  fCutV0sDaughterFilterBit; // (-) V0 daughters filter bit
      Double_t                fCutV0sDaughterPtMin; // (GeV/c) min pT of V0 daughters
      Double_t                fCutV0sDaughterPtMax; // (GeV/c) max pT of V0 daughters
      Double_t                fCutV0sDaughterEtaMax; // (-) max value of Eta of V0 daughters
      Double_t                fCutV0sMotherEtaMax; // (-) max eta value of V0 mother
      Double_t                fCutV0sMotherRapMax; // (-) max rapidity value of V0 mother
      Double_t                fCutV0sCPAK0sMin;    // (-) min cosine of pointing angle of K0s candidate to PV
      Double_t                fCutV0sCPALambdaMin; // (-) min cosine of pointing angle of K0s candidate to PV
      Double_t                fCutV0sNumTauK0sMax; // (c*tau) max number of c*tau (K0s)
      Double_t                fCutV0sNumTauLambdaMax; // (c*tau) max number of c*tau ((A)Lambda)
      Double_t                fCutV0sInvMassK0sMin; // [0.4] (GeV/c2) min inv. mass window for selected K0s candidates
      Double_t                fCutV0sInvMassK0sMax; // [0.6] (GeV/c2) max inv. mass window for selected K0s candidates
      Double_t                fCutV0sInvMassLambdaMin; // [1.08] (GeV/c2) min inv. mass window for selected (Anti)Lambda candidates
      Double_t                fCutV0sInvMassLambdaMax; // [1.16] (GeV/c2) max inv. mass window for selected (Anti)Lambda candidates
      Double_t				        fCutV0sArmenterosAlphaK0sMin; // (alpha) min Armenteros alpha for K0s
      Double_t                fCutV0sArmenterosAlphaLambdaMax; // (alpha) max Armenteros alpha for (Anti)Lambda
      Float_t                 fCutV0sK0sPionNumTPCSigmaMax; // (sigmaTPC) max number of TPC sigmas for kaon PID (K0s candidates)
      Float_t                 fCutV0sLambdaPionNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for pion PID (Lambda candidates)
      Float_t                 fCutV0sLambdaProtonNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for proton PID (Lambda candidates)
      // cuts & selection: phi
      Double_t                fCutPhiMotherEtaMax; // (-) max value of phi candidate pseudorapidity
      Double_t                fCutPhiInvMassMin; // [0.99] (GeV/c2) min inv. mass window for selected phi candidates
      Double_t                fCutPhiInvMassMax; // [1.07] (GeV/c2) min inv. mass window for selected phi candidates

      // output lists
      TList*                  fQAEvents; //! events list
      TList*                  fQACharged; //! charged tracks list
      TList*                  fQAPID; //! pi,K,p list
      TList*                  fQAV0s; //! V0s candidates list
      TList*                  fQAPhi; //! Phi candidates list
      TList*                  fFlowWeights; //! list for flow weights
      TList*                  fListFlow[kUnknown]; //! flow lists


      // histograms & profiles

      // Flow
      THnSparseD*     fhsV0sCandK0s; //! distribution of K0s candidates
      THnSparseD*     fhsV0sCandLambda; //!  distribution of Lambda candidates
      THnSparseD*     fhsPhiCandSig; //!  distribution of Phi candidates
      THnSparseD*     fhsPhiCandBg; //!  distribution of Phi background

      TH2D*           fh2Weights[kUnknown]; //! container for GF weights (phi,eta,pt) (2D)
      TH3D*           fh3Weights[kUnknown]; //! container for GF weights (phi,eta,pt)
      TH3D*           fh3AfterWeights[kUnknown]; //! distribution after applying GF weights (phi,eta,pt)

      // Events
      TH2D*           fhEventSampling; //! distribution of sampled events (based on randomly generated numbers)
      TH1D*           fhEventCentrality; //! distribution of event centrality
      TH2D*           fh2EventCentralityNumRefs; //! distribution of event centrality vs number of selected charged tracks
      TH1D*           fhEventCounter; //! counter following event selection
      // Charged
      TH1D*           fhRefsMult; //!multiplicity distribution of selected RFPs
      TH1D*           fhRefsPt; //! pt distribution of selected RFPs
      TH1D*           fhRefsEta; //! pt distribution of selected RFPs
      TH1D*           fhRefsPhi; //! pt distribution of selected RFPs
      TProfile*       fpRefsMult; //! <multiplicity>
      TH1D*           fhChargedCounter; //! counter following charged track selection
      // PID
      TH1D*           fhPIDPionMult; //! multiplicity distribution of selected pions
      TH1D*           fhPIDPionPt; //! pt distribution of selected pions
      TH1D*           fhPIDPionPhi; //! phi distribution of selected pions
      TH1D*           fhPIDPionEta; //! eta distribution of selected pions
      TH1D*           fhPIDPionCharge; //! charge distribution of selected pions
      TH2D*           fh2PIDPionTPCdEdx; //! TPC dEdx response of selected pions
      TH2D*           fh2PIDPionTOFbeta; //! TOF beta of selected pions
      TH2D*           fh2PIDPionTPCnSigmaPion; //! TPC nSigma vs pT for selected pions (pion hypothesis)
      TH2D*           fh2PIDPionTOFnSigmaPion; //! TOF nSigma vs pT for selected pions (pion hypothesis)
      TH2D*           fh2PIDPionBayesPion; //! Bayesian PID probability vs pT for selected pions (pion hypothesis)
      TH2D*           fh2PIDPionTPCnSigmaKaon; //! TPC nSigma vs pT for selected pions (kaon hypothesis)
      TH2D*           fh2PIDPionTOFnSigmaKaon; //! TOF nSigma vs pT for selected pions (kaon hypothesis)
      TH2D*           fh2PIDPionBayesKaon; //! Bayesian PID probability vs pT for selected pions (kaon hypothesis)
      TH2D*           fh2PIDPionTPCnSigmaProton; //! TPC nSigma vs pT for selected pions (proton hypothesis)
      TH2D*           fh2PIDPionTOFnSigmaProton; //! TOF nSigma vs pT for selected pions (proton hypothesis)
      TH2D*           fh2PIDPionBayesProton; //! Bayesian PID probability vs pT for selected pions (proton hypothesis)
      TH1D*           fhPIDKaonMult; //! multiplicity distribution of selected pions
      TH1D*           fhPIDKaonPt; //! pt distribution of selected kaons
      TH1D*           fhPIDKaonPhi; //! phi distribution of selected kaons
      TH1D*           fhPIDKaonEta; //! eta distribution of selected kaons
      TH1D*           fhPIDKaonCharge; //! charge distribution of selected pions
      TH2D*           fh2PIDKaonTPCdEdx; //! TPC dEdx response of selected pions
      TH2D*           fh2PIDKaonTOFbeta; //! TOF beta of selected pions
      TH2D*           fh2PIDKaonTPCnSigmaPion; //! TPC nSigma vs pT for selected kaons (pion hypothesis)
      TH2D*           fh2PIDKaonTOFnSigmaPion; //! TOF nSigma vs pT for selected kaons (pion hypothesis)
      TH2D*           fh2PIDKaonBayesPion; //! Bayesian PID probability vs pT for selected kaons (pion hypothesis)
      TH2D*           fh2PIDKaonTPCnSigmaKaon; //! TPC nSigma vs pT for selected kaons (kaon hypothesis)
      TH2D*           fh2PIDKaonTOFnSigmaKaon; //! TOF nSigma vs pT for selected kaons (kaon hypothesis)
      TH2D*           fh2PIDKaonBayesKaon; //! Bayesian PID probability vs pT for selected kaons (kaon hypothesis)
      TH2D*           fh2PIDKaonTPCnSigmaProton; //! TPC nSigma vs pT for selected kaons (proton hypothesis)
      TH2D*           fh2PIDKaonTOFnSigmaProton; //! TOF nSigma vs pT for selected kaons (proton hypothesis)
      TH2D*           fh2PIDKaonBayesProton; //! Bayesian PID probability vs pT for selected kaons (proton hypothesis)
      TH1D*           fhPIDProtonMult; //! multiplicity distribution of selected pions
      TH1D*           fhPIDProtonPt; //! pt distribution of selected protons
      TH1D*           fhPIDProtonPhi; //! phi distribution of selected protons
      TH1D*           fhPIDProtonEta; //! eta distribution of selected protons
      TH1D*           fhPIDProtonCharge; //! charge distribution of selected pions
      TH2D*           fh2PIDProtonTPCdEdx; //! TPC dEdx response of selected pions
      TH2D*           fh2PIDProtonTOFbeta; //! TOF beta of selected pions
      TH2D*           fh2PIDProtonTPCnSigmaPion; //! TPC nSigma vs pT for selected protons (pion hypothesis)
      TH2D*           fh2PIDProtonTOFnSigmaPion; //! TOF nSigma vs pT for selected protons (pion hypothesis)
      TH2D*           fh2PIDProtonBayesPion; //! Bayesian PID probability vs pT for selected protons (pion hypothesis)
      TH2D*           fh2PIDProtonTPCnSigmaKaon; //! TPC nSigma vs pT for selected protons (kaon hypothesis)
      TH2D*           fh2PIDProtonTOFnSigmaKaon; //! TOF nSigma vs pT for selected protons (kaon hypothesis)
      TH2D*           fh2PIDProtonBayesKaon; //! Bayesian PID probability vs pT for selected protons (kaon hypothesis)
      TH2D*           fh2PIDProtonTPCnSigmaProton; //! TPC nSigma vs pT for selected protons (proton hypothesis)
      TH2D*           fh2PIDProtonTOFnSigmaProton; //! TOF nSigma vs pT for selected protons (proton hypothesis)
      TH2D*           fh2PIDProtonBayesProton; //! Bayesian PID probability vs pT for selected protons (proton hypothesis)
      TH1D*           fhMCRecoSelectedPionPt; //! pt dist of selected (MC reco) pions
      TH1D*           fhMCRecoSelectedTruePionPt; //! pt dist of selected (MC reco) true (tagged in MC gen) pions
      TH1D*           fhMCRecoAllPionPt; //! pt dist of all (MC reco) pions (i.e. selected charged tracks that are tagged in MC)
      TH1D*           fhMCGenAllPionPt; //! pt dist of all (MC) generated pions
      TH1D*           fhMCRecoSelectedKaonPt; //! pt dist of selected (MC reco) Kaons
      TH1D*           fhMCRecoSelectedTrueKaonPt; //! pt dist of selected (MC reco) true (tagged in MC gen) Kaons
      TH1D*           fhMCRecoAllKaonPt; //! pt dist of all (MC reco) Kaons (i.e. selected charged tracks that are tagged in MC)
      TH1D*           fhMCGenAllKaonPt; //! pt dist of all (MC) generated Kaons
      TH1D*           fhMCRecoSelectedProtonPt; //! pt dist of selected (MC reco) Protons
      TH1D*           fhMCRecoSelectedTrueProtonPt; //! pt dist of selected (MC reco) true (tagged in MC gen) Protons
      TH1D*           fhMCRecoAllProtonPt; //! pt dist of all (MC reco) Protons (i.e. selected charged tracks that are tagged in MC)
      TH1D*           fhMCGenAllProtonPt; //! pt dist of all (MC) generated Protons
      // Phi
      TH1D*           fhPhiCounter; //! counter following phi candidate selection
      TH1D*           fhPhiMult; //! multiplicity distribution of selected phi candidates
      TH1D*           fhPhiBGMult; //! multiplicity distribution of BG candidates
      TH1D*           fhPhiInvMass; //! invariant mass distribution of phi candidates
      TH1D*           fhPhiBGInvMass; //! invariant mass distribution of phi background candidates
      TH1D*           fhPhiCharge; //! charge distribution of selected phi candidates
      TH1D*           fhPhiBGCharge; //! charge distribution of phi BG candidates
      TH1D*           fhPhiPt; //! pt distribution of selected phi candidates
      TH1D*           fhPhiEta; //! eta distribution of selected phi candidates
      TH1D*           fhPhiPhi; //! phi distribution of selected phi candidates
      // V0s
      TH1D*           fhV0sCounter; //! counter following V0s selection
      TH1D*           fhV0sCounterK0s; //! counter following K0s selection
      TH1D*           fhV0sCounterLambda; //! counter following (Anti-)Lambda selection
      TH2D*           fhV0sInvMassK0s; //! 2D inv. mass distiburion (K0s mass vs. Lambda/AntiLambda mass)
      TH2D*           fhV0sInvMassLambda; //! 2D inv. mass distiburion (K0s mass vs. Lambda/AntiLambda mass)
      TH2D*           fhV0sCompetingInvMassK0s; //! dist of InvMass of rejected K0s candidates in (Anti-)Lambda peak
      TH2D*           fhV0sCompetingInvMassLambda; //! dist of InvMass of rejected (Anti-)Lambda candidates in K0s peak


      // QA: events
      TH1D*           fhQAEventsPVz[fiNumIndexQA]; //!
      TH1D*           fhQAEventsNumContrPV[fiNumIndexQA]; //!
      TH1D*           fhQAEventsNumSPDContrPV[fiNumIndexQA]; //!
      TH1D*           fhQAEventsDistPVSPD[fiNumIndexQA]; //!
      TH1D*           fhQAEventsSPDresol[fiNumIndexQA]; //!
      TH2D*           fhQAEventsfMult32vsCentr;
      TH2D*           fhQAEventsMult128vsCentr;
      TH2D*           fhQAEventsfMultTPCvsTOF;
      TH2D*           fhQAEventsfMultTPCvsESD;
      // QA: charged tracks
      TH1D*           fhQAChargedMult[fiNumIndexQA];       //! number of AOD charged tracks distribution
      TH1D*           fhQAChargedPt[fiNumIndexQA];         //! pT dist of charged tracks
      TH1D*           fhQAChargedEta[fiNumIndexQA];        //! eta dist of charged tracks
      TH1D*           fhQAChargedPhi[fiNumIndexQA];        //! phi dist of charged tracks
      TH1D*           fhQAChargedCharge[fiNumIndexQA];     //! charge dist of charged tracks
      TH1D*           fhQAChargedFilterBit[fiNumIndexQA];  //! filter bit distribution of charged tracks
      TH1D*           fhQAChargedNumTPCcls[fiNumIndexQA];  //! dist of track number of TPC clusters
      TH1D*           fhQAChargedDCAxy[fiNumIndexQA];      //! dist of Charged DCA in transverse plane
      TH1D*           fhQAChargedDCAz[fiNumIndexQA];       //! dist of charged DCA in z coordinate
      // QA: PID tracks
      TH1D*           fhQAPIDTPCstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
      TH1D*           fhQAPIDTOFstatus[fiNumIndexQA];  //! based on AliPIDResponse::CheckPIDStatus();
      TH2D*           fhQAPIDTPCdEdx[fiNumIndexQA];    //! TPC PID information
      TH2D*           fhQAPIDTOFbeta[fiNumIndexQA];    //! TOF PID information
      TH3D*           fh3QAPIDnSigmaTPCTOFPtPion[fiNumIndexQA]; //! nSigma TPC vs nSigma TOF vs pt
      TH3D*           fh3QAPIDnSigmaTPCTOFPtKaon[fiNumIndexQA]; //! nSigma TPC vs nSigma TOF vs pt
      TH3D*           fh3QAPIDnSigmaTPCTOFPtProton[fiNumIndexQA]; //! nSigma TPC vs nSigma TOF vs pt
      TH3D*           fh3QAPIDnSigmaBayesElectron[fiNumIndexQA]; //! PID information (nSigma TPC, nSigma TOF, Bayes) for pions
      TH3D*           fh3QAPIDnSigmaBayesMuon[fiNumIndexQA]; //! PID information (nSigma TPC, nSigma TOF, Bayes) for pions
      TH3D*           fh3QAPIDnSigmaBayesPion[fiNumIndexQA]; //! PID information (nSigma TPC, nSigma TOF, Bayes) for pions
      TH3D*           fh3QAPIDnSigmaBayesKaon[fiNumIndexQA]; //! PID information (nSigma TPC, nSigma TOF, Bayes) for kaons
      TH3D*           fh3QAPIDnSigmaBayesProton[fiNumIndexQA]; //! PID information (nSigma TPC, nSigma TOF, Bayes) for proton
      // QA: V0s candidates
      TH1D*			  		fhQAV0sMultK0s[fiNumIndexQA];	//! number of K0s candidates
      TH1D*			  		fhQAV0sMultLambda[fiNumIndexQA];	//! number of Lambda candidates
      TH1D*			  		fhQAV0sMultALambda[fiNumIndexQA];	//! number of Anti-Lambda candidates
      TH1D*			  		fhQAV0sRecoMethod[fiNumIndexQA];	//! offline/online V0 reconstruction method
      TH1D*			  		fhQAV0sDaughterTPCRefit[fiNumIndexQA];	//! Daughters TPC refit true/false
      TH1D*			  		fhQAV0sDaughterKinks[fiNumIndexQA];	//! Daughters kinks true/false
      TH1D*           fhQAV0sDaughterNumTPCCls[fiNumIndexQA]; //! Daughter # of TPC findable clusters
      TH1D*           fhQAV0sDaughterNumTPCFind[fiNumIndexQA]; //! Daughter # of TPC clusters
      TH1D*           fhQAV0sDaughterNumTPCCrossRows[fiNumIndexQA]; //! Daughter # of TPC crossed rows
      TH1D*           fhQAV0sDaughterTPCCrossFindRatio[fiNumIndexQA]; //! Daughter # of TPC cross / # of TPC findable cls ratio
      TH1D*           fhQAV0sDaughterNumTPCClsPID[fiNumIndexQA]; //! Daughter # of TPC findable clusters used for PID
      TH1D*			  		fhQAV0sDCAtoPV[fiNumIndexQA];	//! V0 DCA to PV
      TH1D*			  		fhQAV0sDCADaughters[fiNumIndexQA];	//! DCA between V0 daughters
      TH1D*			  		fhQAV0sDecayRadius[fiNumIndexQA];	//! Distance between PV and Secondary vertex in transverse plane
      TH1D*           fhQAV0sInvMassK0s[fiNumIndexQA];    //! inv. mass dist of V0s (K0s mass hypothesis)
      TH1D*					  fhQAV0sInvMassLambda[fiNumIndexQA];	//! inv. mass dist of V0s ((A)Lambda mass hypothesis)
      TH1D*           fhQAV0sMotherPt[fiNumIndexQA];  //! pT dist of V0s
      TH1D*					  fhQAV0sMotherPhi[fiNumIndexQA];	//! azimuthal dist of V0s
      TH1D*           fhQAV0sMotherEta[fiNumIndexQA]; //! pseudorapidity dist of V0s
      TH1D*           fhQAV0sMotherCharge[fiNumIndexQA]; //! charge distribution of mothers
      TH1D*           fhQAV0sMotherRapK0s[fiNumIndexQA];  //! rapidity dist of V0s (K0s mass hypothesis)
      TH1D*           fhQAV0sMotherRapLambda[fiNumIndexQA]; //! rapidity dist of V0s (Lambda mass hypothesis)
      TH1D*           fhQAV0sDaughterPt[fiNumIndexQA];    //! pT dist of V0 daughters
      TH1D*					  fhQAV0sDaughterPhi[fiNumIndexQA];	//! pT dist of V0 daughters
      TH1D*           fhQAV0sDaughterEta[fiNumIndexQA];   //! pseudorapidity dist of V0 daughters
      TH1D*           fhQAV0sDaughterCharge[fiNumIndexQA]; //! charge distribution of daughters
      TH1D*					  fhQAV0sDaughterTPCstatus[fiNumIndexQA];	//! TPC dEdx vs p of K0s daughters
      TH1D*					  fhQAV0sDaughterTOFstatus[fiNumIndexQA];	//! TPC dEdx vs p of K0s daughters
      TH2D*					  fhQAV0sDaughterTPCdEdxK0s[fiNumIndexQA];	//! TPC dEdx vs p of K0s daughters
      TH2D*					  fhQAV0sDaughterNumSigmaPionK0s[fiNumIndexQA];	//! Number of TPC sigmas (pion) vs mother pT of K0s daughters
      TH2D*					  fhQAV0sDaughterTPCdEdxLambda[fiNumIndexQA];	//! TPC dEdx vs p of Lambda daughters
      TH2D*           fhQAV0sDaughterNumSigmaPionLambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of pion (Lambda candidates)
      TH2D*           fhQAV0sDaughterNumSigmaProtonLambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of proton (Lambda candidates)
      TH2D*           fhQAV0sDaughterNumSigmaPionALambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of pion (Anti-Lambda candidates)
      TH2D*           fhQAV0sDaughterNumSigmaProtonALambda[fiNumIndexQA];  //! number of TPC sigmas vs mother pT of proton (Anti-Lambda candidates)
      TH1D*					  fhQAV0sCPAK0s[fiNumIndexQA];	//! cosine of pointing angle of K0s candidates
      TH1D*					  fhQAV0sCPALambda[fiNumIndexQA];	//! cosine of pointing angle of Lambda candidates
      TH1D*					  fhQAV0sNumTauK0s[fiNumIndexQA];	//! number of c*tau of K0s candidates
      TH1D*					  fhQAV0sNumTauLambda[fiNumIndexQA];	//! number of c*tau of Lambda candidates
      TH2D*				   	fhQAV0sArmenterosK0s[fiNumIndexQA];	//! Armenteros-Podolanski plot for K0s candidates
      TH2D*			  		fhQAV0sArmenterosLambda[fiNumIndexQA];	//! Armenteros-Podolanski plot for Lambda candidates
      TH2D*			  		fhQAV0sArmenterosALambda[fiNumIndexQA];	//! Armenteros-Podolanski plot for ALambda candidates

      AliAnalysisTaskUniFlow(const AliAnalysisTaskUniFlow&); // not implemented
      AliAnalysisTaskUniFlow& operator=(const AliAnalysisTaskUniFlow&); // not implemented

      ClassDef(AliAnalysisTaskUniFlow, 7);
};

#endif
