/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUniFlow_H
#define AliAnalysisTaskUniFlow_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"


class AliAnalysisTaskUniFlow : public AliAnalysisTaskSE
{
    public:
      enum    RunMode {kTest, kFillWeights, kFull}; // task running mode (NOT GRID MODE)
      enum    ColSystem {kPP, kPPb, kPbPb}; // tag for collisional system
      enum    AnalType {kAOD, kESD}; // tag for analysis type
      enum    DataPeriod {kNon, k16k, k16l, k16q, k16r, k16s, k16t}; // tag for data period
      enum    PartSpecies {kUnknown, kCharged, kPion, kKaon, kProton, kK0s, kLambda, kPhi}; // list of all particle species of interest

      AliEventCuts fEventCuts; //

      struct FlowPart // representation of selected particle (species independent) storing only basic properties for flow calculations
      {
                FlowPart(Double_t dPt = 0, Double_t dPhi = 0, Double_t dEta = 0, Short_t iCharge = 0, PartSpecies sSpecies = kUnknown, Double_t dMass = 0, Double_t dPx = 0, Double_t dPy = 0, Double_t dPz = 0) :
                  pt(dPt), px(dPx), py(dPy), pz(dPz), phi(dPhi), eta(dEta), mass(dMass), charge(iCharge), species(sSpecies) {} // constructor

        void    PrintPart() const { printf("pt %g | px %g | py %g| pz %g | phi %g | eta %g | mass %g | charge %d | species %d \n",pt,px,py,pz,phi,eta,mass,charge,species); } // print struct members
        static FlowPart MakeMother(const FlowPart* part1 = 0x0, const FlowPart* part2 = 0x0, const PartSpecies = kUnknown); // reconstruct & fill properties of potential mother from two FlowPart daughters
        static Double_t InvMass(const FlowPart* part1, const FlowPart* part2); // calculate invariant mass from two FlowPart daughters

        Double_t pt,px,py,pz,phi,eta,mass;
        Short_t charge;
        PartSpecies species;
      };
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
      void                    SetSampling(Bool_t sample = kTRUE) { fSampling = sample; }
      void                    SetFillQAhistos(Bool_t fill = kTRUE) { fFillQA = fill; }
      //void                    SetNumberOfSamples(Short_t numSamples = 10) { fNumSamples = numSamples; } // not implemented yet
      void                    SetProcessCharged(Bool_t filter = kTRUE) { fProcessCharged = filter; }
      void                    SetProcessPID(Bool_t filter = kTRUE) { fProcessPID = filter; }
      void                    SetProcessV0s(Bool_t filter = kTRUE) { fProcessV0s = filter; }
      void                    SetProcessPhi(Bool_t filter = kTRUE) { fProcessPhi = filter; }
      // flow related setters
      void                    SetUseFixedMultBins(Bool_t fixed = kTRUE) { fUseFixedMultBins = fixed; }
      void                    SetFlowRFPsPtMin(Float_t pt) { fCutFlowRFPsPtMin = pt; }
      void                    SetFlowRFPsPtMax(Float_t pt) { fCutFlowRFPsPtMax = pt; }
      void                    SetFlowDoFourCorrelations(Bool_t four = kTRUE) { fCutFlowDoFourCorrelations = four; }
      void                    SetFlowFillWeights(Bool_t weights = kTRUE) { fFlowFillWeights = weights; }
      void                    SetUseWeigthsFile(const char* file) { fFlowWeightsPath = file; fFlowUseWeights = kTRUE; }
      // events setters
      void                    SetColisionSystem(ColSystem colSystem = kPP) { fColSystem = colSystem; }
      void                    SetPeriod(DataPeriod period = kNon) { fPeriod = period; }
      void                    SetMultEstimator(const char* mult = "CHARGED") { fMultEstimator = mult; }
      void                    SetTrigger(Short_t trigger = 0) { fTrigger = trigger; }
      void                    SetUseAliEventCuts(Bool_t bUseCuts = kTRUE) { fUseAliEventCuts = bUseCuts; }
      void					          SetPVtxZMax(Double_t z) { fPVtxCutZ = z; }
      // track setters
      void                    SetChargedEtaMax(Double_t eta) { fCutChargedEtaMax = eta; }
      void                    SetChargedPtMax(Double_t pt) { fCutChargedPtMax = pt; }
      void                    SetChargedPtMin(Double_t pt) { fCutChargedPtMin = pt; }
      void                    SetChargedDCAzMax(Double_t dcaz) {  fCutChargedDCAzMax = dcaz; }
      void                    SetChargedDCAxyMax(Double_t dcaxy) {  fCutChargedDCAxyMax = dcaxy; }
      void                    SetChargedNumTPCclsMin(UShort_t tpcCls) { fCutChargedNumTPCclsMin = tpcCls; }
      void                    SetChargedTrackFilterBit(UInt_t filter) { fCutChargedTrackFilterBit = filter; }
      // PID (pi,K,p) setters
      void                    SetPIDUseAntiProtonOnly(Bool_t use = kTRUE) { fCutPIDUseAntiProtonOnly = use; }
      void                    SetPIDNumSigmasPionMax(Double_t numSigmas) { fCutPIDnSigmaPionMax = numSigmas; }
      void                    SetPIDNumSigmasKaonMax(Double_t numSigmas) { fCutPIDnSigmaKaonMax = numSigmas; }
      void                    SetPIDNumSigmasProtonMax(Double_t numSigmas) { fCutPIDnSigmaProtonMax = numSigmas; }
      void                    SetPIDNumSigmasCombinedNoTOFrejection(Bool_t reject = kTRUE) { fCutPIDnSigmaCombinedNoTOFrejection = reject; }
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
      void                    SetV0sDaughterNumTPCrossMin(Float_t cls) { fCutV0sDaughterNumTPCCrossMin = cls; }
      void                    SetV0sDaughterNumTPCFindMin(UShort_t cls) { fCutV0sDaughterNumTPCFindMin = cls; }
      void                    SetV0sDaughterNumTPCClsPIDMin(UShort_t cls) { fCutV0sDaughterNumTPCClsPIDMin = cls; }
      void                    SetV0sDaughterRatioCrossFindMin(Float_t ratio) { fCutV0sDaughterRatioCrossFindMin = ratio; }
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
      void                    SetV0sK0sPionNumTPCSigmaMax(Double_t nSigma) { fCutV0sK0sPionNumTPCSigmaMax = nSigma; }
      void                    SetV0sLambdaPionNumTPCSigmaMax(Double_t nSigma) { fCutV0sLambdaPionNumTPCSigmaMax = nSigma; }
      void                    SetV0sLambdaProtonNumTPCSigmaMax(Double_t nSigma) { fCutV0sLambdaProtonNumTPCSigmaMax = nSigma; }
      // phi setters
      void					          SetPhiMotherEtaMax(Double_t eta) { fCutPhiMotherEtaMax = eta; }
      void					          SetPhiInvMassMin(Double_t mass) { fCutPhiInvMassMin = mass; }
      void					          SetPhiInvMassMax(Double_t mass) { fCutPhiInvMassMax = mass; }


      // related to Alex's code
      Bool_t                  IsV0SelectedK0sAlex(const AliAODv0* v0 = 0x0); //
      Double_t                GetRapidity(Double_t mass, Double_t Pt, Double_t Eta); //
      Bool_t       fDoAlexK0sSelection; //

      // end of Alex's code

    private:
      // array lenghts & constants
      const Double_t          fPDGMassPion; // [DPGMass] DPG mass of charged pion
      const Double_t          fPDGMassKaon; // [DPGMass] DPG mass of charged kaon
      const Double_t          fPDGMassProton; // [DPGMass] DPG mass of proton
      const Double_t          fPDGMassPhi; // [DPGMass] DPG mass of phi (333) meson
      const Double_t          fPDGMassK0s; // [DPGMass] DPG mass of K0s
      const Double_t          fPDGMassLambda; // [DPGMass] DPG mass of (Anti)Lambda
      static const Short_t    fFlowNumHarmonicsMax = 10; // maximum harmonics length of flow vector array
      static const Short_t    fFlowNumWeightPowersMax = 10; // maximum weight power length of flow vector array
      static const Short_t    fFlowPOIsPtNumBins = 200; // number of pT bins for POIs
      const Float_t           fFlowPOIsPtMin; // [0] (GeV/c) min pT treshold for POIs for differential flow
      const Float_t           fFlowPOIsPtMax; // [20] (GeV/c) max pT treshold for POIs for differential flow
      const Int_t             fFlowCentMin; // [0] min range for centrality/multiplicity histos
      const Int_t             fFlowCentMax; // [150] min range for centrality/multiplicity histos
      const Int_t             fFlowCentNumBins; // [150] min range for centrality/multiplicity histos
      static const Short_t    fV0sNumBinsMass = 60; // number of InvMass bins for V0s distribution
      static const Short_t    fPhiNumBinsMass = 60; // number of InvMass bins for phi distribution
      static const Short_t    fiNumIndexQA = 2; // QA indexes: 0: before cuts // 1: after cuts

      const static Short_t    fNumSamples = 10; // overall number of samples (from random sampling) used
      const static Int_t      fNumHarmonics = 1; // number of harmonics
      static Int_t            fHarmonics[fNumHarmonics]; // values of used harmonics
      const static Int_t      fNumEtaGap = 1; // number of harmonics
      static Double_t         fEtaGap[fNumEtaGap]; // values of used harmonics
      const static Int_t      fNumMultBins = 6; // number of multiplicity bins
      static Double_t         fMultBins[fNumMultBins+1]; // multiplicity bins

      Bool_t                  InitializeTask(); // called once on beginning of task (within CreateUserObjects method)
      void                    ListParameters(); // list all task parameters

      Bool_t                  EventSelection(); // main method for event selection (specific event selection is applied within)
      Bool_t                  IsEventSelected_2016(); // event selection for LHC2016 pp & pPb data
      void                    FillEventsQA(const Short_t iQAindex); // filling QA plots related to event selection
      Short_t                 GetSamplingIndex(); // returns sampling index based on sampling selection (number of samples)
      Short_t                 GetCentralityIndex(); // returns centrality index based centrality estimator or number of selected tracks
      Bool_t                  HasTrackPIDTPC(const AliAODTrack* track); // is TPC PID OK for this track ?
      Bool_t                  HasTrackPIDTOF(const AliAODTrack* track); // is TOF PID OK for this track ?

      Bool_t                  ProcessEvent(); // main (envelope) method for processing events passing selection

      void                    Filtering(); // main (envelope) method for filtering all POIs in event
      void                    FilterCharged(); // charged tracks filtering
      Bool_t                  IsChargedSelected(const AliAODTrack* track = 0x0); // charged track selection
      void                    FillQARefs(const Short_t iQAindex, const AliAODTrack* track = 0x0); // filling QA plots for RFPs selection
      void                    FillQACharged(const Short_t iQAindex, const AliAODTrack* track = 0x0); // filling QA plots for charged track selection
      void                    FilterPID(); // pi,K,p filtering
      PartSpecies             IsPIDSelected(const AliAODTrack* track); // PID tracks selections
      void                    FillPIDQA(const Short_t iQAindex, const AliAODTrack* track = 0x0, const PartSpecies species = kUnknown); // filling pi,K,p QA histograms
      void                    FilterV0s(); // K0s, Lambda, ALambda filtering
      Bool_t                  IsV0Selected(const AliAODv0* v0 = 0x0); // general (common) V0 selection
      Bool_t                  IsV0aK0s(const AliAODv0* v0 = 0x0); // V0 selection: K0s specific
      Short_t                 IsV0aLambda(const AliAODv0* v0 = 0x0); // V0 selection: (A)Lambda specific
      void                    FillQAV0s(const Short_t iQAindex, const AliAODv0* v0 = 0x0, const Bool_t bIsK0s = kTRUE, const Short_t bIsLambda = 2); // filling QA plots for V0s candidates
      void                    FilterPhi(); // reconstruction and filtering of Phi meson candidates
      void                    FillQAPhi(const Short_t iQAindex, const FlowPart* part = 0x0); // filling QA plots for V0s candidates
      // Flow related methods
      void                    DoFlowRefs(const Short_t iEtaGapIndex = 0); // Estimate <2> for reference flow
      void                    DoFlowCharged(const Short_t iEtaGapIndex = 0); // Estimate <2'> for pt diff. flow of charged hadrons
      void                    DoFlowPID(const Short_t iEtaGapIndex = 0, const PartSpecies species = kUnknown); // Estimate <2'> for pt diff. flow of PID (pi,K,p) hadrons
      void                    DoFlowPhi(const Short_t iEtaGapIndex = 0, const Short_t iMassIndex = 0); // Estimate <2'> for pt diff. flow of phi particles
      void                    DoFlowV0s(const Short_t iEtaGapIndex = 0, const Short_t iMassIndex = 0, const PartSpecies species = kUnknown); // Estimate <2'> for pt diff. flow of V0 particles
      void                    FillRefsVectors(const Short_t iEtaGapIndex = 0); // fill flow vector Q with RFPs for reference flow
      void                    FillPOIsVectors(const Short_t iEtaGapIndex = 0, const PartSpecies species = kUnknown, const Short_t iMassIndex = 0); // fill flow vectors p,q and s with POIs (for given species) for differential flow calculations
      Short_t                 GetPOIsPtBinIndex(const Double_t pt); // return pT bin index based on momenta value
      void                    ResetRFPsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]); // set values to TComplex(0,0,0) for given array
      void                    ResetPOIsVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]); // set values to TComplex(0,0,0) for given array
      void                    ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]); // printf all values of given Flow vector array

      TComplex                Q(const Short_t n, const Short_t p);
      TComplex                QGapPos(const Short_t n, const Short_t p);
      TComplex                QGapNeg(const Short_t n, const Short_t p);
      TComplex                P(const Short_t n, const Short_t p, const Short_t pt);
      TComplex                PGapPos(const Short_t n, const Short_t p, const Short_t pt);
      TComplex                PGapNeg(const Short_t n, const Short_t p, const Short_t pt);
      TComplex                S(const Short_t n, const Short_t p, const Short_t pt);

      TComplex                Two(const Short_t n1, const Short_t n2); // Two particle reference correlation calculations (no eta gap)
      TComplex                TwoGap(const Short_t n1, const Short_t n2); // Two particle reference correlation calculations (with eta gap)
      TComplex                TwoDiff(const Short_t n1, const Short_t n2, const Short_t pt); // Two particle diff. correlation calculations (no eta gap)
      TComplex                TwoDiffGapPos(const Short_t n1, const Short_t n2, const Short_t pt); // Two particle diff. correlation calculations (with eta gap)
      TComplex                TwoDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t pt); // Two particle diff. correlation calculations (with eta gap)
      TComplex                Four(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4); // Four particle reference correlation calculations (no eta gap)
      TComplex                FourDiff(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4, const Short_t pt); // Four particle reference correlation calculations (no eta gap)

      // properties
      AliAODEvent*            fEventAOD; //! AOD event countainer
      AliPIDResponse*         fPIDResponse; //! AliPIDResponse container
      AliPIDCombined*         fPIDCombined; //! AliPIDCombined container
      TFile*                  fFlowWeightsFile; //! source file containing weights
      Bool_t                  fInit; // initialization check
      Short_t                 fIndexSampling; // sampling index (randomly generated)
      Short_t                 fIndexCentrality; // centrality bin index (based on centrality est. or number of selected tracks)
      Short_t                 fEventCounter; // event counter (used for local test runmode purpose)
      Short_t                 fNumEventsAnalyse; // [50] number of events to be analysed / after passing selection (only in test mode)
      Int_t                   fRunNumber; // [-1] run number obtained from AliVHeader

      TComplex                fFlowVecQpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecQneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
      TComplex                fFlowVecPpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]; // flow vector array for flow calculation
      TComplex                fFlowVecPneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]; // flow vector array for flow calculation
      TComplex                fFlowVecS[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax][fFlowPOIsPtNumBins]; // flow vector array for flow calculation

      // selected POIs containers
      std::vector<FlowPart>*  fVectorCharged; //! container for selected charged particles
      std::vector<FlowPart>*  fVectorPion; //! container for selected pion candidates
      std::vector<FlowPart>*  fVectorKaon; //! container for selected kaon candidates
      std::vector<FlowPart>*  fVectorProton; //! container for selected proton candidates
      std::vector<FlowPart>*  fVectorK0s; //! container for selected K0s candidates
      std::vector<FlowPart>*  fVectorLambda; //! container for selected (Anti)Lambda candidates
      std::vector<FlowPart>*  fVectorPhi; //! container for selected phi candidates (unlike-sign pairs)

      //cuts & selection: analysis
      RunMode                 fRunMode; // running mode (not grid related)
      AnalType                fAnalType; // analysis type: AOD / ESD
      Bool_t                  fSampling;      // Do random sampling ? (estimation of vn stat. uncertanity)
      Bool_t                  fFillQA; //[kTRUE] flag for filling the QA plots
      Bool_t                  fProcessCharged; // flag for processing charged tracks (both RPF and POIs)
      Bool_t                  fProcessPID; // flag for processing PID tracks (pi,K,p)
      Bool_t                  fProcessV0s; // flag for processing V0 candidates (K0s, Lambda/ALambda)
      Bool_t                  fProcessPhi; // flag for processing Phi meson candidates
      // cuts & selection: flow related
      Bool_t                  fUseFixedMultBins; // [kFALSE] setting fixed multiplicity bins
      Float_t                 fCutFlowRFPsPtMin; // [0] (GeV/c) min pT treshold for RFPs particle for reference flow
      Float_t                 fCutFlowRFPsPtMax; // [0] (GeV/c) max pT treshold for RFPs particle for reference flow
      Bool_t                  fCutFlowDoFourCorrelations; // [kTRUE] flag for processing <4>
      Bool_t                  fFlowFillWeights; //[kFALSE] flag for filling weights
      Bool_t                  fFlowUseWeights; //[kFALSE] flag for using the previously filled weights (NOTE: this is turned on only when path to file is applied via fFlowWeightsPath)
      TString                 fFlowWeightsPath; //[] path to source root file with weigthts (if empty unit weights are applied) e.g. "alice/cern.ch/user/k/kgajdoso/EfficienciesWeights/2016/PhiWeight_LHC16kl.root"

      //cuts & selection: events
      ColSystem               fColSystem; // collisional system
      DataPeriod              fPeriod; // period of analysed data sample (e.g. LHC16k, ...)
      Short_t                 fTrigger; // physics selection trigger
      TString                 fMultEstimator; // [''] multiplicity estimator (suported: ''/Charged,VOA,V0C,V0M,CL0,CL1,ZNA,ZNC)
      Bool_t                  fUseAliEventCuts; // use decision of AliEventCuts in event selection
      Float_t                 fPVtxCutZ; // (cm) PV z cut
      //cuts & selection: tracks
      UInt_t                  fCutChargedTrackFilterBit; // (-) tracks filter bit
      UShort_t                fCutChargedNumTPCclsMin;  // (-) Minimal number of TPC clusters used for track reconstruction
      Float_t                 fCutChargedEtaMax; // (-) Maximum pseudorapidity range
      Float_t                 fCutChargedPtMax; // (GeV/c) Maximal track pT
      Float_t                 fCutChargedPtMin; // (GeV/c) Minimal track pT
      Float_t                 fCutChargedDCAzMax; // (cm) Maximal DCA-z cuts for tracks (pile-up rejection suggested for LHC16)
      Float_t                 fCutChargedDCAxyMax; // (cm) Maximal DCA-xy cuts for tracks (pile-up rejection suggested for LHC16)
      // cuts & selection: PID selection
      Bool_t                  fCutPIDUseAntiProtonOnly; // [kFALSE] check proton PID charge to select AntiProtons only
      Double_t                fCutPIDnSigmaPionMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for pion candidates
      Double_t                fCutPIDnSigmaKaonMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for kaon candidates
      Double_t                fCutPIDnSigmaProtonMax; // [3] maximum of nSigmas (TPC or TPC & TOF combined) for proton candidates
      Double_t                fCutPIDnSigmaTPCRejectElectron; // [3] number of TPC nSigma for electron rejection
      Bool_t                  fCutPIDnSigmaCombinedNoTOFrejection; // [kFALSE] flag for rejection candidates in TPC+TOF pt region if TOF is not available (if true and no TOF, only TPC is used)
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
      Float_t                 fCutV0sDaughterNumTPCCrossMin; // min number of crossed TPC rows
      UShort_t                fCutV0sDaughterNumTPCFindMin; // min number of findable TPC clusters
      UShort_t                fCutV0sDaughterNumTPCClsPIDMin; // min number of TPC clusters used for PID
      Float_t                 fCutV0sDaughterRatioCrossFindMin; // min ratio of crossed / findable TPC clusters
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
      Double_t                fCutV0sK0sPionNumTPCSigmaMax; // (sigmaTPC) max number of TPC sigmas for kaon PID (K0s candidates)
      Double_t                fCutV0sLambdaPionNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for pion PID (Lambda candidates)
      Double_t                fCutV0sLambdaProtonNumTPCSigmaMax;    // (sigmaTPC) max number of TPC sigma for proton PID (Lambda candidates)
      // cuts & selection: phi
      Double_t                fCutPhiMotherEtaMax; // (-) max value of phi candidate pseudorapidity
      Double_t                fCutPhiInvMassMin; // [0.99] (GeV/c2) min inv. mass window for selected phi candidates
      Double_t                fCutPhiInvMassMax; // [1.07] (GeV/c2) min inv. mass window for selected phi candidates

      // output lists
      TList*      fQAEvents; //! events list
      TList*      fQACharged; //! charged tracks list
      TList*      fQAPID; //! pi,K,p list
      TList*      fQAV0s; //! V0s candidates list
      TList*      fQAPhi; //! Phi candidates list
      TList*      fFlowWeights; //! list for flow weights
      TList*      fFlowRefs; //! list for flow of reference particles
      TList*      fFlowCharged; //! list for flow of charged particles
      TList*      fFlowPID; //! list for flow of PID (pi,K,p) particles
      TList*      fFlowPhi; //! list for flow of Phi particles
      TList*      fFlowK0s; //! list for flow of K0s particles
      TList*      fFlowLambda; //! list for flow of (anti-)Lambda particles

      // histograms & profiles

      // Flow
      TH3D*           fh3V0sEntriesK0sPos[fNumEtaGap]; //! distribution of K0s candidates (cent, pT, InvMass)
      TH3D*           fh3V0sEntriesK0sNeg[fNumEtaGap]; //! distribution of K0s candidates (cent, pT, InvMass)
      TH3D*           fh3V0sEntriesLambdaPos[fNumEtaGap]; //! distribution of (Anti-)Lambda candidates (cent, pT, InvMass)
      TH3D*           fh3V0sEntriesLambdaNeg[fNumEtaGap]; //! distribution of (Anti-)Lambda candidates (cent, pT, InvMass)
      TH3D*           fh3PhiEntriesSignalPos[fNumEtaGap]; //! distribution of phi candidates / unlike-sign pairs (cent, pT, InvMass)
      TH3D*           fh3PhiEntriesSignalNeg[fNumEtaGap]; //! distribution of phi candidates / unlike-sign pairs (cent, pT, InvMass)
      TH3D*           fh3PhiEntriesBGPos[fNumEtaGap]; //! distribution of phi background candidates / like-sign pairs (cent, pT, InvMass)
      TH3D*           fh3PhiEntriesBGNeg[fNumEtaGap]; //! distribution of phi background candidates / like-sign pairs (cent, pT, InvMass)

      TH3D*           fh3WeightsRefs; //! distribution of Refs particles for estimating weight purpose (phi,eta,pt)
      TH3D*           fh3WeightsCharged; //! distribution of Charged POIs particles for estimating weight purpose (phi,eta,pt)
      TH3D*           fh3WeightsPion; //! distribution of Pion POIs particles for estimating weight purpose (phi,eta,pt)
      TH3D*           fh3WeightsKaon; //! distribution of Kaon POIs particles for estimating weight purpose (phi,eta,pt)
      TH3D*           fh3WeightsProton; //! distribution of Proton POIs particles for estimating weight purpose (phi,eta,pt)
      TH3D*           fh3WeightsPhi; //! distribution of Phi POIs particles for estimating weight purpose (phi,eta,pt)
      TH3D*           fh3WeightsK0s; //! distribution of K0s POIs particles for estimating weight purpose (phi,eta,pt)
      TH3D*           fh3WeightsLambda; //! distribution of Lambda POIs particles for estimating weight purpose (phi,eta,pt)

      TH3D*           fh3AfterWeightsRefs; //! distribution of Refs particles after applying the weights (phi,eta,pt)
      TH3D*           fh3AfterWeightsCharged; //! distribution of Charged POIs particles after applying the weights (phi,eta,pt)
      TH3D*           fh3AfterWeightsPion; //! distribution of Pion POIs particles after applying the weights (phi,eta,pt)
      TH3D*           fh3AfterWeightsKaon; //! distribution of Kaon POIs particles after applying the weights (phi,eta,pt)
      TH3D*           fh3AfterWeightsProton; //! distribution of Proton POIs particles after applying the weights (phi,eta,pt)
      TH3D*           fh3AfterWeightsPhi; //! distribution of Phi POIs particles after applying the weights (phi,eta,pt)
      TH3D*           fh3AfterWeightsK0s; //! distribution of K0s POIs particles after applying the weights (phi,eta,pt)
      TH3D*           fh3AfterWeightsLambda; //! distribution of Lambda POIs particles after applying the weights (phi,eta,pt)

      TH2D*           fh2WeightRefs; //! container for loading weights for given run
      TH2D*           fh2WeightCharged; //! container for loading weights for given run
      TH2D*           fh2WeightPion; //! container for loading weights for given run
      TH2D*           fh2WeightKaon; //! container for loading weights for given run
      TH2D*           fh2WeightProton; //! container for loading weights for given run
      TH2D*           fh2WeightK0s; //! container for loading weights for given run
      TH2D*           fh2WeightLambda; //! container for loading weights for given run
      TH2D*           fh2WeightPhi; //! container for loading weights for given run

      TProfile*       fpMeanQxRefsPos[fNumEtaGap][fNumHarmonics]; //! average of Qx (vs. centrality) for Refs
      TProfile*       fpMeanQxRefsNeg[fNumEtaGap][fNumHarmonics]; //! average of Qx (vs. centrality) for Refs
      TProfile*       fpMeanQyRefsPos[fNumEtaGap][fNumHarmonics]; //! average of Qy (vs. centrality) for Refs
      TProfile*       fpMeanQyRefsNeg[fNumEtaGap][fNumHarmonics]; //! average of Qy (vs. centrality) for Refs

      TProfile*       fpRefsCor2[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2> correlations for RFPs
      TProfile2D*     fp2ChargedCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for Charged tracks POIs: POIs in Eta>0
      TProfile2D*     fp2ChargedCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for Charged tracks POIs: POIs in Eta<0
      TProfile2D*     fp2PionCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for pion POIs: POIs in Eta>0
      TProfile2D*     fp2PionCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for pion POIs: POIs in Eta>0
      TProfile2D*     fp2KaonCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for kaon POIs: POIs in Eta>0
      TProfile2D*     fp2KaonCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for kaon POIs: POIs in Eta>0
      TProfile2D*     fp2ProtonCor2Pos[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for proton POIs: POIs in Eta>0
      TProfile2D*     fp2ProtonCor2Neg[fNumSamples][fNumEtaGap][fNumHarmonics]; //! <2'> correlations for proton POIs: POIs in Eta>0
      TProfile3D*     fp3V0sCorrK0sCor2Pos[fNumEtaGap][fNumHarmonics]; //! <2'> correlations of K0s candidates: POIs in Eta>0  (cent, pT, InvMass)
      TProfile3D*     fp3V0sCorrK0sCor2Neg[fNumEtaGap][fNumHarmonics]; //! <2'> correlations of K0s candidates: POIs in Eta<0  (cent, pT, InvMass)
      TProfile3D*     fp3V0sCorrLambdaCor2Pos[fNumEtaGap][fNumHarmonics]; //! <2'> correlations of (Anti-)Lambda candidates: POIs in Eta>0  (cent, pT, InvMass)
      TProfile3D*     fp3V0sCorrLambdaCor2Neg[fNumEtaGap][fNumHarmonics]; //! <2'> correlations of (Anti-)Lambda candidates: POIs in Eta<0  (cent, pT, InvMass)
      TProfile3D*     fp3PhiCorrCor2Pos[fNumEtaGap][fNumHarmonics]; //! <2'> correlations of phi candidates / unlike-sign pairs: POIs in Eta>0 (cent, pT, InvMass)
      TProfile3D*     fp3PhiCorrCor2Neg[fNumEtaGap][fNumHarmonics]; //! <2'> correlations of phi candidates / unlike-sign pairs: POIs in Eta<0 (cent, pT, InvMass)

      TProfile*       fpRefsCor4[fNumSamples][fNumHarmonics]; //! <4> correlations for RFPs
      TProfile2D*     fp2ChargedCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for Charged tracks POIs
      TProfile2D*     fp2PionCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for pion POIs
      TProfile2D*     fp2KaonCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for kaon POIs
      TProfile2D*     fp2ProtonCor4[fNumSamples][fNumHarmonics]; //! <4'> correlations for proton POIs
      TProfile3D*     fp3V0sCorrK0sCor4[fNumHarmonics]; //! <4'> correlations of K0s candidates (cent, pT, InvMass)
      TProfile3D*     fp3V0sCorrLambdaCor4[fNumHarmonics]; //! <4'> correlations of (Anti-)Lambda candidates (cent, pT, InvMass)
      TProfile3D*     fp3PhiCorrCor4[fNumHarmonics]; //! <4'> correlations of phi candidates / unlike-sign pairs (cent, pT, InvMass)

      // Events
      TH2D*           fhEventSampling; //! distribution of sampled events (based on randomly generated numbers)
      TH1D*           fhEventCentrality; //! distribution of event centrality
      TH2D*           fh2EventCentralityNumSelCharged; //! distribution of event centrality vs number of selected charged tracks
      TH1D*           fhEventCounter; //! counter following event selection
      // Charged
      TH1D*           fhRefsMult; //!multiplicity distribution of selected RFPs
      TH1D*           fhRefsPt; //! pt distribution of selected RFPs
      TH1D*           fhRefsEta; //! pt distribution of selected RFPs
      TH1D*           fhRefsPhi; //! pt distribution of selected RFPs
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


      // Alex V0s selection
      Float_t      fNcrFind;            // number of cross rows over findable clusters
      Float_t      fDCADghtPV;          // DCA daughters to primary vertex
      Float_t      fMaxDCADght;         // DCA daughters
      Float_t      fCosPA;              // cos pointing angle
      Float_t      fMinRad;             // V0 radius cut low
      Float_t      fMaxRad;             // V0 radius cut low
      Bool_t       fArmPodCut;          // flag for Arm-Pod cut
      Bool_t       fMinPtDght;          // flag for min pT cut for daughters
      Double_t     fEtaCut;             // Eta cut used to select particles
      Int_t        fNoClus;	          // No of TPC clusters
      Double_t     fMinPt;              // Min pt - for histogram limits
      Double_t     fMaxPt;              // Max pt - for histogram limits
      Double_t     fNsigCut;            // combined sigma cut value
      Int_t        fNoClusPid;          // no of TPC clusters used for PID
      // end of Alex selection




      AliAnalysisTaskUniFlow(const AliAnalysisTaskUniFlow&); // not implemented
      AliAnalysisTaskUniFlow& operator=(const AliAnalysisTaskUniFlow&); // not implemented

      ClassDef(AliAnalysisTaskUniFlow, 4);
};

#endif
