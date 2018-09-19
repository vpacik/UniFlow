/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// =================================================================================================
// AliAnalysisTaskUniFlow - ALICE Unified Flow framework
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016-2018
// =================================================================================================
//
// ALICE analysis task for universal study of flow via 2-(or multi-)particle correlations
// using Generic Framework notation for calculations including per-particle weights.
//
// Implemented flow calculation for both reference & pt-differential flow
// of both inclusive charged (done anyway) & identified particles (pi,K,p,K0s,Lambda,phi).
//
// Note: So far implemented only for AOD analysis and tuned on Run2 pp & pPb analyses!
//
// PLEASE READ THE INSTRUCTION BELLOW BEFORE RUNNING !!!
//
// =================================================================================================
// Analysis can run in these modes setup via AliAnalysisTaskUniFlow::SetRunMode(RunMode)
//  -- AnalysisTaskUniFlow::kFull : running mode
//      - full scale analysis
//
//  -- AnalysisTaskUniFlow::kTest : development / testing / debugging mode
//      - only limited number of events is processed (AliAnalysisTaskUniFlow::SetNumEventsAnalyse(Int_t))
//
//  -- AnalysisTaskUniFlow::kSkipFlow : usable for QA and(or) weights estimation before full scale running
//      - events are processed whilst skipping correlation calculations, i.e. event loop ends after particles are filtered
//
//  NOTE: by default, all modes includes :
//  - filling QA plots : can be turned-off by AliAnalysisTaskUniFlow::SetFillQAhistos(kFALSE)
//    (might be heavy while running several tasks at once)
//  - filling GF weights : can be turned-off by AliAnalysisTaskUniFlow::SetFlowFillWeights(kFALSE)
//
// =================================================================================================
// Overview of analysis flow (see implementation of corresonding method for details)
// 1) Event selection : EventSelection()
//
// 2) Particle selection & reconstruction : Filtering()
//        - whether or not are particles processed is driven by 'fProcess?' flag setup by AliAnalysisTaskUniFlow::SetProcess?(kTRUE)
//          (except for incl. charged, which are processed anyway)
//        - filling QA plots
//        - filling GF weights
//
//    !!! here the event loop ends in case of running in 'kSkipFlow' mode
//
// 3) Flow / correlation calculations : DoFlowPOIs(), DoFlowRefs()
//        - to setup harmonics to be calculated, modify 'AliAnalysisTaskUniFlow::fHarmonics[]' (in .cxx) AND 'fNumHarmonics' (in .h)
//        - to setup eta gaps (2part) to be calculated, modify 'AliAnalysisTaskUniFlow::fEtaGap[]' (in .cxx) AND 'fNumEtaGap' (in .h)
//              - for NO gap case, fill in value '-1.0'
//        - 2-particle cumulants are run by ©default (in 'kFull' mode)
//        - 4-particle cumulants has to be setup by invoking AliAnalysisTaskUniFlow::SetFlowDoFourCorrelations(kTRUE)
//
// =================================================================================================

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TComplex.h"
#include "TRandom3.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliPicoTrack.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

#include "AliAnalysisTaskUniFlow.h"

class AliAnalysisTaskUniFlow;

ClassImp(AliAnalysisTaskUniFlow);

Double_t AliAnalysisTaskUniFlow::fMultBins[] = {0.,5.,10.,20.,40.,60.,100.};

AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow() : AliAnalysisTaskSE(),
  fEventAOD(0x0),
  fPVz(0.0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fFlowWeightsFile(0x0),
  fInit(kFALSE),
  fMC(kFALSE),
  fArrayMC(0x0),
  fIndexSampling(0),
  fIndexCentrality(-1),
  fEventCounter(0),
  fNumEventsAnalyse(50),
  fRunNumber(-1),
  fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
  fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
  fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),
  fPDGMassPhi(TDatabasePDG::Instance()->GetParticle(333)->Mass()),
  fPDGMassK0s(TDatabasePDG::Instance()->GetParticle(310)->Mass()),
  fPDGMassLambda(TDatabasePDG::Instance()->GetParticle(3122)->Mass()),

  // analysis selection
  fRunMode(kFull),
  fAnalType(kAOD),
  fSampling(kTRUE),
  fFillQA(kTRUE),
  // fNumSamples(10),

  // flow related
  fUseFixedMultBins(kFALSE),
  fCutFlowRFPsPtMin(0.2),
  fCutFlowRFPsPtMax(5.0),
  fFlowPOIsPtMin(0.0),
  fFlowPOIsPtMax(15.0),
  fFlowFillWeights(kTRUE),
  fFlowCentMin(0),
  fFlowCentMax(0),
  fFlowWeightsPath(),
  fFlowRunByRunWeights(kTRUE),
  fFlowUseWeights(kFALSE),
  fFlowUse3Dweights(kFALSE),

  // events selection
  fPVtxCutZ(10.0),
  fColSystem(kPPb),
  fMultEstimator(kV0A),
  fTrigger(AliVEvent::kINT7),
  fEventRejectAddPileUp(kTRUE),

  // charged tracks selection
  fCutChargedEtaMax(0.8),
  fCutChargedDCAzMax(0.),
  fCutChargedDCAxyMax(0.),
  fCutChargedTrackFilterBit(96),
  fCutChargedNumTPCclsMin(70),

  // PID tracks selection
  fCutPIDUseAntiProtonOnly(kFALSE),
  fCutPIDnSigmaPionMax(3.),
  fCutPIDnSigmaKaonMax(3.),
  fCutPIDnSigmaProtonMax(3.),
  fCutPIDnSigmaTPCRejectElectron(3.),
  fCutPIDnSigmaCombinedTOFrejection(kTRUE),
  fCutUseBayesPID(kTRUE),
  fCutPIDBayesPionMin(0.8),
  fCutPIDBayesKaonMin(0.8),
  fCutPIDBayesProtonMin(0.8),
  fCutPIDBayesRejectElectron(0.8),
  fCutPIDBayesRejectMuon(0.8),

  // V0s selection
  fCutV0sOnFly(kFALSE),
  fCutV0srejectKinks(kTRUE),
  fCutV0sDaughterNumTPCClsMin(0),
  fCutV0sDaughterNumTPCCrossMin(70),
  fCutV0sDaughterNumTPCFindMin(1),
  fCutV0sDaughterNumTPCClsPIDMin(0),
  fCutV0sDaughterRatioCrossFindMin(0.8),
  fCutV0srefitTPC(kTRUE),
  fCutV0sCrossMassRejection(kTRUE),
  fCutV0sCrossMassCutK0s(0.005),
  fCutV0sCrossMassCutLambda(0.010),
  fCutV0sCPAK0sMin(0.97),
  fCutV0sCPALambdaMin(0.995),
  fCutV0sDCAtoPVMin(0.06),
  fCutV0sDCAtoPVMax(0.),
  fCutV0sDCAtoPVzMax(0.),
  fCutV0sDCADaughtersMin(0.),
  fCutV0sDCADaughtersMax(1.0),
  fCutV0sDecayRadiusMin(0.5),
  fCutV0sDecayRadiusMax(200.0),
  fCutV0sDaughterFilterBit(0),
  fCutV0sDaughterPtMin(0.),
  fCutV0sDaughterPtMax(0.),
  fCutV0sDaughterEtaMax(0.8),
  fCutV0sMotherEtaMax(0.8),
  fCutV0sMotherRapMax(0.),
  fCutV0sArmenterosAlphaK0sMin(0.2),
  fCutV0sArmenterosAlphaLambdaMax(0.),
  fCutV0sInvMassK0sMin(0.4),
  fCutV0sInvMassK0sMax(0.6),
  fCutV0sInvMassLambdaMin(1.08),
  fCutV0sInvMassLambdaMax(1.16),
  fCutV0sNumTauK0sMax(7.46),
  fCutV0sNumTauLambdaMax(3.8),
  fCutV0sK0sPionNumTPCSigmaMax(5.0),
  fCutV0sLambdaPionNumTPCSigmaMax(5.0),
  fCutV0sLambdaProtonNumTPCSigmaMax(5.0),

  // phi selection
  fCutPhiMotherEtaMax(0.8),
  fCutPhiInvMassMin(0.99),
  fCutPhiInvMassMax(1.07),

  // output lists
  fQAEvents(0x0),
  fQACharged(0x0),
  fQAPID(0x0),
  fQAV0s(0x0),
  fQAPhi(0x0),
  fFlowWeights(0x0),

  // flow histograms & profiles
  fhsV0sCandK0s(0x0),
  fhsV0sCandLambda(0x0),
  fhsPhiCandSig(0x0),
  fhsPhiCandBg(0x0),

  // event histograms
  fhEventSampling(0x0),
  fhEventCentrality(0x0),
  fh2EventCentralityNumRefs(0x0),
  fhEventCounter(0x0),
  fhQAEventsfMult32vsCentr(0x0),
  fhQAEventsMult128vsCentr(0x0),
  fhQAEventsfMultTPCvsTOF(0x0),
  fhQAEventsfMultTPCvsESD(0x0),

  // charged histogram
  fhRefsMult(0x0),
  fhRefsPt(0x0),
  fhRefsEta(0x0),
  fhRefsPhi(0x0),
  fpRefsMult(0x0),
  fhChargedCounter(0x0),

  // PID histogram
  fhMCRecoSelectedPionPt(0x0),
  fhMCRecoSelectedTruePionPt(0x0),
  fhMCRecoAllPionPt(0x0),
  fhMCGenAllPionPt(0x0),
  fhMCRecoSelectedKaonPt(0x0),
  fhMCRecoSelectedTrueKaonPt(0x0),
  fhMCRecoAllKaonPt(0x0),
  fhMCGenAllKaonPt(0x0),
  fhMCRecoSelectedProtonPt(0x0),
  fhMCRecoSelectedTrueProtonPt(0x0),
  fhMCRecoAllProtonPt(0x0),
  fhMCGenAllProtonPt(0x0),
  fhPIDPionMult(0x0),
  fhPIDPionPt(0x0),
  fhPIDPionPhi(0x0),
  fhPIDPionEta(0x0),
  fhPIDPionCharge(0x0),
  fhPIDKaonMult(0x0),
  fhPIDKaonPt(0x0),
  fhPIDKaonPhi(0x0),
  fhPIDKaonEta(0x0),
  fhPIDKaonCharge(0x0),
  fhPIDProtonMult(0x0),
  fhPIDProtonPt(0x0),
  fhPIDProtonPhi(0x0),
  fhPIDProtonEta(0x0),
  fhPIDProtonCharge(0x0),
  fh2PIDPionTPCdEdx(0x0),
  fh2PIDPionTOFbeta(0x0),
  fh2PIDKaonTPCdEdx(0x0),
  fh2PIDKaonTOFbeta(0x0),
  fh2PIDProtonTPCdEdx(0x0),
  fh2PIDProtonTOFbeta(0x0),
  fh2PIDPionTPCnSigmaPion(0x0),
  fh2PIDPionTOFnSigmaPion(0x0),
  fh2PIDPionTPCnSigmaKaon(0x0),
  fh2PIDPionTOFnSigmaKaon(0x0),
  fh2PIDPionTPCnSigmaProton(0x0),
  fh2PIDPionTOFnSigmaProton(0x0),
  fh2PIDPionBayesPion(0x0),
  fh2PIDPionBayesKaon(0x0),
  fh2PIDPionBayesProton(0x0),
  fh2PIDKaonTPCnSigmaPion(0x0),
  fh2PIDKaonTOFnSigmaPion(0x0),
  fh2PIDKaonTPCnSigmaKaon(0x0),
  fh2PIDKaonTOFnSigmaKaon(0x0),
  fh2PIDKaonTPCnSigmaProton(0x0),
  fh2PIDKaonTOFnSigmaProton(0x0),
  fh2PIDKaonBayesPion(0x0),
  fh2PIDKaonBayesKaon(0x0),
  fh2PIDKaonBayesProton(0x0),
  fh2PIDProtonTPCnSigmaPion(0x0),
  fh2PIDProtonTOFnSigmaPion(0x0),
  fh2PIDProtonTPCnSigmaKaon(0x0),
  fh2PIDProtonTOFnSigmaKaon(0x0),
  fh2PIDProtonTPCnSigmaProton(0x0),
  fh2PIDProtonTOFnSigmaProton(0x0),
  fh2PIDProtonBayesPion(0x0),
  fh2PIDProtonBayesKaon(0x0),
  fh2PIDProtonBayesProton(0x0),

  // phi histograms
  fhPhiCounter(0x0),
  fhPhiMult(0x0),
  fhPhiBGMult(0x0),
  fhPhiInvMass(0x0),
  fhPhiBGInvMass(0x0),
  fhPhiCharge(0x0),
  fhPhiBGCharge(0x0),
  fhPhiPt(0x0),
  fhPhiEta(0x0),
  fhPhiPhi(0x0),

  // V0s histogram
  fhV0sCounter(0x0),
  fhV0sCounterK0s(0x0),
  fhV0sCounterLambda(0x0),
  fhV0sInvMassK0s(0x0),
  fhV0sInvMassLambda(0x0),
  fhV0sCompetingInvMassK0s(0x0),
  fhV0sCompetingInvMassLambda(0x0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::AliAnalysisTaskUniFlow(const char* name) : AliAnalysisTaskSE(name),
  fEventAOD(0x0),
  fPVz(0.0),
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fFlowWeightsFile(0x0),
  fInit(kFALSE),
  fMC(kFALSE),
  fArrayMC(0x0),
  fIndexSampling(0),
  fIndexCentrality(-1),
  fEventCounter(0),
  fNumEventsAnalyse(50),
  fRunNumber(-1),
  fPDGMassPion(TDatabasePDG::Instance()->GetParticle(211)->Mass()),
  fPDGMassKaon(TDatabasePDG::Instance()->GetParticle(321)->Mass()),
  fPDGMassProton(TDatabasePDG::Instance()->GetParticle(2212)->Mass()),
  fPDGMassPhi(TDatabasePDG::Instance()->GetParticle(333)->Mass()),
  fPDGMassK0s(TDatabasePDG::Instance()->GetParticle(310)->Mass()),
  fPDGMassLambda(TDatabasePDG::Instance()->GetParticle(3122)->Mass()),

  // analysis selection
  fRunMode(kFull),
  fAnalType(kAOD),
  fSampling(kTRUE),
  fFillQA(kTRUE),
  // fNumSamples(10),

  // flow related
  fUseFixedMultBins(kFALSE),
  fCutFlowRFPsPtMin(0.2),
  fCutFlowRFPsPtMax(5.0),
  fFlowPOIsPtMin(0.0),
  fFlowPOIsPtMax(15.0),
  fFlowFillWeights(kTRUE),
  fFlowCentMin(0),
  fFlowCentMax(0),
  fFlowWeightsPath(),
  fFlowRunByRunWeights(kTRUE),
  fFlowUseWeights(kFALSE),
  fFlowUse3Dweights(kFALSE),

  // events selection
  fPVtxCutZ(10.0),
  fColSystem(kPPb),
  fMultEstimator(kV0A),
  fTrigger(AliVEvent::kINT7),
  fEventRejectAddPileUp(kTRUE),

  // charged tracks selection
  fCutChargedEtaMax(0.8),
  fCutChargedDCAzMax(0.),
  fCutChargedDCAxyMax(0.),
  fCutChargedTrackFilterBit(96),
  fCutChargedNumTPCclsMin(70),

  // PID tracks selection
  fCutPIDUseAntiProtonOnly(kFALSE),
  fCutPIDnSigmaPionMax(3.),
  fCutPIDnSigmaKaonMax(3.),
  fCutPIDnSigmaProtonMax(3.),
  fCutPIDnSigmaTPCRejectElectron(3.),
  fCutPIDnSigmaCombinedTOFrejection(kTRUE),
  fCutUseBayesPID(kTRUE),
  fCutPIDBayesPionMin(0.8),
  fCutPIDBayesKaonMin(0.8),
  fCutPIDBayesProtonMin(0.8),
  fCutPIDBayesRejectElectron(0.8),
  fCutPIDBayesRejectMuon(0.8),

  // V0s selection
  fCutV0sOnFly(kFALSE),
  fCutV0srejectKinks(kTRUE),
  fCutV0sDaughterNumTPCClsMin(0),
  fCutV0sDaughterNumTPCCrossMin(70),
  fCutV0sDaughterNumTPCFindMin(1),
  fCutV0sDaughterNumTPCClsPIDMin(0),
  fCutV0sDaughterRatioCrossFindMin(0.8),
  fCutV0srefitTPC(kTRUE),
  fCutV0sCrossMassRejection(kTRUE),
  fCutV0sCrossMassCutK0s(0.005),
  fCutV0sCrossMassCutLambda(0.010),
  fCutV0sCPAK0sMin(0.97),
  fCutV0sCPALambdaMin(0.995),
  fCutV0sDCAtoPVMin(0.06),
  fCutV0sDCAtoPVMax(0.),
  fCutV0sDCAtoPVzMax(0.),
  fCutV0sDCADaughtersMin(0.),
  fCutV0sDCADaughtersMax(1.0),
  fCutV0sDecayRadiusMin(0.5),
  fCutV0sDecayRadiusMax(200.0),
  fCutV0sDaughterFilterBit(0),
  fCutV0sDaughterPtMin(0.),
  fCutV0sDaughterPtMax(0.),
  fCutV0sDaughterEtaMax(0.8),
  fCutV0sMotherEtaMax(0.8),
  fCutV0sMotherRapMax(0.),
  fCutV0sArmenterosAlphaK0sMin(0.2),
  fCutV0sArmenterosAlphaLambdaMax(0.),
  fCutV0sInvMassK0sMin(0.4),
  fCutV0sInvMassK0sMax(0.6),
  fCutV0sInvMassLambdaMin(1.08),
  fCutV0sInvMassLambdaMax(1.16),
  fCutV0sNumTauK0sMax(7.46),
  fCutV0sNumTauLambdaMax(3.8),
  fCutV0sK0sPionNumTPCSigmaMax(5.0),
  fCutV0sLambdaPionNumTPCSigmaMax(5.0),
  fCutV0sLambdaProtonNumTPCSigmaMax(5.0),

  // phi selection
  fCutPhiMotherEtaMax(0.8),
  fCutPhiInvMassMin(0.99),
  fCutPhiInvMassMax(1.07),

  // output lists
  fQAEvents(0x0),
  fQACharged(0x0),
  fQAPID(0x0),
  fQAV0s(0x0),
  fQAPhi(0x0),
  fFlowWeights(0x0),

  // flow histograms & profiles
  fhsV0sCandK0s(0x0),
  fhsV0sCandLambda(0x0),
  fhsPhiCandSig(0x0),
  fhsPhiCandBg(0x0),

  // event histograms
  fhEventSampling(0x0),
  fhEventCentrality(0x0),
  fh2EventCentralityNumRefs(0x0),
  fhEventCounter(0x0),
  fhQAEventsfMult32vsCentr(0x0),
  fhQAEventsMult128vsCentr(0x0),
  fhQAEventsfMultTPCvsTOF(0x0),
  fhQAEventsfMultTPCvsESD(0x0),

  // charged histogram
  fhRefsMult(0x0),
  fhRefsPt(0x0),
  fhRefsEta(0x0),
  fhRefsPhi(0x0),
  fpRefsMult(0x0),
  fhChargedCounter(0x0),

  // PID histogram
  fhMCRecoSelectedPionPt(0x0),
  fhMCRecoSelectedTruePionPt(0x0),
  fhMCRecoAllPionPt(0x0),
  fhMCGenAllPionPt(0x0),
  fhMCRecoSelectedKaonPt(0x0),
  fhMCRecoSelectedTrueKaonPt(0x0),
  fhMCRecoAllKaonPt(0x0),
  fhMCGenAllKaonPt(0x0),
  fhMCRecoSelectedProtonPt(0x0),
  fhMCRecoSelectedTrueProtonPt(0x0),
  fhMCRecoAllProtonPt(0x0),
  fhMCGenAllProtonPt(0x0),
  fhPIDPionMult(0x0),
  fhPIDPionPt(0x0),
  fhPIDPionPhi(0x0),
  fhPIDPionEta(0x0),
  fhPIDPionCharge(0x0),
  fhPIDKaonMult(0x0),
  fhPIDKaonPt(0x0),
  fhPIDKaonPhi(0x0),
  fhPIDKaonEta(0x0),
  fhPIDKaonCharge(0x0),
  fhPIDProtonMult(0x0),
  fhPIDProtonPt(0x0),
  fhPIDProtonPhi(0x0),
  fhPIDProtonEta(0x0),
  fhPIDProtonCharge(0x0),
  fh2PIDPionTPCdEdx(0x0),
  fh2PIDPionTOFbeta(0x0),
  fh2PIDKaonTPCdEdx(0x0),
  fh2PIDKaonTOFbeta(0x0),
  fh2PIDProtonTPCdEdx(0x0),
  fh2PIDProtonTOFbeta(0x0),
  fh2PIDPionTPCnSigmaPion(0x0),
  fh2PIDPionTOFnSigmaPion(0x0),
  fh2PIDPionTPCnSigmaKaon(0x0),
  fh2PIDPionTOFnSigmaKaon(0x0),
  fh2PIDPionTPCnSigmaProton(0x0),
  fh2PIDPionTOFnSigmaProton(0x0),
  fh2PIDPionBayesPion(0x0),
  fh2PIDPionBayesKaon(0x0),
  fh2PIDPionBayesProton(0x0),
  fh2PIDKaonTPCnSigmaPion(0x0),
  fh2PIDKaonTOFnSigmaPion(0x0),
  fh2PIDKaonTPCnSigmaKaon(0x0),
  fh2PIDKaonTOFnSigmaKaon(0x0),
  fh2PIDKaonTPCnSigmaProton(0x0),
  fh2PIDKaonTOFnSigmaProton(0x0),
  fh2PIDKaonBayesPion(0x0),
  fh2PIDKaonBayesKaon(0x0),
  fh2PIDKaonBayesProton(0x0),
  fh2PIDProtonTPCnSigmaPion(0x0),
  fh2PIDProtonTOFnSigmaPion(0x0),
  fh2PIDProtonTPCnSigmaKaon(0x0),
  fh2PIDProtonTOFnSigmaKaon(0x0),
  fh2PIDProtonTPCnSigmaProton(0x0),
  fh2PIDProtonTOFnSigmaProton(0x0),
  fh2PIDProtonBayesPion(0x0),
  fh2PIDProtonBayesKaon(0x0),
  fh2PIDProtonBayesProton(0x0),

  // phi histograms
  fhPhiCounter(0x0),
  fhPhiMult(0x0),
  fhPhiBGMult(0x0),
  fhPhiInvMass(0x0),
  fhPhiBGInvMass(0x0),
  fhPhiCharge(0x0),
  fhPhiBGCharge(0x0),
  fhPhiPt(0x0),
  fhPhiEta(0x0),
  fhPhiPhi(0x0),

  // V0s histogram
  fhV0sCounter(0x0),
  fhV0sCounterK0s(0x0),
  fhV0sCounterLambda(0x0),
  fhV0sInvMassK0s(0x0),
  fhV0sInvMassLambda(0x0),
  fhV0sCompetingInvMassK0s(0x0),
  fhV0sCompetingInvMassLambda(0x0)
{
  // particle dependent
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    fVector[iSpec] = 0x0;
    fListFlow[iSpec] = 0x0;
    fProcessSpec[iSpec] = kTRUE;
    fh2Weights[iSpec] = 0x0;
    fh3Weights[iSpec] = 0x0;
    fh3AfterWeights[iSpec] = 0x0;
  }
  fVecFlowTask = std::vector<FlowTask*>(); //

  // Flow vectors
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  {
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    {
      fFlowVecQpos[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecQneg[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecQmid[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecPpos[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecPneg[iHarm][iPower] = TComplex(0,0,kFALSE);
      fFlowVecS[iHarm][iPower] = TComplex(0,0,kFALSE);
    }
  }

  // QA histograms
  for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
  {
    // Event histograms
    fhQAEventsPVz[iQA] = 0x0;
    fhQAEventsNumContrPV[iQA] = 0x0;
    fhQAEventsNumSPDContrPV[iQA] = 0x0;
    fhQAEventsDistPVSPD[iQA] = 0x0;
    fhQAEventsSPDresol[iQA] = 0x0;

    // charged
    fhQAChargedMult[iQA] = 0x0;
    fhQAChargedPt[iQA] = 0x0;
    fhQAChargedEta[iQA] = 0x0;
    fhQAChargedPhi[iQA] = 0x0;
    fhQAChargedCharge[iQA] = 0x0;
    fhQAChargedFilterBit[iQA] = 0x0;
    fhQAChargedNumTPCcls[iQA] = 0x0;
    fhQAChargedDCAxy[iQA] = 0x0;
    fhQAChargedDCAz[iQA] = 0x0;

    // PID
    fhQAPIDTPCstatus[iQA] = 0x0;
    fhQAPIDTOFstatus[iQA] = 0x0;
    fhQAPIDTPCdEdx[iQA] = 0x0;
    fhQAPIDTOFbeta[iQA] = 0x0;
    fh3QAPIDnSigmaTPCTOFPtPion[iQA] = 0x0;
    fh3QAPIDnSigmaTPCTOFPtKaon[iQA] = 0x0;
    fh3QAPIDnSigmaTPCTOFPtProton[iQA] = 0x0;
    fh3QAPIDnSigmaBayesElectron[iQA] = 0x0;
    fh3QAPIDnSigmaBayesMuon[iQA] = 0x0;
    fh3QAPIDnSigmaBayesPion[iQA] = 0x0;
    fh3QAPIDnSigmaBayesKaon[iQA] = 0x0;
    fh3QAPIDnSigmaBayesProton[iQA] = 0x0;

    // V0s
    fhQAV0sMultK0s[iQA] = 0x0;
    fhQAV0sMultLambda[iQA] = 0x0;
  	fhQAV0sRecoMethod[iQA] = 0x0;
		fhQAV0sDCAtoPV[iQA] = 0x0;
		fhQAV0sDCADaughters[iQA] = 0x0;
		fhQAV0sDecayRadius[iQA] = 0x0;
    fhQAV0sDaughterTPCRefit[iQA] = 0x0;
    fhQAV0sDaughterKinks[iQA] = 0x0;
    fhQAV0sDaughterNumTPCCls[iQA] = 0x0;
    fhQAV0sDaughterNumTPCFind[iQA] = 0x0;
    fhQAV0sDaughterNumTPCCrossRows[iQA] = 0x0;
    fhQAV0sDaughterTPCCrossFindRatio[iQA] = 0x0;
    fhQAV0sDaughterNumTPCClsPID[iQA] = 0x0;
    fhQAV0sDaughterPt[iQA] = 0x0;
		fhQAV0sDaughterPhi[iQA] = 0x0;
		fhQAV0sDaughterEta[iQA] = 0x0;
    fhQAV0sDaughterCharge[iQA] = 0x0;
    fhQAV0sDaughterTPCdEdxK0s[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionK0s[iQA] = 0x0;
    fhQAV0sDaughterTPCstatus[iQA] = 0x0;
    fhQAV0sDaughterTOFstatus[iQA] = 0x0;
    fhQAV0sDaughterTPCdEdxLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaProtonLambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaPionALambda[iQA] = 0x0;
    fhQAV0sDaughterNumSigmaProtonALambda[iQA] = 0x0;
    fhQAV0sMotherPt[iQA] = 0x0;
		fhQAV0sMotherPhi[iQA] = 0x0;
		fhQAV0sMotherEta[iQA] = 0x0;
    fhQAV0sMotherCharge[iQA] = 0x0;
		fhQAV0sMotherRapK0s[iQA] = 0x0;
		fhQAV0sMotherRapLambda[iQA] = 0x0;
    fhQAV0sInvMassK0s[iQA] = 0x0;
    fhQAV0sInvMassLambda[iQA] = 0x0;
		fhQAV0sCPAK0s[iQA] = 0x0;
		fhQAV0sCPALambda[iQA] = 0x0;
		fhQAV0sNumTauK0s[iQA] = 0x0;
		fhQAV0sNumTauLambda[iQA] = 0x0;
		fhQAV0sArmenterosK0s[iQA] = 0x0;
		fhQAV0sArmenterosLambda[iQA] = 0x0;
		fhQAV0sArmenterosALambda[iQA] = 0x0;
  }

  // defining input/output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());
  DefineOutput(7, TList::Class());
  DefineOutput(8, TList::Class());
  DefineOutput(9, TList::Class());
  DefineOutput(10, TList::Class());
  DefineOutput(11, TList::Class());
  DefineOutput(12, TList::Class());
  DefineOutput(13, TList::Class());
  DefineOutput(14, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::~AliAnalysisTaskUniFlow()
{
  // destructor
  // if(fPIDCombined)
  // {
  //   delete fPIDCombined;
  // }

  // clearing vectors before deleting
  ClearVectors();

  // deleting FlowPart vectors (containers)
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    if(fVector[iSpec]) delete fVector[iSpec];
    if(fListFlow[iSpec]) delete fListFlow[iSpec];
  }

  for(Size_t i(0); i < fVecFlowTask.size(); ++i) { delete fVecFlowTask.at(i); }

  // deleting output lists
  if(fFlowWeights) delete fFlowWeights;
  if(fQAEvents) delete fQAEvents;
  if(fQACharged) delete fQACharged;
  if(fQAPID) delete fQAPID;
  if(fQAPhi) delete fQAPhi;
  if(fQAV0s) delete fQAV0s;
}
//_____________________________________________________________________________
const char* AliAnalysisTaskUniFlow::GetSpeciesName(PartSpecies species)
{
  const char* name;

  switch(species)
  {
    case kRefs: name = "Refs"; break;
    case kCharged: name = "Charged"; break;
    case kPion: name = "Pion"; break;
    case kKaon: name = "Kaon"; break;
    case kProton: name = "Proton"; break;
    case kK0s: name = "K0s"; break;
    case kLambda: name = "Lambda"; break;
    case kPhi: name = "Phi"; break;
    default: name = "Unknown";
  }

  return name;
}
//_____________________________________________________________________________
const char* AliAnalysisTaskUniFlow::GetSpeciesLabel(PartSpecies species)
{
  const char* label;

  switch(species)
  {
    case kRefs: label = "RFP"; break;
    case kCharged: label = "h^{#pm}"; break;
    case kPion: label = "#pi^{#pm}"; break;
    case kKaon: label = "K^{#pm}"; break;
    case kProton: label = "p(#bar{p})"; break;
    case kK0s: label = "K^{0}_{S}"; break;
    case kLambda: label = "#Lambda(#bar{#Lambda})"; break;
    case kPhi: label = "#phi"; break;
    default: label = "NA";
  }

  return label;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListParameters()
{
  // lists all task parameters
  // *************************************************************
  AliInfo("Listing all AliAnalysisTaskUniFlow parameters");
  printf("   -------- Analysis task ---------------------------------------\n");
  printf("      fRunMode: (RunMode) %d\n",    fRunMode);
  printf("      fAnalType: (AnalType) %d\n",    fAnalType);
  printf("      fSampling: (Bool_t) %s\n",    fSampling ? "kTRUE" : "kFALSE");
  printf("      fFillQA: (Bool_t) %s\n",    fFillQA ? "kTRUE" : "kFALSE");
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) { printf("      fProcessSpec[k%s]: (Bool_t) %s\n",   GetSpeciesName(PartSpecies(iSpec)), fProcessSpec[iSpec] ? "kTRUE" : "kFALSE"); }
  printf("   -------- Flow related ----------------------------------------\n");
  printf("      fCutFlowRFPsPtMin: (Double_t) %g (GeV/c)\n",    fCutFlowRFPsPtMin);
  printf("      fCutFlowRFPsPtMax: (Double_t) %g (GeV/c)\n",    fCutFlowRFPsPtMax);
  printf("      fFlowPOIsPtMin: (Double_t) %g (GeV/c)\n",    fFlowPOIsPtMin);
  printf("      fFlowPOIsPtMax: (Double_t) %g (GeV/c)\n",    fFlowPOIsPtMax);
  printf("      fFlowCentMin: (Int_t) %d (GeV/c)\n",    fFlowCentMin);
  printf("      fFlowCentMax: (Int_t) %d (GeV/c)\n",    fFlowCentMax);
  printf("      fFlowUseWeights: (Bool_t) %s\n",    fFlowUseWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowRunByRunWeights: (Bool_t) %s\n",    fFlowRunByRunWeights ? "kTRUE" : "kFALSE");
  printf("      fFlowWeightsPath: (TString) '%s' \n",    fFlowWeightsPath.Data());
  printf("   -------- Events ----------------------------------------------\n");
  printf("      fColSystem: (ColSystem) %d\n",    fColSystem);
  printf("      fTrigger: (Short_t) %d\n",    fTrigger);
  printf("      fMultEstimator: (MultiEst) '%s' (%d)\n",    GetMultiEstimatorLabel(fMultEstimator), fMultEstimator);
  printf("      fPVtxCutZ: (Double_t) %g (cm)\n",    fPVtxCutZ);
  printf("   -------- Charged tracks --------------------------------------\n");
  printf("      fCutChargedTrackFilterBit: (UInt) %d\n",    fCutChargedTrackFilterBit);
  printf("      fCutChargedNumTPCclsMin: (UShort_t) %d\n",    fCutChargedNumTPCclsMin);
  printf("      fCutChargedEtaMax: (Double_t) %g\n",    fCutChargedEtaMax);
  printf("      fCutChargedDCAzMax: (Double_t) %g (cm)\n",    fCutChargedDCAzMax);
  printf("      fCutChargedDCAxyMax: (Double_t) %g (cm)\n",    fCutChargedDCAxyMax);
  printf("   -------- PID (pi,K,p) tracks ---------------------------------\n");
  printf("      fCutPIDUseAntiProtonOnly: (Bool_t) %s\n",  fCutPIDUseAntiProtonOnly ? "kTRUE" : "kFALSE");
  printf("      fCutPIDnSigmaCombinedTOFrejection: (Bool_t) %s\n",  fCutPIDnSigmaCombinedTOFrejection ? "kTRUE" : "kFALSE");
  printf("      fCutPIDnSigmaPionMax: (Float_t) %g\n",    fCutPIDnSigmaPionMax);
  printf("      fCutPIDnSigmaKaonMax: (Float_t) %g\n",    fCutPIDnSigmaKaonMax);
  printf("      fCutPIDnSigmaProtonMax: (Float_t) %g\n",    fCutPIDnSigmaProtonMax);
  printf("      fCutPIDnSigmaTPCRejectElectron: (Float_t) %g\n",    fCutPIDnSigmaTPCRejectElectron);
  printf("      fCutUseBayesPID: (Bool_t) %s\n",    fCutUseBayesPID ? "kTRUE" : "kFALSE");
  printf("      fCutPIDBayesPionMin: (Double_t) %g\n",    fCutPIDBayesPionMin);
  printf("      fCutPIDBayesKaonMin: (Double_t) %g\n",    fCutPIDBayesKaonMin);
  printf("      fCutPIDBayesProtonMin: (Double_t) %g\n",    fCutPIDBayesProtonMin);
  printf("      SetPIDBayesRejectElectron: (Double_t) %g\n",    fCutPIDBayesRejectElectron);
  printf("      SetPIDBayesRejectMuon: (Double_t) %g\n",    fCutPIDBayesRejectMuon);
  printf("   -------- Phi candidates --------------------------------------\n");
  printf("      fCutPhiMotherEtaMax: (Double_t) %g\n",    fCutPhiMotherEtaMax);
  printf("      fCutPhiInvMassMin: (Double_t) %g\n",    fCutPhiInvMassMin);
  printf("      fCutPhiInvMassMax: (Double_t) %g\n",    fCutPhiInvMassMax);
  printf("   -------- V0s candidates --------------------------------------\n");
  printf("      fCutV0sOnFly: (Bool_t) %s\n",    fCutV0sOnFly ? "kTRUE" : "kFALSE");
  printf("      fCutV0srefitTPC: (Bool_t) %s\n",     fCutV0srefitTPC ? "kTRUE" : "kFALSE");
  printf("      fCutV0srejectKinks: (Bool_t) %s\n",     fCutV0srejectKinks ? "kTRUE" : "kFALSE");
  printf("      fCutV0sDaughterNumTPCClsMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCClsMin);
  printf("      fCutV0sDaughterNumTPCClsPIDMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCClsPIDMin);
  printf("      fCutV0sDaughterNumTPCCrossMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCCrossMin);
  printf("      fCutV0sDaughterNumTPCFindMin: (UShort_t) %d\n",     fCutV0sDaughterNumTPCFindMin);
  printf("      fCutV0sDaughterRatioCrossFindMin: (Double_t) %g\n",     fCutV0sDaughterRatioCrossFindMin);
  printf("      fCutV0sCrossMassRejection: (Bool_t) %s\n",     fCutV0sCrossMassRejection ? "kTRUE" : "kFALSE");
  printf("      fCutV0sCrossMassCutK0s: (Double_t) %g (GeV/c2)\n",     fCutV0sCrossMassCutK0s);
  printf("      fCutV0sCrossMassCutLambda: (Double_t) %g (GeV/c2)\n",     fCutV0sCrossMassCutLambda);
  printf("      fCutV0sDCAtoPVMin: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVMin);
  printf("      fCutV0sDCAtoPVMax: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVMax);
  printf("      fCutV0sDCAtoPVzMax: (Double_t) %g (cm)\n",    fCutV0sDCAtoPVzMax);
  printf("      fCutV0sDCADaughtersMin: (Double_t) %g (cm)\n",    fCutV0sDCADaughtersMin);
  printf("      fCutV0sDCADaughtersMax: (Double_t) %g (cm)\n",    fCutV0sDCADaughtersMax);
  printf("      fCutV0sDecayRadiusMin: (Double_t) %g (cm)\n",    fCutV0sDecayRadiusMin);
  printf("      fCutV0sDecayRadiusMax: (Double_t) %g (cm)\n",    fCutV0sDecayRadiusMax);
  printf("      fCutV0sDaughterFilterBit: (UInt) %d\n",    fCutV0sDaughterFilterBit);
  printf("      fCutV0sDaughterPtMin: (Double_t) %g (GeV/c)\n",    fCutV0sDaughterPtMin);
  printf("      fCutV0sDaughterPtMax: (Double_t) %g (GeV/c)\n",    fCutV0sDaughterPtMax);
  printf("      fCutV0sDaughterEtaMax: (Double_t) %g ()\n",    fCutV0sDaughterEtaMax);
  printf("      fCutV0sMotherEtaMax: (Double_t) %g ()\n",    fCutV0sMotherEtaMax);
  printf("      fCutV0sMotherRapMax: (Double_t) %g ()\n",    fCutV0sMotherRapMax);
  printf("      fCutV0sCPAK0sMin: (Double_t) %g ()\n",    fCutV0sCPAK0sMin);
  printf("      fCutV0sCPALambdaMin: (Double_t) %g ()\n",    fCutV0sCPALambdaMin);
  printf("      fCutV0sNumTauK0sMax: (Double_t) %g (c*tau)\n",    fCutV0sNumTauK0sMax);
  printf("      fCutV0sNumTauLambdaMax: (Double_t) %g (c*tau)\n",    fCutV0sNumTauLambdaMax);
  printf("      fCutV0sInvMassK0sMin: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassK0sMin);
  printf("      fCutV0sInvMassK0sMax: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassK0sMax);
  printf("      fCutV0sInvMassLambdaMin: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassLambdaMin);
  printf("      fCutV0sInvMassLambdaMax: (Double_t) %g (GeV/c2)\n",    fCutV0sInvMassLambdaMax);
  printf("      fCutV0sArmenterosAlphaK0sMin: (Double_t) %g (alpha)\n",    fCutV0sArmenterosAlphaK0sMin);
  printf("      fCutV0sArmenterosAlphaLambdaMax: (Double_t) %g (alpha)\n",    fCutV0sArmenterosAlphaLambdaMax);
  printf("      fCutV0sK0sPionNumTPCSigmaMax: (Float_t) %g (n*sigma)\n",    fCutV0sK0sPionNumTPCSigmaMax);
  printf("      fCutV0sLambdaPionNumTPCSigmaMax: (Float_t) %g (n*sigma)\n",    fCutV0sLambdaPionNumTPCSigmaMax);
  printf("      fCutV0sLambdaProtonNumTPCSigmaMax: (Float_t) %g (n*sigma)\n",    fCutV0sLambdaProtonNumTPCSigmaMax);
  printf("=====================================================================\n\n");

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ClearVectors()
{
  // Properly clear all particle vectors (if exists)
  // NOTE: should be called at the end of each event & before vectors deleting
  // *************************************************************

  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    std::vector<AliVTrack*>* vector = fVector[iSpec];
    if(!vector) { continue; }
    if(iSpec == kK0s || iSpec == kLambda || iSpec == kPhi) { for(auto part = vector->begin(); part != vector->end(); ++part) { delete *part; } }
    vector->clear();
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::InitializeTask()
{
  // called once on beginning of task (within UserCreateOutputObjects method)
  // check if task parameters are specified and valid
  // returns kTRUE if succesfull
  // *************************************************************
  AliInfo("Checking task setting");
  if(fAnalType != kAOD)
  {
    AliFatal("Analysis type: not kAOD (not implemented for ESDs yet)! Terminating!");
    return kFALSE;
  }

  // checking PID response
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fPIDResponse = inputHandler->GetPIDResponse();
  if(!fPIDResponse)
  {
    AliFatal("AliPIDResponse object not found! Terminating!");
    return kFALSE;
  }

  fPIDCombined = new AliPIDCombined();
  if(!fPIDCombined)
  {
    AliFatal("AliPIDCombined object not found! Terminating!");
    return kFALSE;
  }
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetSelectedSpecies(5); // all particle species
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF); // setting TPC + TOF mask

  // checking the fFlowNumHarmonicsMax, fFlowNumWeightPowersMax dimensions of p,Q,S vectors
  if(fFlowNumWeightPowersMax < 5) { AliFatal("Low range of flow vector weight dimension! Not enought for <4>!"); return kFALSE; }
  // TODO fFlowNumHarmonicsMax;

  if(fSampling && fNumSamples == 0)
  {
    AliFatal("Sampling used, but number of samples is 0! Terminating!");
    return kFALSE;
  }

  // checking cut setting
  AliInfo("Checking task parameters setting conflicts (ranges, etc)");
  if(fCutFlowRFPsPtMin > 0. && fCutFlowRFPsPtMax > 0. && fCutFlowRFPsPtMin > fCutFlowRFPsPtMax)
  {
    AliFatal("Cut: RFPs Pt range wrong! Terminating!");
    return kFALSE;
  }
  if(fCutV0sInvMassK0sMin > fCutV0sInvMassK0sMax || fCutV0sInvMassK0sMin < 0. || fCutV0sInvMassK0sMax < 0.)
  {
    AliFatal("Cut: InvMass (K0s) range wrong! Terminating! ");
    return kFALSE;
  }
  if(fCutV0sInvMassLambdaMin > fCutV0sInvMassLambdaMax || fCutV0sInvMassLambdaMin < 0. || fCutV0sInvMassLambdaMax < 0.)
  {
    AliFatal("Cut: InvMass (Lambda) range wrong! Terminating!");
    return kFALSE;
  }

  // setting fFlowCentMax according to estimator : 0-100 for %'s and 0-150 for CHARGED
  if(fMultEstimator == kRFP) { fFlowCentMin = 0; fFlowCentMax = 200; }
  else { fFlowCentMin = 0; fFlowCentMax = 100; }

  // increasing fFlowCentMax+1 (just to be sure)
  fFlowCentMax++;

  // checking for weights source file
  if(fFlowUseWeights && !fFlowWeightsPath.EqualTo(""))
  {
    fFlowWeightsFile = TFile::Open(Form("%s",fFlowWeightsPath.Data()));
    if(!fFlowWeightsFile) { AliFatal("Flow weights file not found! Terminating!"); return kFALSE; }
    if(!LoadWeights(kTRUE)) { AliFatal("Initial flow weights not loaded! Terminating!"); return kFALSE; }
  }

  AliInfo("Preparing particle containers (std::vectors)");
  // creating particle vectors & reserving capacity in order to avoid memory re-allocation

  Int_t iReserve = 0;
  switch(fColSystem)
  {
    case kPP :
      iReserve = 150;
      break;

    case kPPb :
      iReserve = 200;
      break;

    default :
      iReserve = 300;
      break;
  }

  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    if(!fProcessSpec[iSpec]) { continue; }
    fVector[iSpec] = new std::vector<AliVTrack*>();
    fVector[iSpec]->reserve(iReserve);
  }

  AliInfo("Initialization succesfull!");
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::UserExec(Option_t *)
{
  // main method called for each event (event loop)
  // *************************************************************

  // check if initialization succesfull (done within UserCreateOutputObjects())
  if(!fInit) { AliFatal("Something went wrong : task not initialized!"); return; }

  // local event counter check: if running in test mode, it runs until the 50 events are succesfully processed
  if(fRunMode == kTest && fEventCounter >= fNumEventsAnalyse) { return; }

  // event selection
  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fEventAOD) { return; }

  // loading array with MC particles
  if(fMC)
  {
    fArrayMC = (TClonesArray*) fEventAOD->FindListObject("mcparticles");
    if(!fArrayMC) { AliFatal("TClonesArray with MC particle not found!"); return; }
  }

  // "valid" events before selection
  fhEventCounter->Fill("Input",1);

  // Fill event QA BEFORE cuts
  if(fFillQA) { FillEventsQA(0); }

  if(!IsEventSelected()) { return; }
  fhEventCounter->Fill("Event OK",1);
  // deprecenated selected by HHTF
  // if( (fColSystem == kPP || fColSystem == kPPb) && !IsEventSelected_oldsmall2016() ) { return; }

  // Filter charged (& Refs) particles to evaluate event multiplcity / N_RFPs
  // NB: clear charged vectors because it might keep particles from previous event (not happen for other species)
  fVector[kRefs]->clear();
  fVector[kCharged]->clear();
  FilterCharged();

  // checking if there is at least 5 particles: needed to "properly" calculate correlations
  if(fVector[kRefs]->empty() || fVector[kRefs]->size() < 5) { return; }
  fhEventCounter->Fill("#RPFs OK",1);

  // estimate centrality & assign indexes (centrality/percentile, sampling, ...)
  fIndexCentrality = GetCentralityIndex();
  if(fIndexCentrality < 0) { return; }
  if(fFlowCentMin > 0 && fIndexCentrality < fFlowCentMin) { return; }
  if(fFlowCentMax > 0 && fIndexCentrality > fFlowCentMax) { return; }
  fhEventCounter->Fill("Multiplicity OK",1);

  // here events are selected
  fhEventCounter->Fill("Selected",1);

  // event sampling
  fIndexSampling = GetSamplingIndex();

  // extract PV-z for weights
  fPVz = fEventAOD->GetPrimaryVertex()->GetZ();

  // Fill event QA AFTER cuts
  if(fFillQA) { FillEventsQA(1); }

  if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton] || fProcessSpec[kPhi]) { FilterPID(); }
  if(fProcessSpec[kK0s] || fProcessSpec[kLambda]) { FilterV0s(); }
  if(fProcessSpec[kPhi]) { FilterPhi(); }

  // processing of selected event
  Bool_t bProcessed = CalculateFlow();

  // should be cleared at the end of processing especially for reconstructed
  // particles (Phi, V0s) because here new AliPicoTracks are created
  ClearVectors();

  // extracting run number here to store run number from previous event (for current run number use info in AliAODEvent)
  fRunNumber = fEventAOD->GetRunNumber();

  if(!bProcessed) return;

  // posting data (mandatory)
  Int_t i = 0;
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) { PostData(++i, fListFlow[iSpec]); }
  PostData(++i, fQAEvents);
  PostData(++i, fQACharged);
  PostData(++i, fQAPID);
  PostData(++i, fQAV0s);
  PostData(++i, fQAPhi);
  PostData(++i, fFlowWeights);

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsEventSelected()
{
  // Event selection for pp & p-Pb collision recorded in Run 2 using AliEventCuts
  // return kTRUE if event passes all criteria, kFALSE otherwise
  // *************************************************************

  // Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  if(!(fSelectMask & fTrigger)) { return kFALSE; }

  // events passing physics && trigger selection
  fhEventCounter->Fill("Physics selection OK",1);

  // events passing AliEventCuts selection
  if(!fEventCuts.AcceptEvent(fEventAOD)) { return kFALSE; }
  fhEventCounter->Fill("EventCuts OK",1);

  // Additional pile-up rejection cuts for LHC15o dataset
  if(fColSystem == kPbPb && fEventRejectAddPileUp && IsEventRejectedAddPileUp()) { return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsEventSelected_oldsmall2016()
{
  // Event selection for small system collision recorder in Run 2 year 2016
  // pp (LHC16kl...), pPb (LHC16rqts)
  // return kTRUE if event passes all criteria, kFALSE otherwise
  // *************************************************************

  AliWarning("Using old / deprecenated event selection for small systems!");

  // Physics selection (trigger)
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  if(!(fSelectMask & fTrigger)) { return kFALSE; }

  // events passing physics selection
  fhEventCounter->Fill("Physics selection OK",1);

  // primary vertex selection
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;
  fhEventCounter->Fill("PV OK",1);

  // SPD vertex selection
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(fEventAOD->GetPrimaryVertexSPD());

  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;
  fhEventCounter->Fill("SPD Vtx OK",1);

  // PileUp rejection included in Physics selection
  // but with values for high mult pp (> 5 contrib) => for low ones: do manually (> 3 contrib)

  /*
  if(fTrigger == 0 && fAOD->IsPileupFromSPD(3,0.8) )
  {
    return kFALSE;
  }
  */

  //fhEventCounter->Fill("Pileup SPD OK",1);

  // pileup rejection from multivertexer
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5);
  utils.SetMaxPlpChi2MV(5);
  utils.SetMinWDistMV(15);
  utils.SetCheckPlpFromDifferentBCMV(kFALSE);
  Bool_t isPileupFromMV = utils.IsPileUpMV(fEventAOD);

  if(isPileupFromMV) return kFALSE;
  fhEventCounter->Fill("Pileup MV OK",1);

  // if(fRejectOutOfBunchPU) // out-of-bunch rejection (provided by Christian)
  // {
  //   //out-of-bunch 11 BC
  //   if (utils.IsOutOfBunchPileUp(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //   fhEventCounter->Fill("OOBPU OK",1);
  //
  //   if (utils.IsSPDClusterVsTrackletBG(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //
  //   fhEventCounter->Fill("SPDClTrBG OK",1);
  //
  //   // SPD pileup
  //   if (utils.IsPileUpSPD(fEventAOD))
  //   {
  //     return kFALSE;
  //   }
  //
  //   fhEventCounter->Fill("SPDPU OK",1);
  // }

  //fhEventCounter->Fill("Utils OK",1);

  // cutting on PV z-distance
  const Double_t aodVtxZ = vtx->GetZ();
  if( TMath::Abs(aodVtxZ) > fPVtxCutZ )
  {
    return kFALSE;
  }
  fhEventCounter->Fill("PV #it{z} OK",1);

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsEventRejectedAddPileUp()
{
  // Check for additional pile-up rejection in Run 2 Pb-Pb collisions (15o, 17n)
  // based on multiplicity correlations
  // ***************************************************************************

  Bool_t bIs17n = kFALSE;
  Bool_t bIs15o = kFALSE;

  Int_t iRunNumber = fEventAOD->GetRunNumber();
  if(iRunNumber >= 244824 && iRunNumber <= 246994) { bIs15o = kTRUE; }
  else if(iRunNumber == 280235 || iRunNumber == 20234) { bIs17n = kTRUE; }
  else { return kFALSE; }

  // recounting multiplcities
  Int_t multTPC32 = 0;
  Int_t multTPC128 = 0;
  Int_t multTOF = 0;
  Int_t multESD = 0;
  Int_t multTrk = 0;
  Double_t multESDTPCdif = 0.0;
  Double_t v0Centr = 0.0;

  const Int_t nTracks = fEventAOD->GetNumberOfTracks();
  for(Int_t it(0); it < nTracks; it++)
  {
    AliAODTrack* track = (AliAODTrack*) fEventAOD->GetTrack(it);
    if(!track) { continue; }

    if(track->TestFilterBit(32))
    {
      multTPC32++;
      if(TMath::Abs(track->GetTOFsignalDz()) <= 10.0 && track->GetTOFsignal() >= 12000.0 && track->GetTOFsignal() <= 25000.0) { multTOF++; }
      if((TMath::Abs(track->Eta())) < fCutChargedEtaMax && (track->GetTPCNcls() >= fCutChargedNumTPCclsMin) && (track->Pt() >= fCutFlowRFPsPtMin) && (track->Pt() < fCutFlowRFPsPtMax)) { multTrk++; }
    }

    if(track->TestFilterBit(128)) { multTPC128++; }
  }

  if(bIs17n)
  {
    multESDTPCdif = multESD - (6.6164 + 3.64583*multTPC128 + 0.000126397*multTPC128*multTPC128);
    if(multESDTPCdif > 1000) { return kTRUE; }
    if( ((AliAODHeader*) fEventAOD->GetHeader())->GetRefMultiplicityComb08() < 0) { return kTRUE; }
  }

  if(bIs15o)
  {
    multESDTPCdif = multESD - 3.38*multTPC128;
    if(multESDTPCdif > 500) { return kTRUE; }

    TF1* fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) < fMultTOFLowCut->Eval(Double_t (multTPC32))) { return kTRUE; }

    TF1* fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) > fMultTOFHighCut->Eval(Double_t (multTPC32))) { return kTRUE; }

    AliMultSelection* multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
    v0Centr = multSelection->GetMultiplicityPercentile("V0M");
    TF1* fMultCentLowCut = new TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCentLowCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
    if(Double_t(multTrk) < fMultCentLowCut->Eval(v0Centr)) { return kTRUE; }
  }

  // QA Plots
  fhQAEventsfMult32vsCentr->Fill(v0Centr, multTrk);
  fhQAEventsMult128vsCentr->Fill(v0Centr, multTPC128);
  fhQAEventsfMultTPCvsTOF->Fill(multTPC32, multTOF);
  fhQAEventsfMultTPCvsESD->Fill(multTPC128, multESD);

  // fCentralityDis->Fill(centrV0);
  // fV0CentralityDis->Fill(cent);
  //
  // fCentralityV0MCL1->Fill(v0Centr, cl1Centr);
  // fCentralityV0MCL0->Fill(v0Centr, cl0Centr);
  // fCentralityCL0CL1->Fill(cl0Centr, cl1Centr);

  return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::LoadWeights(Bool_t init)
{
  // (Re-) Loading of flow vector weights
  // ***************************************************************************
  if(!fFlowWeightsFile) { AliError("File with flow weights not found!"); return kFALSE; }

  TList* listFlowWeights = 0x0;
  if(init || !fFlowRunByRunWeights)
  {
    // information about current run is unknown in Initialization(); load only "averaged" weights
    AliInfo("Loading initial GF weights (run-averaged)");
    listFlowWeights = (TList*) fFlowWeightsFile->Get("weights");
    if(!listFlowWeights) { AliError("TList with flow weights not found."); return kFALSE; }
  }
  else
  {
    listFlowWeights = (TList*) fFlowWeightsFile->Get(Form("%d",fRunNumber));
    if(!listFlowWeights) { AliError(Form("TList with flow weights (run %d) not found.",fRunNumber)); return kFALSE; }
  }

  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    if(!fProcessSpec[iSpec]) { continue; }

    if(fFlowUse3Dweights)
    {
      fh3Weights[iSpec] = (TH3D*) listFlowWeights->FindObject(Form("%s3D",GetSpeciesName(PartSpecies(iSpec))));
      if(!fh3Weights[iSpec]) { AliError(Form("Weight (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    }
    else
    {
      fh2Weights[iSpec] = (TH2D*) listFlowWeights->FindObject(GetSpeciesName(PartSpecies(iSpec)));
      if(!fh2Weights[iSpec]) { AliError(Form("Weight 2D (%s) not found",GetSpeciesName(PartSpecies(iSpec)))); return kFALSE; }
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillEventsQA(const Short_t iQAindex)
{
  // Filling various QA plots related with event selection
  // *************************************************************

  if(iQAindex == 1)
  {
    fhEventCentrality->Fill(fIndexCentrality);
    fh2EventCentralityNumRefs->Fill(fIndexCentrality,fVector[kRefs]->size());
    fpRefsMult->Fill(fIndexCentrality,fVector[kRefs]->size(),1.0);
    fhEventSampling->Fill(fIndexCentrality,fIndexSampling);
  }

  const AliAODVertex* aodVtx = fEventAOD->GetPrimaryVertex();
  const Double_t dVtxZ = aodVtx->GetZ();
  const Int_t iNumContr = aodVtx->GetNContributors();
  const AliAODVertex* spdVtx = fEventAOD->GetPrimaryVertexSPD();
  const Int_t iNumContrSPD = spdVtx->GetNContributors();
  const Double_t spdVtxZ = spdVtx->GetZ();

  fhQAEventsPVz[iQAindex]->Fill(dVtxZ);
  fhQAEventsNumContrPV[iQAindex]->Fill(iNumContr);
  fhQAEventsNumSPDContrPV[iQAindex]->Fill(iNumContrSPD);
  fhQAEventsDistPVSPD[iQAindex]->Fill(TMath::Abs(dVtxZ - spdVtxZ));

  // // event / physics selection criteria
  // AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  // UInt_t fSelectMask = inputHandler->IsEventSelected();
  //
  // if( fSelectMask& AliVEvent::kINT7 ) { fQAEventsTriggerSelection[iQAindex]->Fill("kINT7",1); }
  // else if (fSelectMask& AliVEvent::kHighMultV0) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultV0",1); }
  // else if (fSelectMask& AliVEvent::kHighMultSPD) { fQAEventsTriggerSelection[iQAindex]->Fill("kHighMultSPD",1); }
  // else { fQAEventsTriggerSelection[iQAindex]->Fill("Other",1); }

  // SPD vertexer resolution
  Double_t cov[6] = {0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  fhQAEventsSPDresol[iQAindex]->Fill(zRes);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterCharged()
{
  // Filtering input charged tracks for POIs (stored in fVector[kCharged]) or RFPs (fVector[kRefs])
  // If track passes all requirements its pointer is pushed to relevant vector container
  // *************************************************************

  Int_t iNumTracks = fEventAOD->GetNumberOfTracks();
  if(iNumTracks < 1) return;

  Int_t iNumRefs = 0;
  for(Int_t iTrack(0); iTrack < iNumTracks; iTrack++)
  {
    AliAODTrack* track = static_cast<AliAODTrack*>(fEventAOD->GetTrack(iTrack));
    if(!track) continue;

    if(fFillQA) FillQACharged(0,track); // QA before selection

    if(!IsChargedSelected(track)) continue;

    fVector[kCharged]->push_back(track);
    if(fFillQA) FillQACharged(1,track); // QA after selection

    if(fFlowFillWeights) { fh3Weights[kCharged]->Fill(track->Phi(),track->Eta(),fPVz); }
    if(fFlowUseWeights)
    {
      Double_t weight = 1.0;
      if(fFlowUse3Dweights){ weight = fh3Weights[kCharged]->GetBinContent( fh3Weights[kCharged]->FindFixBin(track->Eta(),track->Phi(),fPVz)); }
      else { weight = fh2Weights[kCharged]->GetBinContent( fh2Weights[kCharged]->FindFixBin(track->Eta(),track->Phi())); }
      fh3AfterWeights[kCharged]->Fill(track->Phi(),track->Eta(),fPVz,weight);
    }

    if(fMC)
    {
      AliAODMCParticle* trackMC = GetMCParticle(track->GetLabel());
      if(trackMC)
      {
        Int_t iPDG = TMath::Abs(trackMC->GetPdgCode());

        // filling info about all (i.e. before selection) reconstructed PID particles
        if(iPDG == 211) { fhMCRecoAllPionPt->Fill(track->Pt()); }
        if(iPDG == 321) { fhMCRecoAllKaonPt->Fill(track->Pt()); }
        if(iPDG == 2212) { fhMCRecoAllProtonPt->Fill(track->Pt()); }
      }
    }

    // Checking if selected track is eligible for Ref. flow
    if(!IsWithinRefs(track)) continue;

    // track is used for Ref. flow
    fVector[kRefs]->push_back(track);
    iNumRefs++;
    FillQARefs(1,track);
    if(fFlowFillWeights) { fh3Weights[kRefs]->Fill(track->Phi(),track->Eta(),fPVz); }
    if(fFlowUseWeights)
    {
      Double_t weight = 1.0;
      if(fFlowUse3Dweights) {weight = fh3Weights[kRefs]->GetBinContent( fh3Weights[kRefs]->FindFixBin(track->Eta(),track->Phi(),fPVz)); }
      else { weight = fh2Weights[kRefs]->GetBinContent( fh2Weights[kRefs]->FindFixBin(track->Eta(),track->Phi()) ); }
      fh3AfterWeights[kRefs]->Fill(track->Phi(),track->Eta(),fPVz,weight);
    }

  }

  // fill QA charged multiplicity
  fhRefsMult->Fill(iNumRefs);

  if(fFillQA)
  {
    fhQAChargedMult[0]->Fill(fEventAOD->GetNumberOfTracks());
    fhQAChargedMult[1]->Fill(fVector[kCharged]->size());
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsChargedSelected(const AliAODTrack* track)
{
  // Selection of charged track
  // returns kTRUE if track pass all requirements, kFALSE otherwise
  // *************************************************************
  if(!track) return kFALSE;
  fhChargedCounter->Fill("Input",1);

  // pt
  if(track->Pt() >= fFlowPOIsPtMax) return kFALSE;
  fhChargedCounter->Fill("Pt",1);

  // pseudorapidity (eta)
  if(fCutChargedEtaMax > 0. && TMath::Abs(track->Eta()) >= fCutChargedEtaMax) return kFALSE;
  fhChargedCounter->Fill("Eta",1);

  // filter bit
  if( !track->TestFilterBit(fCutChargedTrackFilterBit) ) return kFALSE;
  fhChargedCounter->Fill("FB",1);

  // number of TPC clusters (additional check for not ITS-standalone tracks)
  if( track->GetTPCNcls() < fCutChargedNumTPCclsMin && fCutChargedTrackFilterBit != 2) return kFALSE;
  fhChargedCounter->Fill("#TPC-Cls",1);

  // track DCA coordinates
  // note AliAODTrack::XYZAtDCA() works only for constrained tracks
  Double_t dTrackXYZ[3] = {0.,0.,0.};
  Double_t dVertexXYZ[3] = {0.,0.,0.};
  Double_t dDCAXYZ[3] = {0.,0.,0.};
  if( fCutChargedDCAzMax > 0. || fCutChargedDCAxyMax > 0.)
  {
    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(!vertex) return kFALSE; // event does not have a PV

    track->GetXYZ(dTrackXYZ);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++) { dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i]; }
  }

  if(fCutChargedDCAzMax > 0. && TMath::Abs(dDCAXYZ[2]) > fCutChargedDCAzMax) return kFALSE;
  fhChargedCounter->Fill("DCA-z",1);

  if(fCutChargedDCAxyMax > 0. && TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]) > fCutChargedDCAxyMax) return kFALSE;
  fhChargedCounter->Fill("DCA-xy",1);


  // track passing all criteria
  fhChargedCounter->Fill("Selected",1);
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsWithinRefs(const AliAODTrack* track)
{
  // Checking if (preselected) track fulfills criteria for RFPs
  // NOTE: This is not a standalone selection, but complementary check for IsChargedSelected()
  // It is used to selecting RFPs out of selected charged tracks
  // OR for estimating autocorrelations for Charged & PID particles
  // *************************************************************
  if(fCutFlowRFPsPtMin > 0.0 && track->Pt() < fCutFlowRFPsPtMin) return kFALSE;
  if(fCutFlowRFPsPtMax > 0.0 && track->Pt() > fCutFlowRFPsPtMax) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillSparseCand(THnSparse* sparse, AliVTrack* track)
{
  // Fill sparse histogram for inv. mass distribution of candidates (V0s,Phi)
  // *************************************************************

  if(!sparse) { Error("THnSparse not valid!","FillSparseCand"); return; }
  if(!track) { Error("Track not valid!","FillSparseCand"); return; }

  Double_t dValues[SparseCand::kDim] = {0};
  dValues[SparseCand::kCent] = fIndexCentrality;
  dValues[SparseCand::kInvMass] = track->M();
  dValues[SparseCand::kPt] = track->Pt();
  dValues[SparseCand::kEta] = track->Eta();
  sparse->Fill(dValues);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQARefs(const Short_t iQAindex, const AliAODTrack* track)
{
  // Filling various QA plots related to RFPs subset of charged track selection
  // *************************************************************

  if(!track) return;
  if(iQAindex == 0) return; // NOTE implemented only for selected RFPs

  fhRefsPt->Fill(track->Pt());
  fhRefsEta->Fill(track->Eta());
  fhRefsPhi->Fill(track->Phi());

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQACharged(const Short_t iQAindex, const AliAODTrack* track)
{
  // Filling various QA plots related to charged track selection
  // *************************************************************
  if(!track) return;

  // filter bit testing
  for(Short_t i(0); i < 32; i++)
  {
    if(track->TestFilterBit(TMath::Power(2.,i)))
      fhQAChargedFilterBit[iQAindex]->Fill(i);
  }

  // track charge
  fhQAChargedCharge[iQAindex]->Fill(track->Charge());

  // number of TPC clusters
  fhQAChargedNumTPCcls[iQAindex]->Fill(track->GetTPCNcls());

  // track DCA
  Double_t dDCAXYZ[3] = {-999., -999., -999.};
  const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
  if(vertex)
  {
    Double_t dTrackXYZ[3] = {-999., -999., -999.};
    Double_t dVertexXYZ[3] = {-999., -999., -999.};

    track->GetXYZ(dTrackXYZ);
    vertex->GetXYZ(dVertexXYZ);

    for(Short_t i(0); i < 3; i++)
      dDCAXYZ[i] = dTrackXYZ[i] - dVertexXYZ[i];
  }
  fhQAChargedDCAxy[iQAindex]->Fill(TMath::Sqrt(dDCAXYZ[0]*dDCAXYZ[0] + dDCAXYZ[1]*dDCAXYZ[1]));
  fhQAChargedDCAz[iQAindex]->Fill(dDCAXYZ[2]);

  // kinematics
  fhQAChargedPt[iQAindex]->Fill(track->Pt());
  fhQAChargedPhi[iQAindex]->Fill(track->Phi());
  fhQAChargedEta[iQAindex]->Fill(track->Eta());

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterV0s()
{
  // Filtering input V0s candidates (K0s, (Anti)Lambda)
  // If track passes all requirements as defined in IsV0sSelected() (and species dependent one)
  // the relevant properties (pT, eta, phi,mass,species) are stored in a new AliPicoTrack
  // and pushed to relevant vector container.
  // *************************************************************

  Int_t iNumK0sSelected = 0;  // counter for selected K0s candidates
  Int_t iNumLambdaSelected = 0; // counter for selected Lambda candidates
  Int_t iNumALambdaSelected = 0; // counter for selected Anti-Lambda candidates

  Int_t iNumV0s = fEventAOD->GetNumberOfV0s();
  if(iNumV0s < 1) return;

  for(Int_t iV0(0); iV0 < iNumV0s; iV0++)
  {
    AliAODv0* v0 = static_cast<AliAODv0*>(fEventAOD->GetV0(iV0));
    if(!v0) continue;

    if(fFillQA) FillQAV0s(0,v0); // QA BEFORE selection

    if(IsV0Selected(v0))
    {
      Bool_t bIsK0s = IsV0aK0s(v0);
      Short_t iIsLambda = IsV0aLambda(v0);

      if(fFillQA && (bIsK0s || iIsLambda != 0))
      FillQAV0s(1,v0,bIsK0s,iIsLambda); // QA AFTER selection

      if(bIsK0s)
      {
        iNumK0sSelected++;
        fhV0sCounter->Fill("K^{0}_{S}",1);
        fhV0sInvMassK0s->Fill(v0->MassK0Short(),v0->MassLambda());

        AliPicoTrack* pico = new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassK0Short());
        fVector[kK0s]->push_back(pico);
        FillSparseCand(fhsV0sCandK0s, pico);

        if(fFlowFillWeights) { fh3Weights[kK0s]->Fill(v0->Phi(),v0->Eta(),fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3Weights[kK0s]->GetBinContent( fh3Weights[kK0s]->FindFixBin(v0->Eta(),v0->Phi(),fPVz) ); }

          else { weight = fh2Weights[kK0s]->GetBinContent( fh2Weights[kK0s]->FindFixBin(v0->Eta(),v0->Phi()) ); }
          fh3AfterWeights[kK0s]->Fill(v0->Phi(),v0->Eta(),fPVz,weight);
        }
      }

      if(iIsLambda == 1) // lambda
      {
        iNumLambdaSelected++;
        fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
        fhV0sInvMassLambda->Fill(v0->MassK0Short(),v0->MassLambda());

        AliPicoTrack* pico = new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassLambda());
        fVector[kLambda]->push_back(pico);
        FillSparseCand(fhsV0sCandLambda, pico);

        if(fFlowFillWeights) { fh3Weights[kLambda]->Fill(v0->Phi(),v0->Eta(),fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3Weights[kLambda]->GetBinContent( fh3Weights[kLambda]->FindFixBin(v0->Eta(),v0->Phi(),fPVz) ); }
          else { weight = fh2Weights[kLambda]->GetBinContent( fh2Weights[kLambda]->FindFixBin(v0->Eta(),v0->Phi()) ); }
          fh3AfterWeights[kLambda]->Fill(v0->Phi(),v0->Eta(),fPVz,weight);
        }
      }

      if(iIsLambda == -1) // anti-lambda
      {
        iNumALambdaSelected++;
        fhV0sCounter->Fill("#Lambda/#bar{#Lambda}",1);
        fhV0sInvMassLambda->Fill(v0->MassK0Short(),v0->MassAntiLambda());

        AliPicoTrack* pico = new AliPicoTrack(v0->Pt(),v0->Eta(),v0->Phi(),v0->Charge(),0,0,0,0,0,0,v0->MassAntiLambda());
        fVector[kLambda]->push_back(pico);
        FillSparseCand(fhsV0sCandLambda, pico);

        if(fFlowFillWeights) { fh3Weights[kLambda]->Fill(v0->Phi(),v0->Eta(),fPVz); }
        if(fFlowUseWeights)
        {
          Double_t weight = 1.0;
          if(fFlowUse3Dweights) { weight = fh3Weights[kLambda]->GetBinContent( fh3Weights[kLambda]->FindFixBin(v0->Eta(),v0->Phi(), fPVz) );}
          else { weight = fh2Weights[kLambda]->GetBinContent( fh2Weights[kLambda]->FindFixBin(v0->Eta(),v0->Phi()) ); }
          fh3AfterWeights[kLambda]->Fill(v0->Phi(),v0->Eta(),fPVz,weight);
        }
      }

      if(bIsK0s && iIsLambda != 0)
      fhV0sCounter->Fill("K^{0}_{S} && #Lambda/#bar{#Lambda}",1);
    }
  }

  // fill QA charged multiplicity
  if(fFillQA)
  {
    fhQAV0sMultK0s[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultLambda[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultALambda[0]->Fill(fEventAOD->GetNumberOfV0s());
    fhQAV0sMultK0s[1]->Fill(iNumK0sSelected);
    fhQAV0sMultLambda[1]->Fill(iNumLambdaSelected);
    fhQAV0sMultALambda[1]->Fill(iNumALambdaSelected);
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0aK0s(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for K0s candidates
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;
  fhV0sCounterK0s->Fill("Input",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0. && ( TMath::Abs(v0->RapK0Short()) > fCutV0sMotherRapMax ) ) return kFALSE;
  fhV0sCounterK0s->Fill("#it{y}",1);

  // inv. mass window
  Double_t dMass = v0->MassK0Short();
  if( dMass < fCutV0sInvMassK0sMin || dMass > fCutV0sInvMassK0sMax ) return kFALSE;
  fhV0sCounterK0s->Fill("InvMass",1);

  // cosine of pointing angle (CPA)
  if(fCutV0sCPAK0sMin > 0.)
  {
    Double_t dCPA = v0->CosPointingAngle(fEventAOD->GetPrimaryVertex());
    if( dCPA < fCutV0sCPAK0sMin ) return kFALSE;
  }
  fhV0sCounterK0s->Fill("CPA",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaK0sMin > 0.)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm < (fCutV0sArmenterosAlphaK0sMin * TMath::Abs(dAlpha))) return kFALSE;
  }
  fhV0sCounterK0s->Fill("Armenteros-Podolanski",1);

  // proper life-time
  if( fCutV0sNumTauK0sMax > 0. )
  {
    Double_t dPrimVtxCoor[3] = {0.}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0.}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0.}; // decay vector coor {xyz}
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 3; i++) { dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i]; }

    // implementation in xy plane
    // Double_t dPropLife = ( (fPDGMassK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLife = ( (fPDGMassK0s / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]) );
    if( dPropLife > (fCutV0sNumTauK0sMax * 2.68) ) return kFALSE;
  }
  fhV0sCounterK0s->Fill("c#tau",1);

  // Daughter PID
  if(fCutV0sK0sPionNumTPCSigmaMax > 0.)
  {
    const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
    const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

    if(!HasTrackPIDTPC(daughterPos) || !HasTrackPIDTPC(daughterNeg)) return kFALSE;

    if (daughterPos->GetTPCsignalN() < fCutV0sDaughterNumTPCClsPIDMin || daughterNeg->GetTPCsignalN() < fCutV0sDaughterNumTPCClsPIDMin) return kFALSE;
    Float_t nSigmaPiPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterPos, AliPID::kPion));
    Float_t nSigmaPiNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(daughterNeg, AliPID::kPion));
    if(nSigmaPiPos > fCutV0sK0sPionNumTPCSigmaMax || nSigmaPiNeg > fCutV0sK0sPionNumTPCSigmaMax) return kFALSE;
  }
  fhV0sCounterK0s->Fill("Daughters PID",1);

  // competing V0 rejection based on InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassLambda = v0->MassLambda();
    Double_t dMassALambda = v0->MassAntiLambda();

    // K0s candidate is within 10 MeV of (Anti)Lambda InvMass physSelTask
    if(TMath::Abs(dMassLambda - fPDGMassLambda) < fCutV0sCrossMassCutK0s)
    {
      // in Lambda peak
      fhV0sCompetingInvMassK0s->Fill(dMass,dMassLambda);
      return kFALSE;
    }

    if(TMath::Abs(dMassALambda - fPDGMassLambda) < fCutV0sCrossMassCutK0s)
    {
      // in Anti-Lambda peak
      fhV0sCompetingInvMassK0s->Fill(dMass,dMassALambda);
      return kFALSE;
    }
  }
  fhV0sCounterK0s->Fill("Competing InvMass",1);

  // passing all criteria
  fhV0sCounterK0s->Fill("Selected",1);
  return kTRUE;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::IsV0aLambda(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // specific for Lambda candidates
  // return 0 if candidate does not fullfill any Lambda or Anti-Lambda requirements;
  // return 1 if a candidate fulfill all Lambda requirements;
  // return -1 if a candidate fullfill all Anti-Lambda requirements;
  // return 2 if a candidate fulfill all both Lambda & Anti-Lambda requirements
  // *************************************************************
  if(!v0) return 0;
  fhV0sCounterLambda->Fill("Input",1);

  // rapidity selection
  if(fCutV0sMotherRapMax > 0. && ( TMath::Abs(v0->RapLambda()) > fCutV0sMotherRapMax) ) return 0;
  fhV0sCounterLambda->Fill("#it{y}",1);

  // particle species dependent
  Bool_t bIsLambda = kFALSE;
  Bool_t bIsALambda = kFALSE;

  // inv. mass window
  Double_t dMassLambda = v0->MassLambda();
  Double_t dMassALambda = v0->MassAntiLambda();
  if( dMassLambda > fCutV0sInvMassLambdaMin && dMassLambda < fCutV0sInvMassLambdaMax) { bIsLambda = kTRUE; }
  if( dMassALambda > fCutV0sInvMassLambdaMin && dMassALambda < fCutV0sInvMassLambdaMax) { bIsALambda = kTRUE; }

  if(!bIsLambda && !bIsALambda) return 0;
  fhV0sCounterLambda->Fill("InvMass",1);

  // cosine of pointing angle (CPA)
  if( fCutV0sCPALambdaMin > 0. )
  {
    Double_t dCPA = v0->CosPointingAngle(fEventAOD->GetPrimaryVertex());
    if( dCPA < fCutV0sCPALambdaMin ) return 0;
  }
  fhV0sCounterLambda->Fill("CPA",1);

  // Armenteros-Podolaski plot
  if(fCutV0sArmenterosAlphaLambdaMax > 0.)
  {
    Double_t dPtArm = v0->PtArmV0();
    Double_t dAlpha = v0->AlphaV0();
    if(dPtArm > (fCutV0sArmenterosAlphaLambdaMax * TMath::Abs(dAlpha))) return 0;
  }
  fhV0sCounterLambda->Fill("Armenteros-Podolanski",1);

  // // Armenteros-Podolanski for candidates fullfilling both Lambda and Anti-Lambda selection
  // if(bIsLambda && bIsALambda)
  // {
  //   if(v0->AlphaV0() < 0.) bIsLambda = kFALSE;
  //   if(v0->AlphaV0() > 0.) bIsALambda = kFALSE;
  // }

  // proper life-time
  if( fCutV0sNumTauLambdaMax > 0. )
  {
    Double_t dPrimVtxCoor[3] = {0.}; // primary vertex position {x,y,z}
    Double_t dSecVtxCoor[3] = {0.}; // secondary vertex position {x,y,z}
    Double_t dDecayCoor[3] = {0.}; // decay vector coor {xyz}
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    primVtx->GetXYZ(dPrimVtxCoor);
    v0->GetSecondaryVtx(dSecVtxCoor);

    for(Int_t i(0); i < 3; i++) { dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i]; }

    // implementation in xy plane
    // Double_t dPropLife = ( (fPDGMassLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    Double_t dPropLife = ( (fPDGMassLambda / (v0->P() + 1e-10) ) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1] + dDecayCoor[2]*dDecayCoor[2]) );
    if( dPropLife > (fCutV0sNumTauLambdaMax * 7.89) ) return 0;
  }
  fhV0sCounterLambda->Fill("c#tau",1);

  // daughter PID of Lambda Candidates
  if( fCutV0sLambdaProtonNumTPCSigmaMax > 0. || fCutV0sLambdaPionNumTPCSigmaMax > 0. )
  {
    const AliAODTrack* trackDaughterPos = (AliAODTrack*) v0->GetDaughter(0); // positive charge
    const AliAODTrack* trackDaughterNeg = (AliAODTrack*) v0->GetDaughter(1); // negative charge

    Bool_t bIsPosOK = HasTrackPIDTPC(trackDaughterPos);
    Bool_t bIsNegOK = HasTrackPIDTPC(trackDaughterNeg);

    Float_t dSigmaPos = 999.;
    Float_t dSigmaNeg = 999.;

    if(fCutV0sLambdaPionNumTPCSigmaMax > 0.) // check pions
    {
      if(bIsPosOK) dSigmaPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kPion)); else dSigmaPos = 999.;
      if(bIsNegOK) dSigmaNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kPion)); else dSigmaNeg = 999.;

      if(bIsLambda && (!bIsNegOK || dSigmaNeg > fCutV0sLambdaPionNumTPCSigmaMax) ) bIsLambda = kFALSE;
      if(bIsALambda && (!bIsPosOK || dSigmaPos > fCutV0sLambdaPionNumTPCSigmaMax) ) bIsALambda = kFALSE;
    }

    if(fCutV0sLambdaProtonNumTPCSigmaMax > 0.) // check protons
    {
      if(bIsPosOK) dSigmaPos = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterPos, AliPID::kProton)); else dSigmaPos = 999.;
      if(bIsNegOK) dSigmaNeg = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackDaughterNeg, AliPID::kProton)); else dSigmaNeg = 999.;

      if(bIsLambda && (!bIsPosOK || dSigmaPos > fCutV0sLambdaProtonNumTPCSigmaMax) ) bIsLambda = kFALSE;
      if(bIsALambda && (!bIsNegOK || dSigmaNeg > fCutV0sLambdaProtonNumTPCSigmaMax) ) bIsALambda = kFALSE;
    }

    if(!bIsLambda && !bIsALambda) return 0;
  }
  fhV0sCounterLambda->Fill("Daughter PID",1);

  // Lambda(AntiLamda) candidate is within fCutV0sCrossMassCutLambda of K0s InvMass
  if(fCutV0sCrossMassRejection)
  {
    Double_t dMassK0s = v0->MassK0Short();
    if(TMath::Abs(dMassK0s - fPDGMassK0s) < fCutV0sCrossMassCutLambda)
    {
      if(bIsLambda) fhV0sCompetingInvMassLambda->Fill(dMassK0s,dMassLambda);
      if(bIsALambda) fhV0sCompetingInvMassLambda->Fill(dMassK0s,dMassALambda);
      return 0;
    }
  }
  fhV0sCounterLambda->Fill("Competing InvMass",1);

  // passing all criteria
  fhV0sCounterLambda->Fill("Selected",1);

  if(bIsLambda && bIsALambda) { fhV0sCounterLambda->Fill("#Lambda && #bar{#Lambda}",1); return 2; } // both Lambda & Anti-Lambda
  if(bIsLambda) { fhV0sCounterLambda->Fill("only #Lambda",1); return 1; } // only Lambda
  if(bIsALambda) { fhV0sCounterLambda->Fill("only #bar{#Lambda}",1); return -1; } // only Anti-Lambda
  return 0; // neither
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskUniFlow::GetRapidity(Double_t mass, Double_t Pt, Double_t Eta)
{
    Double_t rapid = TMath::Log( (TMath::Sqrt(mass*mass + Pt*Pt*TMath::CosH(Eta)*TMath::CosH(Eta)) + Pt*TMath::SinH(Eta)) / TMath::Sqrt(mass*mass + Pt*Pt) );
    return rapid;
}
//_____________________________________________________________________________
AliAODMCParticle* AliAnalysisTaskUniFlow::GetMCParticle(Int_t label)
{
  if(!fArrayMC) { AliError("fArrayMC not found!"); return 0x0; }
  if(label < 0) { /*AliWarning("MC label negative");*/ return 0x0; }

  AliAODMCParticle* mcTrack = (AliAODMCParticle*) fArrayMC->At(label);
  if(!mcTrack) { AliError("Corresponding MC track not found!"); return 0x0; }
  return mcTrack;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::IsV0Selected(const AliAODv0* v0)
{
  // Topological reconstruction and selection of V0 candidates
  // common for both K0s and (Anti)-Lambdas
  // return kTRUE if a candidate fulfill all requirements, kFALSE otherwise
  // *************************************************************
  if(!v0) return kFALSE;
  fhV0sCounter->Fill("Input",1);

  const AliAODTrack* daughterPos = (AliAODTrack*) v0->GetDaughter(0);
  const AliAODTrack* daughterNeg = (AliAODTrack*) v0->GetDaughter(1);

  // daughter track check
  if(!daughterPos || !daughterNeg) return kFALSE;
  fhV0sCounter->Fill("Daughters OK",1);

  // acceptance checks
  if(fCutV0sMotherEtaMax > 0. && TMath::Abs(v0->Eta()) >= fCutV0sMotherEtaMax ) return kFALSE;
  if(v0->Pt() >= fFlowPOIsPtMax) return kFALSE;
  fhV0sCounter->Fill("Mother acceptance",1);

  if(fCutV0sDaughterPtMin > 0. && (daughterPos->Pt() <= fCutV0sDaughterPtMin  || daughterNeg->Pt() <= fCutV0sDaughterPtMin) ) return kFALSE;
  if(fCutV0sDaughterPtMax > 0. && (daughterPos->Pt() >= fCutV0sDaughterPtMax  || daughterNeg->Pt() >= fCutV0sDaughterPtMax) ) return kFALSE;
  if(fCutV0sDaughterEtaMax > 0. && ( (TMath::Abs(daughterNeg->Eta()) >= fCutV0sDaughterEtaMax) || (TMath::Abs(daughterPos->Eta()) >= fCutV0sDaughterEtaMax) ) ) return kFALSE;
  fhV0sCounter->Fill("Daughter acceptance",1);

  // daughters & mother charge checks
  if(v0->Charge() != 0) return kFALSE;
  if(daughterPos->Charge() == daughterNeg->Charge()) return kFALSE;
  if( (TMath::Abs(daughterPos->Charge()) != 1) || (TMath::Abs(daughterNeg->Charge()) != 1) ) return kFALSE;
  fhV0sCounter->Fill("Charge",1);

  // reconstruction method: online (on-the-fly) OR offline
  if(v0->GetOnFlyStatus() != fCutV0sOnFly) return kFALSE;
  fhV0sCounter->Fill("Reconstruction method",1);

  // TPC refit
  if(fCutV0srefitTPC && ( !daughterPos->IsOn(AliAODTrack::kTPCrefit) || !daughterNeg->IsOn(AliAODTrack::kTPCrefit) ) ) return kFALSE;
  fhV0sCounter->Fill("TPC refit",1);

  // filter bit
  if( fCutV0sDaughterFilterBit > 0 && (!daughterPos->TestFilterBit(fCutV0sDaughterFilterBit) || !daughterNeg->TestFilterBit(fCutV0sDaughterFilterBit) ) ) return kFALSE;
  fhV0sCounter->Fill("Daughter FB",1);

  // Kinks
  const AliAODVertex* prodVtxDaughterPos = (AliAODVertex*) daughterPos->GetProdVertex(); // production vertex of the positive daughter track
  const AliAODVertex* prodVtxDaughterNeg = (AliAODVertex*) daughterNeg->GetProdVertex(); // production vertex of the negative daughter track
  if(fCutV0srejectKinks && ( (prodVtxDaughterPos->GetType() == AliAODVertex::kKink ) || (prodVtxDaughterNeg->GetType() == AliAODVertex::kKink ) ) ) return kFALSE;
  fhV0sCounter->Fill("Kinks",1);

  // Daughter track quality
  if(daughterPos->GetTPCNcls() < fCutV0sDaughterNumTPCClsMin || daughterNeg->GetTPCNcls() < fCutV0sDaughterNumTPCClsMin) return kFALSE;
  if(daughterPos->GetTPCNCrossedRows() < fCutV0sDaughterNumTPCCrossMin || daughterNeg->GetTPCNCrossedRows() < fCutV0sDaughterNumTPCCrossMin) return kFALSE;
  if(daughterPos->GetTPCNclsF() < fCutV0sDaughterNumTPCFindMin || daughterNeg->GetTPCNclsF() < fCutV0sDaughterNumTPCFindMin) return kFALSE;
  if(fCutV0sDaughterRatioCrossFindMin > -1.)
  {
    if(daughterPos->GetTPCNclsF() < 1 || daughterNeg->GetTPCNclsF() < 1) return kFALSE; // at least 1 findable cls for proper division
    Double_t dRatioCrossFindPos = (Double_t) daughterPos->GetTPCNCrossedRows() / (Double_t) daughterPos->GetTPCNclsF();
    Double_t dRatioCrossFindNeg = (Double_t) daughterNeg->GetTPCNCrossedRows() / (Double_t) daughterNeg->GetTPCNclsF();
    if( dRatioCrossFindPos < fCutV0sDaughterRatioCrossFindMin || dRatioCrossFindNeg < fCutV0sDaughterRatioCrossFindMin) return kFALSE;
  }
  fhV0sCounter->Fill("Daughters track quality",1);

  // Daughters DCA to PV
  Double_t dDCAPosToPV = TMath::Abs(v0->DcaPosToPrimVertex());
  Double_t dDCANegToPV = TMath::Abs(v0->DcaNegToPrimVertex());
  if(fCutV0sDCAtoPVMin > 0. && ( dDCAPosToPV < fCutV0sDCAtoPVMin || dDCANegToPV < fCutV0sDCAtoPVMin ) ) return kFALSE;
  if(fCutV0sDCAtoPVMax > 0. && ( dDCAPosToPV > fCutV0sDCAtoPVMax || dDCANegToPV > fCutV0sDCAtoPVMax ) ) return kFALSE;

  // note AliAODTrack::XYZAtDCA() works only for constrained tracks
  if(fCutV0sDCAtoPVzMax > 0.)
  {
    Double_t dVertexXYZ[3] = {0.};
    Double_t dTrackXYZpos[3] = {0.};
    Double_t dTrackXYZneg[3] = {0.};
    Double_t dDCAXYZpos[3] = {0.};
    Double_t dDCAXYZneg[3] = {0.};

    const AliAODVertex* vertex = fEventAOD->GetPrimaryVertex();
    if(!vertex) return kFALSE; // event does not have a PV

    vertex->GetXYZ(dVertexXYZ);
    daughterPos->GetXYZ(dTrackXYZpos);
    daughterNeg->GetXYZ(dTrackXYZneg);

    for(Short_t i(0); i < 3; i++)
    {
      dDCAXYZpos[i] = dTrackXYZpos[i] - dVertexXYZ[i];
      dDCAXYZneg[i] = dTrackXYZneg[i] - dVertexXYZ[i];
    }

    if( TMath::Abs(dDCAXYZpos[2]) > fCutV0sDCAtoPVzMax || TMath::Abs(dDCAXYZneg[2]) > fCutV0sDCAtoPVzMax ) return kFALSE;
  }
  fhV0sCounter->Fill("DCA to PV",1);

  // Daughter DCA among themselves
  Double_t dDCADaughters = v0->DcaV0Daughters();
  if(fCutV0sDCADaughtersMin > 0. && TMath::Abs(dDCADaughters) < fCutV0sDCADaughtersMin) return kFALSE;
  if(fCutV0sDCADaughtersMax > 0. && TMath::Abs(dDCADaughters) > fCutV0sDCADaughtersMax) return kFALSE;
  fhV0sCounter->Fill("Daughters DCA",1);

  // radius of decay vertex in transverse plane
  Double_t dDecayRadius = v0->RadiusV0();
  if( fCutV0sDecayRadiusMin > 0. && (dDecayRadius < fCutV0sDecayRadiusMin) ) return kFALSE;
  if( fCutV0sDecayRadiusMax > 0. && (dDecayRadius > fCutV0sDecayRadiusMax) ) return kFALSE;
  fhV0sCounter->Fill("Decay radius",1);

  // passing all common criteria
  fhV0sCounter->Fill("Common passed",1);
  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAV0s(const Short_t iQAindex, const AliAODv0* v0, const Bool_t bIsK0s, const Short_t bIsLambda)
{
  // Filling various QA plots related to V0 candidate selection
  // *************************************************************
  // checking mother & daughters
  if(!v0) return;
  AliAODTrack* trackDaughter[2] = {(AliAODTrack*) v0->GetDaughter(0), (AliAODTrack*) v0->GetDaughter(1)};
  if(!trackDaughter[0] || !trackDaughter[1]) return;

  // setting internal flags for Lambdas and Anti-Lambdas
  Bool_t bCandLambda;
  Bool_t bCandAntiLambda;

  switch (bIsLambda)
  {
    case 1:
    {
      bCandLambda = kTRUE;
      bCandAntiLambda = kFALSE;
      break;
    }
    case -1:
    {
      bCandLambda = kFALSE;
      bCandAntiLambda = kTRUE;
      break;
    }
    case 2:
    {
      bCandLambda = kTRUE;
      bCandAntiLambda = kTRUE;
      break;
    }
    default:
    {
      bCandLambda = kFALSE;
      bCandAntiLambda = kFALSE;
    }
  }

  // reconstruction method
  fhQAV0sRecoMethod[iQAindex]->Fill(v0->GetOnFlyStatus());

  // DCA between daughters and PV
  fhQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaPosToPrimVertex());
  fhQAV0sDCAtoPV[iQAindex]->Fill(v0->DcaNegToPrimVertex());

  // Daughter DCA among themselves
  fhQAV0sDCADaughters[iQAindex]->Fill(v0->DcaV0Daughters());

  // charge
  fhQAV0sMotherCharge[iQAindex]->Fill(v0->Charge());

  // radius of decay vertex in transverse plane
  Double_t dSecVtxCoor[3] = {0};
  v0->GetSecondaryVtx(dSecVtxCoor);
  Double_t dDecayRadius = TMath::Sqrt(dSecVtxCoor[0]*dSecVtxCoor[0] + dSecVtxCoor[1]*dSecVtxCoor[1]);
  fhQAV0sDecayRadius[iQAindex]->Fill(dDecayRadius);

  // mother kinematics
  fhQAV0sMotherPt[iQAindex]->Fill(v0->Pt());
  fhQAV0sMotherPhi[iQAindex]->Fill(v0->Phi());
  fhQAV0sMotherEta[iQAindex]->Fill(v0->Eta());

  // proper lifetime preparation (to be filled in particle dependent if scope)
  Double_t dPrimVtxCoor[3] = {0};
  Double_t dDecayCoor[3] = {0};
  AliAODVertex* primVtx2 = fEventAOD->GetPrimaryVertex();
  primVtx2->GetXYZ(dPrimVtxCoor);
  for(Int_t i(0); i < 2; i++)
    dDecayCoor[i] = dSecVtxCoor[i] - dPrimVtxCoor[i];

  // particle dependent
  if(bIsK0s)
  {
    // K0s
    fhQAV0sMotherRapK0s[iQAindex]->Fill(v0->RapK0Short());
    fhQAV0sInvMassK0s[iQAindex]->Fill(v0->MassK0Short());

    // CPA
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    fhQAV0sCPAK0s[iQAindex]->Fill(v0->CosPointingAngle(primVtx));

    // Armenteros-Podolanski
    fhQAV0sArmenterosK0s[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    // proper lifetime
    Double_t dMassPDGK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass();
    Double_t dPropLifeK0s = ( (dMassPDGK0s / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    fhQAV0sNumTauK0s[iQAindex]->Fill(dPropLifeK0s);
  }
  if(bCandLambda || bCandAntiLambda)
  {
    // (Anti)Lambda
    fhQAV0sMotherRapLambda[iQAindex]->Fill(v0->RapLambda());
    fhQAV0sInvMassLambda[iQAindex]->Fill(v0->MassLambda());
    fhQAV0sInvMassLambda[iQAindex]->Fill(v0->MassAntiLambda());

    // CPA
    AliAODVertex* primVtx = fEventAOD->GetPrimaryVertex();
    fhQAV0sCPALambda[iQAindex]->Fill(v0->CosPointingAngle(primVtx));

    // Armenteros-Podolanski
    if(bCandLambda)
      fhQAV0sArmenterosLambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    if(bCandAntiLambda)
      fhQAV0sArmenterosALambda[iQAindex]->Fill(v0->AlphaV0(), v0->PtArmV0());

    // proper lifetime
    Double_t dMassPDGLambda = TDatabasePDG::Instance()->GetParticle(kLambda0)->Mass();
    Double_t dPropLifeLambda = ( (dMassPDGLambda / v0->Pt()) * TMath::Sqrt(dDecayCoor[0]*dDecayCoor[0] + dDecayCoor[1]*dDecayCoor[1]) );
    fhQAV0sNumTauLambda[iQAindex]->Fill(dPropLifeLambda);
  }

  AliPIDResponse::EDetPidStatus pidStatusTPC;
  AliPIDResponse::EDetPidStatus pidStatusTOF;
  UShort_t numTPCfindable = 0;
  Float_t numTPCcrossed = 0;

  // daughters properties
  AliAODVertex* prodVtxDaughter = 0x0;
  for(Short_t i(0); i < 2; i++)
  {
    // TPC refit
    fhQAV0sDaughterTPCRefit[iQAindex]->Fill(trackDaughter[i]->IsOn(AliAODTrack::kTPCrefit));

    // kinks
    prodVtxDaughter = (AliAODVertex*) trackDaughter[i]->GetProdVertex();
    fhQAV0sDaughterKinks[iQAindex]->Fill(prodVtxDaughter->GetType() == AliAODVertex::kKink);

    // track quality
    numTPCcrossed = trackDaughter[i]->GetTPCNCrossedRows();
    numTPCfindable = trackDaughter[i]->GetTPCNclsF();
    fhQAV0sDaughterNumTPCCls[iQAindex]->Fill(trackDaughter[i]->GetTPCNcls());
    fhQAV0sDaughterNumTPCClsPID[iQAindex]->Fill(trackDaughter[i]->GetTPCsignalN());
    fhQAV0sDaughterNumTPCFind[iQAindex]->Fill(numTPCfindable);
    fhQAV0sDaughterNumTPCCrossRows[iQAindex]->Fill(numTPCcrossed);
    if(numTPCfindable > 0.) fhQAV0sDaughterTPCCrossFindRatio[iQAindex]->Fill(numTPCcrossed/numTPCfindable); else fhQAV0sDaughterTPCCrossFindRatio[iQAindex]->Fill(-5.);

    // detector status
    pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[i]);
    pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, trackDaughter[i]);
    fhQAV0sDaughterTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );
    fhQAV0sDaughterTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );

    // daughter kinematics
    fhQAV0sDaughterPt[iQAindex]->Fill(trackDaughter[i]->Pt());
    fhQAV0sDaughterPhi[iQAindex]->Fill(trackDaughter[i]->Phi());
    fhQAV0sDaughterEta[iQAindex]->Fill(trackDaughter[i]->Eta());

    // daughter charge
    fhQAV0sDaughterCharge[iQAindex]->Fill(trackDaughter[i]->Charge());
  }

  AliPIDResponse::EDetPidStatus pidStatusTPCpos;
  AliPIDResponse::EDetPidStatus pidStatusTPCneg;

  // PID checks
  if(fPIDResponse)
  {
    // checking the detector status
    pidStatusTPCpos = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[0]);
    pidStatusTPCneg = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, trackDaughter[1]);

    if(pidStatusTPCpos == AliPIDResponse::kDetPidOk && pidStatusTPCneg == AliPIDResponse::kDetPidOk)
    {
      if(bIsK0s)
      {
        // daughter PID
        fhQAV0sDaughterTPCdEdxK0s[iQAindex]->Fill(trackDaughter[0]->P(), trackDaughter[0]->GetTPCsignal());
        fhQAV0sDaughterTPCdEdxK0s[iQAindex]->Fill(trackDaughter[1]->P(), trackDaughter[1]->GetTPCsignal());

        // Pion PID for daughters
        fhQAV0sDaughterNumSigmaPionK0s[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kPion));
        fhQAV0sDaughterNumSigmaPionK0s[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kPion));
      }

      if(bCandLambda || bCandAntiLambda)
      {
        // daughter PID
        fhQAV0sDaughterTPCdEdxLambda[iQAindex]->Fill(trackDaughter[0]->P(), trackDaughter[0]->GetTPCsignal());
        fhQAV0sDaughterTPCdEdxLambda[iQAindex]->Fill(trackDaughter[1]->P(), trackDaughter[1]->GetTPCsignal());

        if(bCandLambda)
        {
          fhQAV0sDaughterNumSigmaProtonLambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kProton));
          fhQAV0sDaughterNumSigmaPionLambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kPion));
        }

        if(bCandAntiLambda)
        {
          fhQAV0sDaughterNumSigmaProtonALambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[1], AliPID::kProton));
          fhQAV0sDaughterNumSigmaPionALambda[iQAindex]->Fill(v0->Pt(), fPIDResponse->NumberOfSigmasTPC(trackDaughter[0], AliPID::kPion));
        }
      }
    }
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterPhi()
{
  // Reconstruction and filtering of Phi meson candidates out of selected Kaon sample
  // If track passes all requirements, the relevant properties (pT, eta, phi) are stored
  // in FlowPart struct  and pushed to relevant vector container.
  // *************************************************************
  Int_t iNumKaons = (Int_t) fVector[kKaon]->size();
  // check if there are at least 2 selected kaons in event (minimum for phi reconstruction)
  if(iNumKaons < 2) return;

  // start Phi reconstruction
  Int_t iNumBG = 0;
  for(Int_t iKaon1(0); iKaon1 < iNumKaons; iKaon1++)
  {
    AliAODTrack* kaon1 = dynamic_cast<AliAODTrack*>(fVector[kKaon]->at(iKaon1));
    if(!kaon1) continue;

    for(Int_t iKaon2(iKaon1+1); iKaon2 < iNumKaons; iKaon2++)
    {
      AliAODTrack* kaon2 = dynamic_cast<AliAODTrack*>(fVector[kKaon]->at(iKaon2));
      if(!kaon2) continue;

      AliPicoTrack* mother = MakeMother(kaon1,kaon2);
      fhPhiCounter->Fill("Input",1);

      // filling QA BEFORE selection
      if(fFillQA) FillQAPhi(0,mother);

      if(mother->M() < fCutPhiInvMassMin || mother->M() > fCutPhiInvMassMax) { delete mother; continue; }
      fhPhiCounter->Fill("InvMass",1);

      if(mother->Pt() < fFlowPOIsPtMin || mother->Pt() > fFlowPOIsPtMax) { delete mother; continue; }
      fhPhiCounter->Fill("Pt",1);

      if(fCutPhiMotherEtaMax > 0 && TMath::Abs(mother->Eta()) > fCutPhiMotherEtaMax) { delete mother; continue; }
      fhPhiCounter->Fill("Eta",1);

      // mother (phi) candidate passing all criteria
      fhPhiCounter->Fill("Selected",1);

      // filling QA AFTER selection
      if(fFillQA) FillQAPhi(1,mother);

      // filling weights
      if(fFlowFillWeights) { fh3Weights[kPhi]->Fill(mother->Phi(), mother->Eta(), fPVz); }
      if(fFlowUseWeights)
      {
        Double_t weight = 1.0;
        if(fFlowUse3Dweights) { weight =  fh3Weights[kPhi]->GetBinContent(  fh3Weights[kPhi]->FindFixBin(mother->Eta(),mother->Phi(), fPVz) ); }
        else { weight =  fh2Weights[kPhi]->GetBinContent(  fh2Weights[kPhi]->FindFixBin(mother->Eta(),mother->Phi()) ) ; }
        fh3AfterWeights[kPhi]->Fill(mother->Phi(),mother->Eta(),fPVz,weight);
      }

      if(mother->Charge() == 0)
      {
        // opposite-sign combination (signal+background)
        fhPhiCounter->Fill("Unlike-sign",1);
        FillSparseCand(fhsPhiCandSig, mother);
        fVector[kPhi]->push_back(mother);
      }

      if(TMath::Abs(mother->Charge()) == 2)
      {
        // like-sign combination (background)
        FillSparseCand(fhsPhiCandBg, mother);
        iNumBG++;
      }

    } // endfor {iKaon2} : second kaon
  } // endfor {iKaon1} : first Kaon

  // filling multiplicity distribution
  fhPhiMult->Fill(fVector[kPhi]->size());
  fhPhiBGMult->Fill(iNumBG);

  return;
}
//____________________________________________________________________
AliPicoTrack* AliAnalysisTaskUniFlow::MakeMother(const AliAODTrack* part1, const AliAODTrack* part2)
{
  // Reconstructing mother particle from two prongs and fill its properties.
  // return ptr to created mother particle
  // *************************************************************

  if(!part1 || !part2) return 0x0;

  // combining momenta
  TVector3 mom1 = TVector3( part1->Px(), part1->Py(), part1->Pz() );
  TVector3 mom2 = TVector3( part2->Px(), part2->Py(), part2->Pz() );
  TVector3 mom = mom1 + mom2;

  Byte_t iCharge = part1->Charge() + part2->Charge();

  // calculating inv. mass
  Double_t dMass = -999.;
  Double_t dE1 = TMath::Sqrt( mom1.Mag2() + TMath::Power(fPDGMassKaon,2) );
  Double_t dE2 = TMath::Sqrt( mom2.Mag2() + TMath::Power(fPDGMassKaon,2) );

  Double_t dMassSq = TMath::Power((dE1+dE2),2) - mom.Mag2();
  if(dMassSq >= 0.) dMass = TMath::Sqrt(dMassSq);

  return new AliPicoTrack(mom.Pt(),mom.Eta(),mom.Phi(),iCharge,0,0,0,0,0,0,dMass);
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAPhi(const Short_t iQAindex, const AliPicoTrack* part)
{
  if(!part) return;

  if(iQAindex == 0) return; // TODO not implemented (do not know what)

  if(part->Charge() == 0)
  {
    // phi candidate (unlike-sign pair)
    fhPhiInvMass->Fill(part->M());
    fhPhiCharge->Fill(part->Charge());
    fhPhiPt->Fill(part->Pt());
    fhPhiEta->Fill(part->Eta());
    fhPhiPhi->Fill(part->Phi());
  }

  if(TMath::Abs(part->Charge()) == 2)
  {
    // phi candidate (unlike-sign pair)
    fhPhiBGInvMass->Fill(part->M());
    fhPhiBGCharge->Fill(part->Charge());
  }

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FilterPID()
{
  // Filtering input PID tracks (pi,K,p)
  // If track passes all requirements as defined in IsPIDSelected() (and species dependent),
  // the relevant properties (pT, eta, phi) are stored in FlowPart struct
  // and pushed to relevant vector container.
  // return kFALSE if any complications occurs
  // *************************************************************

  for(auto part = fVector[kCharged]->begin(); part != fVector[kCharged]->end(); ++part)
  {
    AliAODTrack* track = static_cast<AliAODTrack*>(*part);
    if(!track) continue;

    if(fFillQA) FillQAPID(0,track,kUnknown);   // filling QA for tracks before selection (but after charged criteria applied)

    // PID tracks are subset of selected charged tracks (same quality requirements)
    if(!IsChargedSelected(track)) continue;

    // PID track selection (return most favourable species)
    PartSpecies species = IsPIDSelected(track);
    if(species != kPion && species != kKaon && species != kProton) { continue; }
    // check if only protons should be used
    if(fCutPIDUseAntiProtonOnly && species == kProton && track->Charge() == 1) { species = kUnknown; }

    fVector[species]->push_back(track);
    if(fFlowFillWeights) { fh3Weights[species]->Fill(track->Phi(), track->Eta(), fPVz); }
    if(fFlowUseWeights)
    {
      Double_t weight = 1.0;
      if(fFlowUse3Dweights) { weight = fh3Weights[species]->GetBinContent( fh2Weights[species]->FindFixBin(track->Eta(),track->Phi(), fPVz)); }
      else { weight = fh2Weights[species]->GetBinContent( fh2Weights[species]->FindFixBin(track->Eta(),track->Phi()) ); }
      fh3AfterWeights[species]->Fill(track->Phi(),track->Eta(),fPVz,weight);
    }

    if(fFillQA) FillQAPID(1,track,species); // filling QA for tracks AFTER selection

    // process MC data
    if(fMC)
    {
      AliAODMCParticle* trackMC = GetMCParticle(track->GetLabel());
      if(trackMC)
      {
        Int_t iPDG = TMath::Abs(trackMC->GetPdgCode());

        switch (species)
        {
          case kPion:
            fhMCRecoSelectedPionPt->Fill(track->Pt());
            if(iPDG == 211) { fhMCRecoSelectedTruePionPt->Fill(track->Pt()); }
          break;
          case kKaon:
            fhMCRecoSelectedKaonPt->Fill(track->Pt());
            if(iPDG == 321) { fhMCRecoSelectedTrueKaonPt->Fill(track->Pt()); }
          break;
          case kProton:
            fhMCRecoSelectedProtonPt->Fill(track->Pt());
            if(iPDG == 2212) { fhMCRecoSelectedTrueProtonPt->Fill(track->Pt()); }
          break;

          default :
            continue;
        }

      }
    }

  }

  fhPIDPionMult->Fill(fVector[kPion]->size());
  fhPIDKaonMult->Fill(fVector[kKaon]->size());
  fhPIDProtonMult->Fill(fVector[kProton]->size());

  return;
}
//_____________________________________________________________________________
AliAnalysisTaskUniFlow::PartSpecies AliAnalysisTaskUniFlow::IsPIDSelected(const AliAODTrack* track)
{
  // Selection of PID tracks (pi,K,p) - track identification
  // Based on fCutUseBayesPID flag, either Bayes PID or nSigma cutting is used
  // returns AliAnalysisTaskUniFlow::PartSpecies enum : kPion, kKaon, kProton if any of this passed kUnknown otherwise
  // *************************************************************

  // checking detector statuses
  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  if(!bIsTPCok) return kUnknown;
  // TODO: TOF check???

  if(fCutUseBayesPID)
  {
    // use Bayesian PID
    Double_t dProbPID[5] = {0.}; // array for Bayes PID probabilities:  0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dProbPID); // filling probabilities to dPropPID array
    Double_t dMaxProb = TMath::MaxElement(5,dProbPID);

    // printf("PID Prob: e %g | mu %g | pi %g | K %g | p %g ||| MAX %g \n",dProbPID[0],dProbPID[1],dProbPID[2],dProbPID[3],dProbPID[4],dMaxProb);

    // electron and mion rejection
    if(dProbPID[0] >= fCutPIDBayesRejectElectron || dProbPID[1] >= fCutPIDBayesRejectMuon) return kUnknown;

    // checking the PID probability
    // TODO: think about: if Pion has maximum probibility < fCutBayesPIDPion, track is rejected -> is it good?
    if(dMaxProb == dProbPID[2] && dProbPID[2] >= fCutPIDBayesPionMin) return kPion;
    if(dMaxProb == dProbPID[3] && dProbPID[3] >= fCutPIDBayesKaonMin) return kKaon;
    if(dMaxProb == dProbPID[4] && dProbPID[4] >= fCutPIDBayesProtonMin) return kProton;
  }
  else
  {
    // use nSigma cuts (based on combination of TPC / TOF nSigma cuts
    const Double_t dPt = track->Pt();

    Float_t dNumSigmaTPC[5] = {-99.,-99.,-99.,-99.,-99.}; // TPC nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
    Float_t dNumSigmaTOF[5] = {-99.,-99.,-99.,-99.,-99.}; // TOF nSigma array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton

    // filling nSigma arrays
    if(bIsTPCok) // should be anyway
    {
      dNumSigmaTPC[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron));
      dNumSigmaTPC[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon));
      dNumSigmaTPC[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
      dNumSigmaTPC[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
      dNumSigmaTPC[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    }

    if(bIsTOFok) // should be anyway
    {
      dNumSigmaTOF[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron));
      dNumSigmaTOF[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon));
      dNumSigmaTOF[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion));
      dNumSigmaTOF[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon));
      dNumSigmaTOF[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton));
    }

    // TPC nSigma cuts
    if(dPt <= 0.4)
    {
      Float_t dMinSigmasTPC = TMath::MinElement(5,dNumSigmaTPC);

      // electron rejection
      if(dMinSigmasTPC == dNumSigmaTPC[0] && dNumSigmaTPC[0] <= fCutPIDnSigmaTPCRejectElectron) return kUnknown;
      if(dMinSigmasTPC == dNumSigmaTPC[2] && dNumSigmaTPC[2] <= fCutPIDnSigmaPionMax) return kPion;
      if(dMinSigmasTPC == dNumSigmaTPC[3] && dNumSigmaTPC[3] <= fCutPIDnSigmaKaonMax) return kKaon;
      if(dMinSigmasTPC == dNumSigmaTPC[4] && dNumSigmaTPC[4] <= fCutPIDnSigmaProtonMax) return kProton;
    }

    // combined TPC + TOF nSigma cuts
    // NB: for testing of efficiency removed the upper limmit for PID
    // if(dPt > 0.4 && dPt < 4.0)
    if(dPt > 0.4)
    {
      Float_t dNumSigmaCombined[5] = {-99.,-99.,-99.,-99.,-99.};

      // discard candidates if no TOF is available if cut is on
      if(fCutPIDnSigmaCombinedTOFrejection && !bIsTOFok) return kUnknown;

      // calculating combined nSigmas
      for(Short_t i(0); i < 5; i++)
      {
        if(bIsTOFok) { dNumSigmaCombined[i] = TMath::Sqrt(dNumSigmaTPC[i]*dNumSigmaTPC[i] + dNumSigmaTOF[i]*dNumSigmaTOF[i]); }
        else { dNumSigmaCombined[i] = dNumSigmaTPC[i]; }
      }

      Float_t dMinSigmasCombined = TMath::MinElement(5,dNumSigmaCombined);

      // electron rejection
      if(dMinSigmasCombined == dNumSigmaCombined[0] && dNumSigmaCombined[0] <= fCutPIDnSigmaTPCRejectElectron) return kUnknown;
      if(dMinSigmasCombined == dNumSigmaCombined[2] && dNumSigmaCombined[2] <= fCutPIDnSigmaPionMax) return kPion;
      if(dMinSigmasCombined == dNumSigmaCombined[3] && dNumSigmaCombined[3] <= fCutPIDnSigmaKaonMax) return kKaon;
      if(dMinSigmasCombined == dNumSigmaCombined[4] && dNumSigmaCombined[4] <= fCutPIDnSigmaProtonMax) return kProton;
    }

    // if(dPt >= 4.0)
    // {
      // NOTE: in this pt range, nSigmaTPC is not enought to distinquish well between species
      // all three values are close to each other -> minimum difference is not applied, just nSigma cut

      // TPC dEdx parametrisation (dEdx - <dEdx>)
      // TODO: TPC dEdx parametrisation cuts
      // if(dPt > 3.)
      //
    // }
  }

  return kUnknown;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillQAPID(const Short_t iQAindex, const AliAODTrack* track, const PartSpecies species)
{
  // Filling various QA plots related to PID (pi,K,p) track selection
  // *************************************************************
  if(!track) return;

  if(!fPIDResponse || !fPIDCombined)
  {
    AliError("AliPIDResponse or AliPIDCombined object not found!");
    return;
  }

  // TPC & TOF statuses & measures
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);

  fhQAPIDTOFstatus[iQAindex]->Fill((Int_t) pidStatusTOF );
  fhQAPIDTPCstatus[iQAindex]->Fill((Int_t) pidStatusTPC );

  Bool_t bIsTPCok = HasTrackPIDTPC(track);
  Bool_t bIsTOFok = HasTrackPIDTOF(track);

  Double_t dNumSigmaTPC[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
  Double_t dNumSigmaTOF[5] = {-11}; // array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton
  Double_t dBayesProb[5] = {-0.1}; // Bayesian probability | array: 0: electron / 1: muon / 2: pion / 3: kaon / 4: proton

  Double_t dTPCdEdx = -5; // TPC dEdx for selected particle
  Double_t dTOFbeta = -0.05; //TOF beta for selected particle

  // filling Bayesian PID probabilities to dBayesProb array
  UInt_t iDetUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, dBayesProb);

  Double_t dP = track->P();
  Double_t dPt = track->Pt();

  // detector status dependent
  if(bIsTPCok)
  {
    dNumSigmaTPC[0] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    dNumSigmaTPC[1] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon);
    dNumSigmaTPC[2] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    dNumSigmaTPC[3] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    dNumSigmaTPC[4] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    dTPCdEdx = track->GetTPCsignal();
    fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), dTPCdEdx);
  }
  else // TPC status not OK
  {
    dNumSigmaTPC[0] = -11.;
    dNumSigmaTPC[1] = -11.;
    dNumSigmaTPC[2] = -11.;
    dNumSigmaTPC[3] = -11.;
    dNumSigmaTPC[4] = -11.;

    fhQAPIDTPCdEdx[iQAindex]->Fill(track->P(), -5.);
  }

  if(bIsTOFok)
  {
    dNumSigmaTOF[0] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
    dNumSigmaTOF[1] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon);
    dNumSigmaTOF[2] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    dNumSigmaTOF[3] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    dNumSigmaTOF[4] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    Double_t dTOF[5];
    track->GetIntegratedTimes(dTOF);
    dTOFbeta = dTOF[0] / track->GetTOFsignal();
    fhQAPIDTOFbeta[iQAindex]->Fill(dP,dTOFbeta);
  }
  else // TOF status not OK
  {
    dNumSigmaTOF[0] = -11.;
    dNumSigmaTOF[1] = -11.;
    dNumSigmaTOF[2] = -11.;
    dNumSigmaTOF[3] = -11.;
    dNumSigmaTOF[4] = -11.;

    fhQAPIDTOFbeta[iQAindex]->Fill(track->P(),-0.05);
  }

  fh3QAPIDnSigmaTPCTOFPtPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],track->Pt());
  fh3QAPIDnSigmaTPCTOFPtKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],track->Pt());
  fh3QAPIDnSigmaTPCTOFPtProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],track->Pt());

  fh3QAPIDnSigmaBayesElectron[iQAindex]->Fill(dNumSigmaTPC[0],dNumSigmaTOF[0],dBayesProb[0]);
  fh3QAPIDnSigmaBayesMuon[iQAindex]->Fill(dNumSigmaTPC[1],dNumSigmaTOF[1],dBayesProb[1]);
  fh3QAPIDnSigmaBayesPion[iQAindex]->Fill(dNumSigmaTPC[2],dNumSigmaTOF[2],dBayesProb[2]);
  fh3QAPIDnSigmaBayesKaon[iQAindex]->Fill(dNumSigmaTPC[3],dNumSigmaTOF[3],dBayesProb[3]);
  fh3QAPIDnSigmaBayesProton[iQAindex]->Fill(dNumSigmaTPC[4],dNumSigmaTOF[4],dBayesProb[4]);

  // species dependent QA
  switch (species)
  {
    case kPion:
      fhPIDPionPt->Fill(track->Pt());
      fhPIDPionPhi->Fill(track->Phi());
      fhPIDPionEta->Fill(track->Eta());
      fhPIDPionCharge->Fill(track->Charge());
      fh2PIDPionTPCdEdx->Fill(dPt,dTPCdEdx);
      fh2PIDPionTOFbeta->Fill(dPt,dTOFbeta);
      fh2PIDPionTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDPionTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDPionTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDPionTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDPionTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDPionTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDPionBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDPionBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDPionBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    case kKaon:
      fhPIDKaonPt->Fill(track->Pt());
      fhPIDKaonPhi->Fill(track->Phi());
      fhPIDKaonEta->Fill(track->Eta());
      fhPIDKaonCharge->Fill(track->Charge());
      fh2PIDKaonTPCdEdx->Fill(dP,dTPCdEdx);
      fh2PIDKaonTOFbeta->Fill(dP,dTOFbeta);
      fh2PIDKaonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDKaonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDKaonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDKaonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDKaonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDKaonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDKaonBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDKaonBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDKaonBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    case kProton:
      fhPIDProtonPt->Fill(track->Pt());
      fhPIDProtonPhi->Fill(track->Phi());
      fhPIDProtonEta->Fill(track->Eta());
      fhPIDProtonCharge->Fill(track->Charge());
      fh2PIDProtonTPCdEdx->Fill(dP,dTPCdEdx);
      fh2PIDProtonTOFbeta->Fill(dP,dTOFbeta);
      fh2PIDProtonTPCnSigmaPion->Fill(dPt,dNumSigmaTPC[2]);
      fh2PIDProtonTOFnSigmaPion->Fill(dPt,dNumSigmaTOF[2]);
      fh2PIDProtonTPCnSigmaKaon->Fill(dPt,dNumSigmaTPC[3]);
      fh2PIDProtonTOFnSigmaKaon->Fill(dPt,dNumSigmaTOF[3]);
      fh2PIDProtonTPCnSigmaProton->Fill(dPt,dNumSigmaTPC[4]);
      fh2PIDProtonTOFnSigmaProton->Fill(dPt,dNumSigmaTOF[4]);
      fh2PIDProtonBayesPion->Fill(dPt,dBayesProb[2]);
      fh2PIDProtonBayesKaon->Fill(dPt,dBayesProb[3]);
      fh2PIDProtonBayesProton->Fill(dPt,dBayesProb[4]);
      break;

    default:
      break;
  }


  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::ProcessFlowTask(FlowTask* task)
{
  if(!task) { AliError("FlowTask does not exists!"); return kFALSE; }
  // task->Print();

  Int_t iNumHarm = task->fiNumHarm;
  Int_t iNumGaps = task->fiNumGaps;

  if(iNumGaps > 1) { AliError("Too many gaps! Not implemented yet!"); return kFALSE; }
  if(iNumHarm > 4) { AliError("Too many harmonics! Not implemented yet!"); return kFALSE; }

  Double_t dGap = -1.0;
  if(iNumGaps > 0) { dGap = task->fdGaps[0]; }

  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    // check if FlowTask should be done for all flow particles (RFP/POI/Both)
    if(task->fPart == FlowTask::kRFP && iSpec != kRefs) { continue; }
    if(task->fPart == FlowTask::kPOI && iSpec == kRefs) { continue; }

    if(iSpec == kRefs) {
      FillRefsVectors(dGap);
      CalculateCorrelations(task, PartSpecies(iSpec));
      continue; // NO need to go further
    }

    // loading (generic) profile to acess axes and bins
    TH1* genProf = (TH1*) fListFlow[iSpec]->FindObject(Form("%s_Pos_sample0",task->fsName.Data()));
    if(!genProf) { AliError(Form("Generic Profile '%s' not found", task->fsName.Data())); return kFALSE; }

    TAxis* axisPt = genProf->GetYaxis(); if(!axisPt) { AliError("Pt axis object not found!"); return kFALSE; }
    Int_t iNumPtBins = axisPt->GetNbins();

    Int_t iNumMassBins = 1;
    TAxis* axisMass = 0x0;

    // check for 'massive' species
    Bool_t bHasMass = kFALSE;
    if(iSpec == kK0s || iSpec == kLambda || iSpec == kPhi) {
      bHasMass = kTRUE;
      axisMass = genProf->GetZaxis();
      if(!axisMass) { AliError("Mass axis object not found!"); return kFALSE; }
      iNumMassBins = axisMass->GetNbins();
    }

    for(Int_t iMass(1); iMass < iNumMassBins+1; ++iMass)
    {
      Double_t dMass = 0.0;
      Double_t dMassLow = 0.0;
      Double_t dMassHigh = 0.0;

      if(bHasMass)
      {
        dMass = axisMass->GetBinCenter(iMass);
        dMassLow = axisMass->GetBinLowEdge(iMass);
        dMassHigh = axisMass->GetBinUpEdge(iMass);
      }

      for(Int_t iPt(1); iPt < iNumPtBins+1; ++iPt)
      {
        Double_t dPt = axisPt->GetBinCenter(iPt);
        Double_t dPtLow = axisPt->GetBinLowEdge(iPt);
        Double_t dPtHigh = axisPt->GetBinUpEdge(iPt);

        // filling POIs (P,S) flow vectors
        FillPOIsVectors(dGap,PartSpecies(iSpec),dPtLow,dPtHigh,dMassLow,dMassHigh);
        CalculateCorrelations(task, PartSpecies(iSpec),dPt,dMass);
      }
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::CalculateCorrelations(FlowTask* task, PartSpecies species, Double_t dPt, Double_t dMass)
{
  if(!task) { AliError("FlowTask does not exists!"); return; }
  if(species >= kUnknown) { AliError(Form("Invalid species: %s!", GetSpeciesName(species))); return; }

  Bool_t bHasGap = task->HasGap();
  Int_t iNumHarm = task->fiNumHarm;
  Bool_t bDiff = kTRUE; if(species == kRefs) { bDiff = kFALSE; }

  TComplex cNom = TComplex(0.0,0.0,kFALSE);
  TComplex cDenom = TComplex(0.0,0.0,kFALSE);
  TComplex cNomNeg = TComplex(0.0,0.0,kFALSE);
  TComplex cDenomNeg = TComplex(0.0,0.0,kFALSE);

  // calculating correlations
  switch(iNumHarm)
  {
    case 2 : {
      if(!bHasGap) { // no gap
        if(bDiff) {
          cDenom = TwoDiff(0,0);
          cNom = TwoDiff(task->fiHarm[0],task->fiHarm[1]);
        }
        else {
          cDenom = Two(0,0);
          cNom = Two(task->fiHarm[0],task->fiHarm[1]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = TwoDiffGapPos(0,0);
          cDenomNeg = TwoDiffGapNeg(0,0);
          cNom = TwoDiffGapPos(task->fiHarm[0],task->fiHarm[1]);
          cNomNeg = TwoDiffGapNeg(task->fiHarm[0],task->fiHarm[1]);
        }
        else {
          cDenom = TwoGap(0,0);
          cNom = TwoGap(task->fiHarm[0],task->fiHarm[1]); }
        }
        break;
      }

    case 3 : {
      if(!bHasGap) { // no gap
        if(bDiff) {
          cDenom = ThreeDiff(0,0,0);
          cNom = ThreeDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
        }
        else {
          cDenom = Three(0,0,0);
          cNom = Three(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = ThreeDiffGapPos(0,0,0);
          cDenomNeg = ThreeDiffGapNeg(0,0,0);
          cNom = ThreeDiffGapPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
          cNomNeg = ThreeDiffGapNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]);
        }
        else {
          AliWarning("ThreeGap() not implemented!");
          return;
          // cDenom = ThreeGap(0,0,0);
          // cNom = ThreeGap(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2]); }
        }
      }
      break;
    }

    case 4 : {
      if(!bHasGap) { // no gap
        if(bDiff) {
          cDenom = FourDiff(0,0,0,0);
          cNom = FourDiff(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
        }
        else {
          cDenom = Four(0,0,0,0);
          cNom = Four(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
        }
      }
      else { // has gap
        if(bDiff) {
          cDenom = FourDiffGapPos(0,0,0,0);
          cDenomNeg = FourDiffGapNeg(0,0,0,0);
          cNom = FourDiffGapPos(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
          cNomNeg = FourDiffGapNeg(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
        }
        else {
          cDenom = FourGap(0,0,0,0);
          cNom = FourGap(task->fiHarm[0],task->fiHarm[1],task->fiHarm[2],task->fiHarm[3]);
        }
      }
      break;
    }

    default:
      return;
  }

  Double_t dNom = cNom.Re();
  Double_t dDenom = cDenom.Re();
  Double_t dNomNeg = cNomNeg.Re();
  Double_t dDenomNeg = cDenomNeg.Re();

  Double_t dValue = 0.0;
  Double_t dValueNeg = 0.0;

  Bool_t bFillPos = kFALSE;
  Bool_t bFillNeg = kFALSE;

  // check which profiles should be filled (POS/NEG)
  if(dDenom > 0.0) { bFillPos = kTRUE; dValue = dNom / dDenom; }
  if(bFillPos && TMath::Abs(dValue) > 1.0) { bFillPos = kFALSE; }

  if(bHasGap && dDenomNeg > 0.0) { bFillNeg = kTRUE; dValueNeg = dNomNeg / dDenomNeg; }
  if(bFillNeg && TMath::Abs(dValueNeg) > 1.0) { bFillNeg = kFALSE; }

  if(!bFillPos && !bFillNeg) { return; } // To save some CPU time

  // Filling corresponding profiles
  switch(species)
  {
    case kRefs:
    {
      if(bFillPos)
      {
        TProfile* prof = (TProfile*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d",task->fsName.Data(),fIndexSampling));
        if(!prof) { AliError(Form("Profile '%s_Pos_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        prof->Fill(fIndexCentrality, dValue, dDenom);
      }
      break;
    }

    case kCharged:
    case kPion:
    case kKaon:
    case kProton:
    {
      if(bFillPos)
      {
        TProfile2D* prof = (TProfile2D*) fListFlow[species]->FindObject(Form("%s_Pos_sample%d",task->fsName.Data(),fIndexSampling));
        if(!prof) { AliError(Form("Profile '%s_Pos_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        prof->Fill(fIndexCentrality, dPt, dValue, dDenom);
      }

      if(bFillNeg)
      {
        TProfile2D* profNeg = (TProfile2D*) fListFlow[species]->FindObject(Form("%s_Neg_sample%d",task->fsName.Data(),fIndexSampling));
        if(!profNeg) { AliError(Form("Profile '%s_Neg_sample%d' not found!", task->fsName.Data(),fIndexSampling)); return; }
        profNeg->Fill(fIndexCentrality, dPt, dValueNeg, dDenomNeg);
      }
      break;
    }


    case kK0s:
    case kLambda:
    case kPhi:
    {
      if(bFillPos)
      {
        TProfile3D* prof = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Pos_sample0",task->fsName.Data()));
        if(!prof) { AliError(Form("Profile '%s_Pos_sample0' not found!", task->fsName.Data())); return; }
        prof->Fill(fIndexCentrality, dPt, dMass, dValue, dDenom);
      }

      if(bFillNeg)
      {
        TProfile3D* profNeg = (TProfile3D*) fListFlow[species]->FindObject(Form("%s_Neg_sample0",task->fsName.Data()));
        if(!profNeg) { AliError(Form("Profile '%s_Neg_sample0' not found!", task->fsName.Data())); return; }
        profNeg->Fill(fIndexCentrality, dPt, dMass, dValueNeg, dDenomNeg);
      }
      break;
    }

    case kUnknown:
      return;
  }

  return;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::CalculateFlow()
{
  // main (envelope) method for flow calculations in selected events
  // returns kFALSE if something failes (with error), kTRUE otherwise
  // *************************************************************

  // if running in kSkipFlow mode, skip the remaining part
  if(fRunMode == kSkipFlow) { fEventCounter++; return kTRUE; }

  // checking the run number for aplying weights & loading TList with weights
  if(fFlowUseWeights && fFlowRunByRunWeights && fRunNumber != fEventAOD->GetRunNumber() && !LoadWeights(kFALSE)) { return kFALSE; }

  // >>>> Using FlowTask <<<<<

  Int_t iNumTasks = fVecFlowTask.size();
  for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
  {
    Bool_t process = ProcessFlowTask(fVecFlowTask.at(iTask));
    if(!process) { AliError("FlowTask processing failed!\n"); fVecFlowTask.at(iTask)->Print(); return kFALSE; }
  }

  fEventCounter++; // counter of processed events

  return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillRefsVectors(const Double_t dGap)
{
  // Filling Q flow vector with RFPs
  // return kTRUE if succesfull (i.e. no error occurs), kFALSE otherwise
  // *************************************************************
  Double_t dEtaGap = dGap;
  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasGap = kFALSE;
  Bool_t bHas3sub = kFALSE;
  if(dEtaGap > -1.0) { bHasGap = kTRUE; }

  TH2D* h2Weights = fh2Weights[kRefs];
  TH3D* h3Weights = fh3Weights[kRefs];
  if(fFlowUseWeights)
  {
    if(fFlowUse3Dweights && !h3Weights) { AliError("Histogram with Refs weights not found."); return; }
    if(!fFlowUse3Dweights && !h2Weights) { AliError("Histogram with Refs weights not found."); return; }
  }

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecQpos);
  ResetFlowVector(fFlowVecQneg);
  if(bHas3sub) { ResetFlowVector(fFlowVecQmid); }

  for (auto part = fVector[kRefs]->begin(); part != fVector[kRefs]->end(); part++)
  {
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();

    Double_t dWeight = 1.0;

    // loading weights if needed
    if(fFlowUseWeights)
    {
      if(fFlowUse3Dweights) { dWeight = h3Weights->GetBinContent(h3Weights->FindFixBin(dEta,dPhi,fPVz)); }
      else { dWeight = h2Weights->GetBinContent(h2Weights->FindFixBin(dEta,dPhi)); }
      if(dWeight <= 0.0) dWeight = 1.0;
    }

    if(!bHasGap) // no eta gap
    {
      for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
    }
    else
    {
      // RFP in positive eta acceptance
      if(dEta > dEtaLimit)
      {
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }
      // RFP in negative eta acceptance
      if(dEta < -dEtaLimit)
      {
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }

      // RFP in middle (for 3sub) if gap > 0
      if(bHas3sub && (TMath::Abs(dEta) < dEtaLimit) )
      {
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecQmid[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
      }

    } // endif {dEtaGap}
  } // endfor {tracks} particle loop

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::FillPOIsVectors(const Double_t dEtaGap, const PartSpecies species, const Double_t dPtLow, const Double_t dPtHigh, const Double_t dMassLow, const Double_t dMassHigh)
{
  // Filling p,q and s flow vectors with POIs (given by species) for differential flow calculation
  // *************************************************************
  Double_t dEtaLimit = dEtaGap / 2.0;
  Bool_t bHasGap = kFALSE; if(dEtaGap > -1.0) { bHasGap = kTRUE; }

  std::vector<AliVTrack*>* vector = fVector[species];
  TH2D* h2Weights = fh2Weights[species];
  TH3D* h3Weights = fh3Weights[species];
  Bool_t bHasMass = kFALSE; if(species == kK0s || species == kLambda || species == kPhi) { bHasMass = kTRUE; }

  if(!vector) { AliError("Vector with selected POIs not found."); return; }
  if(fFlowUseWeights)
  {
    if(!fFlowUse3Dweights && !h2Weights) { AliError("Histogram with weights not found."); return; }
    if(fFlowUse3Dweights && !h3Weights) { AliError("Histogram with weights not found."); return; }
  }

  // clearing output (global) flow vectors
  ResetFlowVector(fFlowVecPpos);
  if(bHasGap) { ResetFlowVector(fFlowVecPneg); } else { ResetFlowVector(fFlowVecS); }

  for(auto part = vector->begin(); part != vector->end(); ++part)
  {
    Double_t dPt = (*part)->Pt();
    Double_t dPhi = (*part)->Phi();
    Double_t dEta = (*part)->Eta();
    Double_t dMass = 0.0;

    // checking if pt is within pt (bin) range
    if(dPt < dPtLow || dPt >= dPtHigh) { continue; }

    // checking if mass is within mass (bin) range
    if(bHasMass) { dMass = (*part)->M(); if(dMass < dMassLow || dMass >= dMassHigh) { continue; } }

    // at this point particles corresponding to this pt (& mass) bin survive

    Double_t dWeight = 1.0;
    if(fFlowUseWeights)
    {
      if(fFlowUse3Dweights) { dWeight = h3Weights->GetBinContent(h3Weights->FindFixBin(dEta,dPhi,fPVz)); }
      else { dWeight = h2Weights->GetBinContent(h2Weights->FindFixBin(dEta,dPhi)); }
      if(dWeight <= 0.0) dWeight = 1.0;
    }

    if(!bHasGap) // no eta gap
    {
      for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
        for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
        {
          Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
          Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
          fFlowVecPpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);

          // check if track (passing criteria) is overlapping with RFPs pT region; if so, fill S (q) vector
          // in case of charged, pions, kaons or protons (one witout mass)
          if(!bHasMass && IsWithinRefs(static_cast<const AliAODTrack*>(*part)))
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecS[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
        }

    }
    else // with eta gap
    {
      if(dEta > dEtaLimit) // particle in positive eta acceptance
      {
        for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
          for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
          {
            Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
            Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
            fFlowVecPpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
          }
       }
       if(dEta < -dEtaLimit) // particle in negative eta acceptance
       {
         for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
           for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
           {
             Double_t dCos = TMath::Power(dWeight,iPower) * TMath::Cos(iHarm * dPhi);
             Double_t dSin = TMath::Power(dWeight,iPower) * TMath::Sin(iHarm * dPhi);
             fFlowVecPneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
           }
       }
     } // endif {dEtaGap}
   } // endfor {tracks}
   return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
      array[iHarm][iPower] = TComplex(0,0,kFALSE);
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::ListFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax])
{
  // List all values of given flow vector TComplex array
  // *************************************************************
  printf(" ### Listing (TComplex) flow vector array ###########################\n");
  for(Short_t iHarm(0); iHarm < fFlowNumHarmonicsMax; iHarm++)
  {
    printf("Harm %d (power):",iHarm);
    for(Short_t iPower(0); iPower < fFlowNumWeightPowersMax; iPower++)
    {
        printf("|(%d) %g+%g(i)",iPower, array[iHarm][iPower].Re(), array[iHarm][iPower].Im());
    }
    printf("\n");
  }
  return;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetSamplingIndex()
{
  // Assessing sampling index based on generated random number
  // returns centrality index
  // *************************************************************

  Short_t index = 0;

  if(fSampling && fNumSamples > 1)
  {
    TRandom3 rr(0);
    Double_t ranNum = rr.Rndm(); // getting random number in (0,1)
    Double_t generated = ranNum * fNumSamples; // getting random number in range (0, fNumSamples)

    // finding right index for sampling based on generated number and total number of samples
    for(Short_t i(0); i < fNumSamples; i++)
    {
      if(generated < (i+1) )
      {
        index = i;
        break;
      }
    }
  }

  return index;
}
//_____________________________________________________________________________
Short_t AliAnalysisTaskUniFlow::GetCentralityIndex()
{
  // Estimating centrality percentile based on selected estimator.
  // (Default) If no multiplicity estimator is specified, percentile is estimated as number of selected / filtered charged tracks (NRFP).
  // If a valid multiplicity estimator is specified, centrality percentile is estimated via AliMultSelection
  // otherwise -1 is returned (and event is skipped)
  // *************************************************************

  Short_t iCentralityIndex = -1;

  // assigning centrality based on number of selected charged tracks
  if(fMultEstimator == kRFP)
  { iCentralityIndex = fVector[kRefs]->size(); }
  else
  {
    // checking AliMultSelection
    AliMultSelection* multSelection = (AliMultSelection*) fEventAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }

    Float_t dPercentile = multSelection->GetMultiplicityPercentile(GetMultiEstimatorLabel(fMultEstimator));
    if(dPercentile > 100 || dPercentile < 0)
    { AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); return -1; }

    iCentralityIndex = (Short_t) dPercentile;
  }

  // transfering centrality percentile to multiplicity bin (if fixed size bins are used)
  if(fUseFixedMultBins)
  {
    for(Int_t multIndex(0); multIndex < fNumMultBins; multIndex++)
    {
      if(iCentralityIndex >= fMultBins[multIndex] && iCentralityIndex < fMultBins[multIndex+1]) { iCentralityIndex = multIndex; break; }
    }
  }

  return iCentralityIndex;
}
//_____________________________________________________________________________
const char* AliAnalysisTaskUniFlow::GetMultiEstimatorLabel(MultiEst est)
{
  // Return string with estimator name or 'n/a' if not available
  // *************************************************************

  switch (est)
  {
    case kRFP: return "N_{RFP}";
    case kV0A: return "V0A";
    case kV0C: return "V0C";
    case kV0M: return "V0M";
    case kCL0: return "CL1";
    case kCL1: return "CL2";
    case kZNA: return "ZNA";
    case kZNC: return "ZNC";
    default: return "n/a";
  }
  return "n/a";
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::HasTrackPIDTPC(const AliAODTrack* track)
{
  // Checks if the track has ok PID information from TPC
  // *************************************************************
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
  return (pidStatusTPC == AliPIDResponse::kDetPidOk);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskUniFlow::HasTrackPIDTOF(const AliAODTrack* track)
{
  // Checks if the track has ok PID information from TOF
  // *************************************************************
  if(!track || !fPIDResponse) return kFALSE;
  AliPIDResponse::EDetPidStatus pidStatusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  return ((pidStatusTOF == AliPIDResponse::kDetPidOk) && (track->GetStatus()& AliVTrack::kTOFout) && (track->GetStatus()& AliVTrack::kTIME));
}
//_____________________________________________________________________________
void AliAnalysisTaskUniFlow::Terminate(Option_t* option)
{
  // called on end of task, after all events are processed
  // *************************************************************
  AliAnalysisTaskSE::Terminate(option);
  return;
}
//_____________________________________________________________________________
// Set of methods returning given complex flow vector based on flow harmonics (n) and weight power indexes (p)
// a la General Framework implementation.
// Q: flow vector of RFPs (with/out eta gap)
// P: flow vector of POIs (with/out eta gap) (in usual notation p)
// S: flow vector of overlaping RFPs and POIs (in usual notation q)

TComplex AliAnalysisTaskUniFlow::Q(const Short_t n, const Short_t p)
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::QGapPos(const Short_t n, const Short_t p)
{
  if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
  else return fFlowVecQpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::QGapNeg(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecQneg[-n][p]);
  else return fFlowVecQneg[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::QGapMid(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecQmid[-n][p]);
  else return fFlowVecQmid[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::P(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::PGapPos(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::PGapNeg(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPneg[-n][p]);
  else return fFlowVecPneg[n][p];
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::S(const Short_t n, const Short_t p)
{
  if(n < 0) return TComplex::Conjugate(fFlowVecS[-n][p]);
  else return fFlowVecS[n][p];
}
//____________________________________________________________________

// Set of flow calculation methods for cumulants of different orders with/out eta gap

TComplex AliAnalysisTaskUniFlow::Two(const Short_t n1, const Short_t n2)
{
  TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoGap(const Short_t n1, const Short_t n2)
{
  TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiff(const Short_t n1, const Short_t n2)
{
  TComplex formula = P(n1,1)*Q(n2,1) - S(n1+n2,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiffGapPos(const Short_t n1, const Short_t n2)
{
  TComplex formula = PGapPos(n1,1)*QGapNeg(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::TwoDiffGapNeg(const Short_t n1, const Short_t n2)
{
  TComplex formula = PGapNeg(n1,1)*QGapPos(n2,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::Three(const Short_t n1, const Short_t n2, const Short_t n3)
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
 		                 - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::Four(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*Q(n1+n3+n4,3)+2.0*Q(n1,1)*Q(n2+n3+n4,3)-6.0*Q(n1+n2+n3+n4,4);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::FourGap(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)-QGapPos(n1+n2,2)*QGapNeg(n3,1)*QGapNeg(n4,1)
                    -QGapPos(n1,1)*QGapPos(n2,1)*QGapNeg(n3+n4,2)+QGapPos(n1+n2,2)*QGapNeg(n3+n4,2);
	return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::Four3sub(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  // left = neg, middle = mid; rigth = pos
  TComplex formula = QGapMid(n1,1)*QGapMid(n2,1)*QGapNeg(n3,1)*QGapPos(n4,1)-QGapMid(n1+n2,2)*QGapNeg(n3,1)*QGapPos(n4,1);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::ThreeDiff(const Short_t n1, const Short_t n2, const Short_t n3)
{
  TComplex formula = P(n1,1)*Q(n2,1)*Q(n3,1)-S(n1+n2,2)*Q(n3,1)-S(n1+n3,2)*Q(n2,1)
 		                 - P(n1,1)*Q(n2+n3,2)+2.0*S(n1+n2+n3,3);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::ThreeDiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3)
{
  TComplex formula = PGapPos(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1) - PGapPos(n1,1)*QGapNeg(n2+n3,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::ThreeDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3)
{
  TComplex formula = PGapNeg(n1,1)*QGapPos(n2,1)*QGapPos(n3,1) - PGapNeg(n1,1)*QGapPos(n2+n3,2);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::FourDiff(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = P(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-S(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*S(n1+n3,2)*Q(n4,1)
                    - P(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*S(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2)
                    + Q(n2+n3,2)*S(n1+n4,2)-P(n1,1)*Q(n3,1)*Q(n2+n4,2)+S(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*S(n1+n2+n4,3)-P(n1,1)*Q(n2,1)*Q(n3+n4,2)+S(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*S(n1+n3+n4,3)+2.0*P(n1,1)*Q(n2+n3+n4,3)-6.0*S(n1+n2+n3+n4,4);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::FourDiffGapPos(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = PGapPos(n1,1)*QGapNeg(n2,1)*QGapNeg(n3,1)*QGapNeg(n4,1)
                     - PGapPos(n1,1)*QGapNeg(n2+n3,2)*QGapNeg(n4,1) - PGapPos(n1,1)*QGapNeg(n3,1)*QGapNeg(n2+n4,2)
                     - PGapPos(n1,1)*QGapNeg(n2,1)*QGapNeg(n3+n4,2) + 2.0*PGapPos(n1,1)*QGapNeg(n2+n3+n4,3);
  return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskUniFlow::FourDiffGapNeg(const Short_t n1, const Short_t n2, const Short_t n3, const Short_t n4)
{
  TComplex formula = PGapNeg(n1,1)*QGapPos(n2,1)*QGapPos(n3,1)*QGapPos(n4,1)
                     - PGapNeg(n1,1)*QGapPos(n2+n3,2)*QGapPos(n4,1) - PGapNeg(n1,1)*QGapPos(n3,1)*QGapPos(n2+n4,2)
                     - PGapNeg(n1,1)*QGapPos(n2,1)*QGapPos(n3+n4,2) + 2.0*PGapNeg(n1,1)*QGapPos(n2+n3+n4,3);
  return formula;
}
//____________________________________________________________________
void AliAnalysisTaskUniFlow::UserCreateOutputObjects()
{
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // *************************************************************

  // task initialization
  fInit = InitializeTask();
  if(!fInit) return;

  // list all parameters used in this analysis
  ListParameters();


  // creating output lists
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    fListFlow[iSpec] = new TList();
    fListFlow[iSpec]->SetOwner(kTRUE);
    fListFlow[iSpec]->SetName(Form("fFlow%s",GetSpeciesName(PartSpecies(iSpec))));
  }
  fFlowWeights = new TList();
  fFlowWeights->SetOwner(kTRUE);
  fFlowWeights->SetName("fFlowWeights");

  fQAEvents = new TList();
  fQAEvents->SetOwner(kTRUE);
  fQACharged = new TList();
  fQACharged->SetOwner(kTRUE);
  fQAPID = new TList();
  fQAPID->SetOwner(kTRUE);
  fQAPhi = new TList();
  fQAPhi->SetOwner(kTRUE);
  fQAV0s = new TList();
  fQAV0s->SetOwner(kTRUE);

  fEventCuts.AddQAplotsToList(fQAEvents);

  // setting number of bins based on set range with fixed width
  const Int_t iPOIsPtNumBins = (Int_t) (fFlowPOIsPtMax - fFlowPOIsPtMin) / 0.1 + 0.5; // fixed width 0.1 GeV/c
  const Int_t iMultNumBins = fFlowCentMax - fFlowCentMin; // fixed unit percentile or Nch bin width

  // creating output correlations profiles based on FlowTasks
  Int_t iNumTasks = fVecFlowTask.size();
  for(Int_t iTask(0); iTask < iNumTasks; ++iTask)
  {
    FlowTask* task = fVecFlowTask.at(iTask);
    if(!task) { fInit = kFALSE; AliError(Form("FlowTask %d does not exists\n",iTask)); return; }

    Bool_t bHasGap = task->HasGap();
    const char* corName = task->fsName.Data();
    const char* corLabel = task->fsLabel.Data();

    for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
    {
      // check if FlowTask should be done for all flow particles (RFP/POI/Both)
      if(task->fPart == FlowTask::kRFP && iSpec != kRefs) { continue; }
      if(task->fPart == FlowTask::kPOI && iSpec == kRefs) { continue; }

      if(!fProcessSpec[iSpec]) { continue; }

      for(Int_t iSample(0); iSample < fNumSamples; ++iSample)
      {
        if(!fSampling && iSample > 0) { break; }
        if(iSample > 0 && (iSpec == kK0s || iSpec == kLambda || iSpec == kPhi)) { break; } // reconstructed are not sampled

        TH1* profile = 0x0;
        TH1* profileNeg = 0x0;

        switch(iSpec)
        {
          case kRefs :
          {
            profile = new TProfile(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s; %s",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax);
            break;
          }

          case kCharged :
          case kPion :
          case kKaon :
          case kProton :
          {
            profile = new TProfile2D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
            if(bHasGap) { profileNeg = new TProfile2D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c})",GetSpeciesLabel(PartSpecies(iSpec)), corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax); }
            break;
          }

          case kK0s:
          {
            profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
            if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax); }
            break;
          }

          case kLambda:
          {
            profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
            if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fV0sNumBinsMass,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax); }
            break;
          }

          case kPhi:
          {
            profile = new TProfile3D(Form("%s_Pos_sample%d",corName,iSample), Form("%s: %s (Pos); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax);
            if(bHasGap) { profileNeg = new TProfile3D(Form("%s_Neg_sample%d",corName,iSample), Form("%s: %s (Neg); %s; #it{p}_{T} (GeV/#it{c}); #it{m}_{inv} (GeV/#it{c}^{2})",GetSpeciesLabel(PartSpecies(iSpec)),corLabel,GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax, fPhiNumBinsMass,fCutPhiInvMassMin,fCutPhiInvMassMax); }
            break;
          }
        }

        if(!profile) { fInit = kFALSE; AliError("Profile (Pos) NOT created!"); task->Print(); return; }

        // check if same profile does not exists already
        if(fListFlow[iSpec]->FindObject(profile->GetName())) {
          AliError(Form("FlowTask %d : Profile '%s' already exists! Please check run macro for FlowTask duplicates!",iTask,profile->GetName()));
          fInit = kFALSE;
          task->Print();
          delete profile;
          return;
        }

        profile->Sumw2();
        fListFlow[iSpec]->Add(profile);

        if(bHasGap && iSpec != kRefs) { // Refs does not distinquish Pos/Neg
          if(!profileNeg) { fInit = kFALSE; AliError("Profile (Neg) NOT created!"); task->Print(); return; }

          // same for Neg
          if(fListFlow[iSpec]->FindObject(profileNeg->GetName())) {
            AliError(Form("FlowTask %d : Profile '%s' already exists! Please check run macro for FlowTask duplicates!",iTask,profile->GetName()));
            fInit = kFALSE;
            task->Print();
            delete profileNeg;
            return;
          }

          profileNeg->Sumw2();
          fListFlow[iSpec]->Add(profileNeg);
        }

      } // end-for {iSample}
    } // end-for {iSpec}
  } // end-for {iTask}

  // creating GF weights
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec)
  {
    if(!fProcessSpec[iSpec]) { continue; }

    if(fFlowFillWeights)
    {
      const char* weightName = Form("fh3Weights%s",GetSpeciesName(PartSpecies(iSpec)));
      const char* weightLabel = Form("Weights: %s; #varphi; #eta; PV-z (cm)", GetSpeciesName(PartSpecies(iSpec)));
      fh3Weights[iSpec] = new TH3D(weightName, weightLabel, 100,0.0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3Weights[iSpec]->Sumw2();
      fFlowWeights->Add(fh3Weights[iSpec]);
    }

    if(fFlowUseWeights)
    {
      const char* weightName = Form("fh3AfterWeights%s",GetSpeciesName(PartSpecies(iSpec)));
      const char* weightLabel = Form("Weights (after): %s; #varphi; #eta; PV-z (cm)",GetSpeciesLabel(PartSpecies(iSpec)));
      fh3AfterWeights[iSpec] = new TH3D(weightName,weightLabel, 100,0,TMath::TwoPi(), 80,-1.0,1.0, 2*fPVtxCutZ,-fPVtxCutZ,fPVtxCutZ);
      fh3AfterWeights[iSpec]->Sumw2();
      fFlowWeights->Add(fh3AfterWeights[iSpec]);
    }
  } // end-for {iSpec}


  // NB: bellow NOT re-factorized

  // creating QA histos
    // event histogram
    fhEventSampling = new TH2D("fhEventSampling",Form("Event sampling; %s; sample index", GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, fNumSamples,0,fNumSamples);
    fQAEvents->Add(fhEventSampling);
    fhEventCentrality = new TH1D("fhEventCentrality",Form("Event centrality (%s); %s", GetMultiEstimatorLabel(fMultEstimator), GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax);
    fQAEvents->Add(fhEventCentrality);
    fh2EventCentralityNumRefs = new TH2D("fh2EventCentralityNumRefs",Form("Event centrality (%s) vs. N_{RFP}; %s; N_{RFP}",GetMultiEstimatorLabel(fMultEstimator), GetMultiEstimatorLabel(fMultEstimator)), iMultNumBins,fFlowCentMin,fFlowCentMax, 150,0,150);
    fQAEvents->Add(fh2EventCentralityNumRefs);

    TString sEventCounterLabel[] = {"Input","Physics selection OK","EventCuts OK","Event OK","#RPFs OK","Multiplicity OK","Selected"};
    const Short_t iEventCounterBins = sizeof(sEventCounterLabel)/sizeof(sEventCounterLabel[0]);
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    for(Short_t i(0); i < iEventCounterBins; i++) fhEventCounter->GetXaxis()->SetBinLabel(i+1, sEventCounterLabel[i].Data() );
    fQAEvents->Add(fhEventCounter);

    if(fEventRejectAddPileUp)
    {
      fhQAEventsfMult32vsCentr = new TH2D("fhQAEventsfMult32vsCentr", "; centrality V0M; TPC multiplicity (FB32)", 100, 0, 100, 1000, 0, 4000);
      fQAEvents->Add(fhQAEventsfMult32vsCentr);
      fhQAEventsMult128vsCentr = new TH2D("fhQAEventsfMult128vsCentr", "; centrality V0M; TPC multiplicity (FB128)", 100, 0, 100, 1000, 0, 5000);
      fQAEvents->Add(fhQAEventsMult128vsCentr);
      fhQAEventsfMultTPCvsTOF = new TH2D("fhQAEventsfMultTPCvsTOF", "; TPC FB32 multiplicity; TOF multiplicity", 400, 0, 4000, 200, 0, 2000);
      fQAEvents->Add(fhQAEventsfMultTPCvsTOF);
      fhQAEventsfMultTPCvsESD = new TH2D("fhQAEventsfMultTPCvsESD", "; TPC FB128 multiplicity; ESD multiplicity", 700, 0, 7000, 1000, -1000, 35000);
      fQAEvents->Add(fhQAEventsfMultTPCvsESD);
    }

    // Making THnSparse distribution of candidates
    // species independent
    TString sLabelCand[SparseCand::kDim];
    sLabelCand[SparseCand::kInvMass] = "#it{m}_{inv} (GeV/#it{c}^{2})";
    sLabelCand[SparseCand::kCent] = GetMultiEstimatorLabel(fMultEstimator);
    sLabelCand[SparseCand::kPt] = "#it{p}_{T} (GeV/c)";
    sLabelCand[SparseCand::kEta] = "#eta";
    TString sAxes = TString(); for(Int_t i(0); i < SparseCand::kDim; ++i) { sAxes += Form("%s; ",sLabelCand[i].Data()); }

    Int_t iNumBinsCand[SparseCand::kDim]; Double_t dMinCand[SparseCand::kDim]; Double_t dMaxCand[SparseCand::kDim];
    iNumBinsCand[SparseCand::kCent] = iMultNumBins; dMinCand[SparseCand::kCent] = fFlowCentMin; dMaxCand[SparseCand::kCent] = fFlowCentMax;
    iNumBinsCand[SparseCand::kPt] = iPOIsPtNumBins; dMinCand[SparseCand::kPt] = fFlowPOIsPtMin; dMaxCand[SparseCand::kPt] = fFlowPOIsPtMax;
    iNumBinsCand[SparseCand::kEta] = 2*fCutV0sMotherEtaMax/0.05; dMinCand[SparseCand::kEta] = -fCutV0sMotherEtaMax; dMaxCand[SparseCand::kEta] = fCutV0sMotherEtaMax;

    // species dependent
    if(fProcessSpec[kK0s] || fProcessSpec[kLambda])
    {
      iNumBinsCand[SparseCand::kInvMass] = fV0sNumBinsMass; dMinCand[SparseCand::kInvMass] = fCutV0sInvMassK0sMin; dMaxCand[SparseCand::kInvMass] = fCutV0sInvMassK0sMax;
      fhsV0sCandK0s = new THnSparseD("fhsV0sCandK0s",Form("K_{S}^{0}: Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
      fhsV0sCandK0s->Sumw2();
      fListFlow[kK0s]->Add(fhsV0sCandK0s);

      iNumBinsCand[SparseCand::kInvMass] = fV0sNumBinsMass; dMinCand[SparseCand::kInvMass] = fCutV0sInvMassLambdaMin; dMaxCand[SparseCand::kInvMass] = fCutV0sInvMassLambdaMax;
      fhsV0sCandLambda = new THnSparseD("fhsV0sCandLambda",Form("#Lambda: Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
      fhsV0sCandLambda->Sumw2();
      fListFlow[kLambda]->Add(fhsV0sCandLambda);
    }

    if(fProcessSpec[kPhi])
    {
      iNumBinsCand[SparseCand::kEta] = 2*fCutPhiMotherEtaMax/0.05; dMinCand[SparseCand::kEta] = -fCutPhiMotherEtaMax; dMaxCand[SparseCand::kEta] = fCutPhiMotherEtaMax;
      iNumBinsCand[SparseCand::kInvMass] = fPhiNumBinsMass; dMinCand[SparseCand::kInvMass] = fCutPhiInvMassMin; dMaxCand[SparseCand::kInvMass] = fCutPhiInvMassMax;

      fhsPhiCandSig = new THnSparseD("fhsPhiCandSig",Form("#phi (Sig): Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
      fhsPhiCandSig->Sumw2();
      fListFlow[kPhi]->Add(fhsPhiCandSig);

      fhsPhiCandBg = new THnSparseD("fhsPhiCandBg",Form("#phi (Bg): Distribution; %s;", sAxes.Data()), SparseCand::kDim, iNumBinsCand, dMinCand, dMaxCand);
      fhsPhiCandBg->Sumw2();
      fListFlow[kPhi]->Add(fhsPhiCandBg);
    }

    // charged (tracks) histograms
    fhRefsMult = new TH1D("fhRefsMult","RFPs: Multiplicity; multiplicity", 1000,0,1000);
    fQACharged->Add(fhRefsMult);
    fhRefsPt = new TH1D("fhRefsPt","RFPs: #it{p}_{T};  #it{p}_{T} (GeV/#it{c})", 300,0,30);
    fQACharged->Add(fhRefsPt);
    fhRefsEta = new TH1D("fhRefsEta","RFPs: #eta; #eta", 151,-1.5,1.5);
    fQACharged->Add(fhRefsEta);
    fhRefsPhi = new TH1D("fhRefsPhi","RFPs: #varphi; #varphi", 100,0,TMath::TwoPi());
    fQACharged->Add(fhRefsPhi);
    fpRefsMult = new TProfile("fpRefsMult","Ref mult; %s", iMultNumBins,fFlowCentMin,fFlowCentMax);
    fpRefsMult->Sumw2();
    fQACharged->Add(fpRefsMult);

    TString sChargedCounterLabel[] = {"Input","Pt","Eta","FB","#TPC-Cls","DCA-z","DCA-xy","Selected"};
    const Short_t iNBinsChargedCounter = sizeof(sChargedCounterLabel)/sizeof(sChargedCounterLabel[0]);
    fhChargedCounter = new TH1D("fhChargedCounter","Charged tracks: Counter",iNBinsChargedCounter,0,iNBinsChargedCounter);
    for(Short_t i(0); i < iNBinsChargedCounter; i++) fhChargedCounter->GetXaxis()->SetBinLabel(i+1, sChargedCounterLabel[i].Data() );
    fQACharged->Add(fhChargedCounter);

    // PID tracks histograms
    if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton] || fProcessSpec[kPhi])
    {
      fhPIDPionMult = new TH1D("fhPIDPionMult","PID: #pi: Multiplicity; multiplicity", 200,0,200);
      fQAPID->Add(fhPIDPionMult);
      fhPIDKaonMult = new TH1D("fhPIDKaonMult","PID: K: Multiplicity; multiplicity", 100,0,100);
      fQAPID->Add(fhPIDKaonMult);
      fhPIDProtonMult = new TH1D("fhPIDProtonMult","PID: p: Multiplicity; multiplicity", 100,0,100);
      fQAPID->Add(fhPIDProtonMult);

      if(fMC)
      {
        fhMCRecoSelectedPionPt = new TH1D("fhMCRecoSelectedPionPt","fhMCRecoSelectedPionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedPionPt);
        fhMCRecoSelectedTruePionPt = new TH1D("fhMCRecoSelectedTruePionPt","fhMCRecoSelectedTruePionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedTruePionPt);
        fhMCRecoAllPionPt = new TH1D("fhMCRecoAllPionPt","fhMCRecoAllPionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoAllPionPt);
        fhMCGenAllPionPt = new TH1D("fhMCGenAllPionPt","fhMCGenAllPionPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCGenAllPionPt);

        fhMCRecoSelectedKaonPt = new TH1D("fhMCRecoSelectedKaonPt","fhMCRecoSelectedKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedKaonPt);
        fhMCRecoSelectedTrueKaonPt = new TH1D("fhMCRecoSelectedTrueKaonPt","fhMCRecoSelectedTrueKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedTrueKaonPt);
        fhMCRecoAllKaonPt = new TH1D("fhMCRecoAllKaonPt","fhMCRecoAllKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoAllKaonPt);
        fhMCGenAllKaonPt = new TH1D("fhMCGenAllKaonPt","fhMCGenAllKaonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCGenAllKaonPt);

        fhMCRecoSelectedProtonPt = new TH1D("fhMCRecoSelectedProtonPt","fhMCRecoSelectedProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedProtonPt);
        fhMCRecoSelectedTrueProtonPt = new TH1D("fhMCRecoSelectedTrueProtonPt","fhMCRecoSelectedTrueProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoSelectedTrueProtonPt);
        fhMCRecoAllProtonPt = new TH1D("fhMCRecoAllProtonPt","fhMCRecoAllProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCRecoAllProtonPt);
        fhMCGenAllProtonPt = new TH1D("fhMCGenAllProtonPt","fhMCGenAllProtonPt; p_{T} (GeV/c); Counts", iPOIsPtNumBins,fFlowPOIsPtMin,fFlowPOIsPtMax);
        fQAPID->Add(fhMCGenAllProtonPt);
      }


      if(fFillQA)
      {
        fhPIDPionPt = new TH1D("fhPIDPionPt","PID: #pi: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDPionPt);
        fhPIDPionPhi = new TH1D("fhPIDPionPhi","PID: #pi: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDPionPhi);
        fhPIDPionEta = new TH1D("fhPIDPionEta","PID: #pi: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDPionEta);
        fhPIDPionCharge = new TH1D("fhPIDPionCharge","PID: #pi: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDPionCharge);
        fhPIDKaonPt = new TH1D("fhPIDKaonPt","PID: K: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDKaonPt);
        fhPIDKaonPhi = new TH1D("fhPIDKaonPhi","PID: K: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDKaonPhi);
        fhPIDKaonEta = new TH1D("fhPIDKaonEta","PID: K: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDKaonEta);
        fhPIDKaonCharge = new TH1D("fhPIDKaonCharge","PID: K: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDKaonCharge);
        fhPIDProtonPt = new TH1D("fhPIDProtonPt","PID: p: #it{p}_{T}; #it{p}_{T}", 150,0.,30.);
        fQAPID->Add(fhPIDProtonPt);
        fhPIDProtonPhi = new TH1D("fhPIDProtonPhi","PID: p: #varphi; #varphi", 100,0,TMath::TwoPi());
        fQAPID->Add(fhPIDProtonPhi);
        fhPIDProtonEta = new TH1D("fhPIDProtonEta","PID: p: #eta; #eta", 151,-1.5,1.5);
        fQAPID->Add(fhPIDProtonEta);
        fhPIDProtonCharge = new TH1D("fhPIDProtonCharge","PID: p: charge; charge", 3,-1.5,1.5);
        fQAPID->Add(fhPIDProtonCharge);
        fh2PIDPionTPCdEdx = new TH2D("fh2PIDPionTPCdEdx","PID: #pi: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDPionTPCdEdx);
        fh2PIDPionTOFbeta = new TH2D("fh2PIDPionTOFbeta","PID: #pi: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDPionTOFbeta);
        fh2PIDKaonTPCdEdx = new TH2D("fh2PIDKaonTPCdEdx","PID: K: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDKaonTPCdEdx);
        fh2PIDKaonTOFbeta = new TH2D("fh2PIDKaonTOFbeta","PID: K: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDKaonTOFbeta);
        fh2PIDProtonTPCdEdx = new TH2D("fh2PIDProtonTPCdEdx","PID: p: TPC dE/dx; #it{p} (GeV/#it{c}); TPC dE/dx", 200,0,20, 131,-10,1000);
        fQAPID->Add(fh2PIDProtonTPCdEdx);
        fh2PIDProtonTOFbeta = new TH2D("fh2PIDProtonTOFbeta","PID: p: TOF #beta; #it{p} (GeV/#it{c});TOF #beta", 200,0,20, 101,-0.1,1.5);
        fQAPID->Add(fh2PIDProtonTOFbeta);
        fh2PIDPionTPCnSigmaPion = new TH2D("fh2PIDPionTPCnSigmaPion","PID: #pi: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaPion);
        fh2PIDPionTOFnSigmaPion = new TH2D("fh2PIDPionTOFnSigmaPion","PID: #pi: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaPion);
        fh2PIDPionTPCnSigmaKaon = new TH2D("fh2PIDPionTPCnSigmaKaon","PID: #pi: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaKaon);
        fh2PIDPionTOFnSigmaKaon = new TH2D("fh2PIDPionTOFnSigmaKaon","PID: #pi: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaKaon);
        fh2PIDPionTPCnSigmaProton = new TH2D("fh2PIDPionTPCnSigmaProton","PID: #pi: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTPCnSigmaProton);
        fh2PIDPionTOFnSigmaProton = new TH2D("fh2PIDPionTOFnSigmaProton","PID: #pi: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDPionTOFnSigmaProton);
        fh2PIDPionBayesPion = new TH2D("fh2PIDPionBayesPion","PID: #pi: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesPion);
        fh2PIDPionBayesKaon = new TH2D("fh2PIDPionBayesKaon","PID: #pi: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesKaon);
        fh2PIDPionBayesProton = new TH2D("fh2PIDPionBayesProton","PID: #pi: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDPionBayesProton);
        fh2PIDKaonTPCnSigmaPion = new TH2D("fh2PIDKaonTPCnSigmaPion","PID: K: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaPion);
        fh2PIDKaonTOFnSigmaPion = new TH2D("fh2PIDKaonTOFnSigmaPion","PID: K: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaPion);
        fh2PIDKaonTPCnSigmaKaon = new TH2D("fh2PIDKaonTPCnSigmaKaon","PID: K: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaKaon);
        fh2PIDKaonTOFnSigmaKaon = new TH2D("fh2PIDKaonTOFnSigmaKaon","PID: K: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaKaon);
        fh2PIDKaonTPCnSigmaProton = new TH2D("fh2PIDKaonTPCnSigmaProton","PID: K: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTPCnSigmaProton);
        fh2PIDKaonTOFnSigmaProton = new TH2D("fh2PIDKaonTOFnSigmaProton","PID: K: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDKaonTOFnSigmaProton);
        fh2PIDKaonBayesPion = new TH2D("fh2PIDKaonBayesPion","PID: K: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesPion);
        fh2PIDKaonBayesKaon = new TH2D("fh2PIDKaonBayesKaon","PID: K: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesKaon);
        fh2PIDKaonBayesProton = new TH2D("fh2PIDKaonBayesProton","PID: K: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDKaonBayesProton);
        fh2PIDProtonTPCnSigmaPion = new TH2D("fh2PIDProtonTPCnSigmaPion","PID: p: TPC n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaPion);
        fh2PIDProtonTOFnSigmaPion = new TH2D("fh2PIDProtonTOFnSigmaPion","PID: p: TOF n#sigma (#pi hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaPion);
        fh2PIDProtonTPCnSigmaKaon = new TH2D("fh2PIDProtonTPCnSigmaKaon","PID: p: TPC n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaKaon);
        fh2PIDProtonTOFnSigmaKaon = new TH2D("fh2PIDProtonTOFnSigmaKaon","PID: p: TOF n#sigma (K hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaKaon);
        fh2PIDProtonTPCnSigmaProton = new TH2D("fh2PIDProtonTPCnSigmaProton","PID: p: TPC n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TPC n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTPCnSigmaProton);
        fh2PIDProtonTOFnSigmaProton = new TH2D("fh2PIDProtonTOFnSigmaProton","PID: p: TOF n#sigma (p hyp.); #it{p}_{T} (GeV/#it{c}); TOF n#sigma", 200,0,20, 21,-11,10);
        fQAPID->Add(fh2PIDProtonTOFnSigmaProton);
        fh2PIDProtonBayesPion = new TH2D("fh2PIDProtonBayesPion","PID: p: Bayes probability (#pi hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesPion);
        fh2PIDProtonBayesKaon = new TH2D("fh2PIDProtonBayesKaon","PID: p: Bayes probability (K hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesKaon);
        fh2PIDProtonBayesProton = new TH2D("fh2PIDProtonBayesProton","PID: p: Bayes probability (p hyp.); #it{p}_{T} (GeV/#it{c}); Bayes prob.", 200,0,20, 50,0,1);
        fQAPID->Add(fh2PIDProtonBayesProton);

      }

    } //endif {fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton]}

    if(fProcessSpec[kPhi])
    {
      TString sPhiCounterLabel[] = {"Input","InvMass","Pt","Eta","Selected","Unlike-sign","BG"};
      const Short_t iNBinsPhiCounter = sizeof(sPhiCounterLabel)/sizeof(sPhiCounterLabel[0]);
      fhPhiCounter = new TH1D("fhPhiCounter","#phi: Counter",iNBinsPhiCounter,0,iNBinsPhiCounter);
      for(Short_t i(0); i < iNBinsPhiCounter; i++) fhPhiCounter->GetXaxis()->SetBinLabel(i+1, sPhiCounterLabel[i].Data() );
      fQAPhi->Add(fhPhiCounter);

      fhPhiMult = new TH1D("fhPhiMult","#phi: Multiplicity; Multiplicity", 150,0,150);
      fQAPhi->Add(fhPhiMult);
      fhPhiBGMult = new TH1D("fhPhiBGMult","#phi (BG): Multiplicity; Multiplicity", 150,0,150);
      fQAPhi->Add(fhPhiBGMult);
      fhPhiInvMass = new TH1D("fhPhiInvMass","#phi: InvMass; #it{m}_{inv} (GeV/#it{c}^{2})", 90,fCutPhiInvMassMin,fCutPhiInvMassMax);
      fQAPhi->Add(fhPhiInvMass);
      fhPhiBGInvMass = new TH1D("fhPhiBGInvMass","#phi (BG): InvMass; #it{m}_{inv} (GeV/#it{c}^{2})", 90,fCutPhiInvMassMin,fCutPhiInvMassMax);
      fQAPhi->Add(fhPhiBGInvMass);
      fhPhiCharge = new TH1D("fhPhiCharge","#phi: charge; charge", 5,-2.5,2.5);
      fQAPhi->Add(fhPhiCharge);
      fhPhiBGCharge = new TH1D("fhPhiBGCharge","#phi (BG): charge; charge", 5,-2.5,2.5);
      fQAPhi->Add(fhPhiBGCharge);
      fhPhiPt = new TH1D("fhPhiPt","#phi: #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 150,0.,30.);
      fQAPhi->Add(fhPhiPt);
      fhPhiEta = new TH1D("fhPhiEta","#phi: #eta; #eta", 151,-1.5,1.5);
      fQAPhi->Add(fhPhiEta);
      fhPhiPhi = new TH1D("fhPhiPhi","#phi: #varphi; #varphi", 100,-TMath::Pi(),TMath::Pi());
      fQAPhi->Add(fhPhiPhi);
    } //endif {fProcessPhi}

    // V0 candidates histograms
    if(fProcessSpec[kK0s] || fProcessSpec[kLambda])
    {
      TString sV0sCounterLabel[] = {"Input","Daughters OK","Mother acceptance","Daughter acceptance","Charge","Reconstruction method","TPC refit","Kinks","Daughters track quality","DCA to PV","Daughters DCA","Decay radius","Common passed","K^{0}_{S}","#Lambda/#bar{#Lambda}","K^{0}_{S} && #Lambda/#bar{#Lambda}"};
      const Short_t iNBinsV0sCounter = sizeof(sV0sCounterLabel)/sizeof(sV0sCounterLabel[0]);
      fhV0sCounter = new TH1D("fhV0sCounter","V^{0}: Counter",iNBinsV0sCounter,0,iNBinsV0sCounter);
      for(Short_t i(0); i < iNBinsV0sCounter; i++) fhV0sCounter->GetXaxis()->SetBinLabel(i+1, sV0sCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounter);

      TString sV0sK0sCounterLabel[] = {"Input","#it{y}","InvMass","CPA","Armenteros-Podolanski","c#tau","Daughters PID","Competing InvMass","Selected"};
      const Short_t iNBinsV0sK0sCounter = sizeof(sV0sK0sCounterLabel)/sizeof(sV0sK0sCounterLabel[0]);
      fhV0sCounterK0s = new TH1D("fhV0sCounterK0s","V^{0}: K^{0}_{S} Counter",iNBinsV0sK0sCounter,0,iNBinsV0sK0sCounter);
      for(Short_t i(0); i < iNBinsV0sK0sCounter; i++) fhV0sCounterK0s->GetXaxis()->SetBinLabel(i+1, sV0sK0sCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounterK0s);

      TString sV0sLambdaCounterLabel[] = {"Input","#it{y}","InvMass","CPA","Armenteros-Podolanski","c#tau","Daughter PID","Competing InvMass","Selected","only #Lambda","only #bar{#Lambda}","#Lambda && #bar{#Lambda}"};
      const Short_t iNBinsV0sLambdaCounter = sizeof(sV0sLambdaCounterLabel)/sizeof(sV0sLambdaCounterLabel[0]);
      fhV0sCounterLambda = new TH1D("fhV0sCounterLambda","V^{0}: #Lambda/#bar{#Lambda} Counter",iNBinsV0sLambdaCounter,0,iNBinsV0sLambdaCounter);
      for(Short_t i(0); i < iNBinsV0sLambdaCounter; i++) fhV0sCounterLambda->GetXaxis()->SetBinLabel(i+1, sV0sLambdaCounterLabel[i].Data() );
      fQAV0s->Add(fhV0sCounterLambda);

      fhV0sInvMassK0s = new TH2D("fhV0sInvMassK0s","V^{0}: K^{0}_{S}: InvMass (selected); K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sInvMassK0s);
      fhV0sInvMassLambda = new TH2D("fhV0sInvMassLambda","V^{0}: #Lambda/#bar{#Lambda}: InvMass (selected); K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sInvMassLambda);
      fhV0sCompetingInvMassK0s = new TH2D("fhV0sCompetingInvMassK0s","V^{0}: K^{0}_{S}: Competing InvMass rejection; K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sCompetingInvMassK0s);
      fhV0sCompetingInvMassLambda = new TH2D("fhV0sCompetingInvMassLambda","V^{0}: #Lambda/#bar{#Lambda}: Competing InvMass rejection; K^{0}_{S} #it{m}_{inv} (GeV/#it{c}^{2}); #Lambda/#bar{#Lambda} #it{m}_{inv} (GeV/#it{c}^{2})", 110,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax, 50,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
      fQAV0s->Add(fhV0sCompetingInvMassLambda);
    } // endif {fProcessSpec[kK0s] || fProcessSpec[kLambda]}

    const Short_t iNBinsPIDstatus = 4;
    TString sPIDstatus[iNBinsPIDstatus] = {"kDetNoSignal","kDetPidOk","kDetMismatch","kDetNoParams"};
    const Short_t iNFilterMapBinBins = 32;

    // QA histograms
    if(fFillQA)
    {
      TString sQAindex[fiNumIndexQA] = {"Before", "After"};
      for(Short_t iQA(0); iQA < fiNumIndexQA; iQA++)
      {
        // EVENTs QA histograms
        fhQAEventsPVz[iQA] = new TH1D(Form("fhQAEventsPVz_%s",sQAindex[iQA].Data()), "QA Events: PV-#it{z}", 101,-50,50);
        fQAEvents->Add(fhQAEventsPVz[iQA]);
        fhQAEventsNumContrPV[iQA] = new TH1D(Form("fhQAEventsNumContrPV_%s",sQAindex[iQA].Data()), "QA Events: Number of contributors to AOD PV", 20,0,20);
        fQAEvents->Add(fhQAEventsNumContrPV[iQA]);
        fhQAEventsNumSPDContrPV[iQA] = new TH1D(Form("fhQAEventsNumSPDContrPV_%s",sQAindex[iQA].Data()), "QA Events: SPD contributors to PV", 20,0,20);
        fQAEvents->Add(fhQAEventsNumSPDContrPV[iQA]);
        fhQAEventsDistPVSPD[iQA] = new TH1D(Form("fhQAEventsDistPVSPD_%s",sQAindex[iQA].Data()), "QA Events: PV SPD vertex", 50,0,5);
        fQAEvents->Add(fhQAEventsDistPVSPD[iQA]);
        fhQAEventsSPDresol[iQA] = new TH1D(Form("fhQAEventsSPDresol_%s",sQAindex[iQA].Data()), "QA Events: SPD resolution", 150,0,15);
        fQAEvents->Add(fhQAEventsSPDresol[iQA]);

        // Charged tracks QA
        fhQAChargedMult[iQA] = new TH1D(Form("fhQAChargedMult_%s",sQAindex[iQA].Data()),"QA Charged: Number of Charged in selected events; #it{N}^{Charged}", 1500,0,1500);
        fQACharged->Add(fhQAChargedMult[iQA]);
        fhQAChargedCharge[iQA] = new TH1D(Form("fhQAChargedCharge_%s",sQAindex[iQA].Data()),"QA Charged: Track charge; charge;", 3,-1.5,1.5);
        fQACharged->Add(fhQAChargedCharge[iQA]);
        fhQAChargedPt[iQA] = new TH1D(Form("fhQAChargedPt_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{p}_{T}; #it{p}_{T} (GeV/#it{c})", 300,0.,30.);
        fQACharged->Add(fhQAChargedPt[iQA]);
        fhQAChargedEta[iQA] = new TH1D(Form("fhQAChargedEta_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#eta}; #it{#eta}", 151,-1.5,1.5);
        fQACharged->Add(fhQAChargedEta[iQA]);
        fhQAChargedPhi[iQA] = new TH1D(Form("fhQAChargedPhi_%s",sQAindex[iQA].Data()),"QA Charged: Track #it{#varphi}; #it{#varphi}", 100,0.,TMath::TwoPi());
        fQACharged->Add(fhQAChargedPhi[iQA]);
        fhQAChargedFilterBit[iQA] = new TH1D(Form("fhQAChargedFilterBit_%s",sQAindex[iQA].Data()), "QA Charged: Filter bit",iNFilterMapBinBins,0,iNFilterMapBinBins);
        for(Int_t j = 0x0; j < iNFilterMapBinBins; j++) fhQAChargedFilterBit[iQA]->GetXaxis()->SetBinLabel(j+1, Form("%g",TMath::Power(2,j)));
        fQACharged->Add(fhQAChargedFilterBit[iQA]);
        fhQAChargedNumTPCcls[iQA] = new TH1D(Form("fhQAChargedNumTPCcls_%s",sQAindex[iQA].Data()),"QA Charged: Track number of TPC clusters; #it{N}^{TPC clusters}", 160,0,160);
        fQACharged->Add(fhQAChargedNumTPCcls[iQA]);
        fhQAChargedDCAxy[iQA] = new TH1D(Form("fhQAChargedDCAxy_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-xy; DCA_{#it{xy}} (cm)", 100,0.,10);
        fQACharged->Add(fhQAChargedDCAxy[iQA]);
        fhQAChargedDCAz[iQA] = new TH1D(Form("fhQAChargedDCAz_%s",sQAindex[iQA].Data()),"QA Charged: Track DCA-z; DCA_{#it{z}} (cm)", 200,-10.,10.);
        fQACharged->Add(fhQAChargedDCAz[iQA]);

        // PID tracks QA
        if(fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton] || fProcessSpec[kPhi])
        {
          fhQAPIDTPCstatus[iQA] = new TH1D(Form("fhQAPIDTPCstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAPID->Add(fhQAPIDTPCstatus[iQA]);
          fhQAPIDTPCdEdx[iQA] = new TH2D(Form("fhQAPIDTPCdEdx_%s",sQAindex[iQA].Data()),"QA PID: TPC PID information; #it{p} (GeV/#it{c}); TPC dEdx (au)", 100,0,10, 131,-10,1000);
          fQAPID->Add(fhQAPIDTPCdEdx[iQA]);
          fhQAPIDTOFstatus[iQA] = new TH1D(Form("fhQAPIDTOFstatus_%s",sQAindex[iQA].Data()),"QA PID: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAPID->Add(fhQAPIDTOFstatus[iQA]);
          fhQAPIDTOFbeta[iQA] = new TH2D(Form("fhQAPIDTOFbeta_%s",sQAindex[iQA].Data()),"QA PID: TOF #beta information; #it{p} (GeV/#it{c}); TOF #beta", 100,0,10, 101,-0.1,1.5);
          fQAPID->Add(fhQAPIDTOFbeta[iQA]);

          fh3QAPIDnSigmaTPCTOFPtPion[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtPion_%s",sQAindex[iQA].Data()), "QA PID: nSigma Pion vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, iPOIsPtNumBins, fFlowPOIsPtMin, fFlowPOIsPtMax);
          fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtPion[iQA]);
          fh3QAPIDnSigmaTPCTOFPtKaon[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtKaon_%s",sQAindex[iQA].Data()), "QA PID: nSigma Kaon vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, iPOIsPtNumBins, fFlowPOIsPtMin, fFlowPOIsPtMax);
          fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtKaon[iQA]);
          fh3QAPIDnSigmaTPCTOFPtProton[iQA] = new TH3D(Form("fh3QAPIDnSigmaTPCTOFPtProton_%s",sQAindex[iQA].Data()), "QA PID: nSigma Proton vs. p_{T}; n#sigma TPC; n#sigma TOF; p_{T} (GeV/c)", 21,-11,10, 21,-11,10, iPOIsPtNumBins, fFlowPOIsPtMin, fFlowPOIsPtMax);
          fQAPID->Add(fh3QAPIDnSigmaTPCTOFPtProton[iQA]);

          fh3QAPIDnSigmaBayesElectron[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesElectron_%s",sQAindex[iQA].Data()),"QA PID: e; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesElectron[iQA]);
          fh3QAPIDnSigmaBayesMuon[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesMuon_%s",sQAindex[iQA].Data()),"QA PID: #mu; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesMuon[iQA]);
          fh3QAPIDnSigmaBayesPion[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesPion_%s",sQAindex[iQA].Data()),"QA PID: #pi; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesPion[iQA]);
          fh3QAPIDnSigmaBayesKaon[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesKaon_%s",sQAindex[iQA].Data()),"QA PID: K; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesKaon[iQA]);
          fh3QAPIDnSigmaBayesProton[iQA] = new TH3D(Form("fh3QAPIDnSigmaBayesProton_%s",sQAindex[iQA].Data()),"QA PID: p; n#sigma^{TPC}; n#sigma^{TOF}; Bayes prob.", 21,-11,10, 21,-11,10, 22,-0.1,1);
          fQAPID->Add(fh3QAPIDnSigmaBayesProton[iQA]);

          for(Int_t j = 0x0; j < iNBinsPIDstatus; j++)
          {
            fhQAPIDTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
            fhQAPIDTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          }
        } // endif {fProcessSpec[kPion] || fProcessSpec[kKaon] || fProcessSpec[kProton]}

        // V0s QA
        if(fProcessSpec[kK0s] || fProcessSpec[kLambda])
        {
          fhQAV0sMultK0s[iQA] = new TH1D(Form("fhQAV0sMultK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of K^{0}_{S} candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultK0s[iQA]);
          fhQAV0sMultLambda[iQA] = new TH1D(Form("fhQAV0sMultLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #Lambda candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultLambda[iQA]);
          fhQAV0sMultALambda[iQA] = new TH1D(Form("fhQAV0sMultALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #bar{#Lambda} candidates", 1000,0,1000);
          fQAV0s->Add(fhQAV0sMultALambda[iQA]);
          fhQAV0sRecoMethod[iQA] = new TH1D(Form("fhQAV0sRecoMethod_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Reconstruction method", 2,-0.5,1.5);
          fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(1, "offline");
          fhQAV0sRecoMethod[iQA]->GetXaxis()->SetBinLabel(2, "online (on-the-fly)");
          fQAV0s->Add(fhQAV0sRecoMethod[iQA]);
          fhQAV0sDCAtoPV[iQA] = new TH1D(Form("fhQAV0sDCAtoPV_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter DCA to PV; daughter DCA^{PV} (cm)", 250,0.0,5.0);
          fQAV0s->Add(fhQAV0sDCAtoPV[iQA]);
          fhQAV0sDCADaughters[iQA] = new TH1D(Form("fhQAV0sDCADaughters_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: DCA among daughters; DCA^{daughters} (cm)", 200,0.,20.);
          fQAV0s->Add(fhQAV0sDCADaughters[iQA]);
          fhQAV0sDecayRadius[iQA] = new TH1D(Form("fhQAV0sDecayRadius_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Decay radius; #it{r_{xy}}^{decay} (cm)", 400,0.0,200.0);
          fQAV0s->Add(fhQAV0sDecayRadius[iQA]);
          fhQAV0sCPAK0s[iQA] = new TH1D(Form("fhQAV0sCPAK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: CPA; CPA^{K0s}", 100,0.9,1.);
          fQAV0s->Add(fhQAV0sCPAK0s[iQA]);
          fhQAV0sCPALambda[iQA] = new TH1D(Form("fhQAV0sCPALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: CPA; CPA^{#Lambda}", 100, 0.9,1.);
          fQAV0s->Add(fhQAV0sCPALambda[iQA]);
          fhQAV0sNumTauK0s[iQA] = new TH1D(Form("fhQAV0sNumTauK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{K0s} (cm)", 100, 0.0,60.0);
          fQAV0s->Add(fhQAV0sNumTauK0s[iQA]);
          fhQAV0sNumTauLambda[iQA] = new TH1D(Form("fhQAV0sNumTauLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Number of #it{c#tau}; #it{c#tau}^{#Lambda} (cm)", 100, 0.0,60.0);
          fQAV0s->Add(fhQAV0sNumTauLambda[iQA]);
          fhQAV0sArmenterosK0s[iQA] = new TH2D(Form("fhQAV0sArmenterosK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}:  K^{0}_{S}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosK0s[iQA]);
          fhQAV0sArmenterosLambda[iQA] = new TH2D(Form("fhQAV0sArmenterosLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosLambda[iQA]);
          fhQAV0sArmenterosALambda[iQA] = new TH2D(Form("fhQAV0sArmenterosALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Armenteros-Podolaski plot; #alpha; #it{p}_{T}^{Arm} (GeV/#it{c});", 100,-1.,1., 100,0.,0.3);
          fQAV0s->Add(fhQAV0sArmenterosALambda[iQA]);
          fhQAV0sInvMassK0s[iQA] = new TH1D(Form("fhQAV0sInvMassK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: InvMass; #it{m}_{inv} (GeV/#it{c}^{2});", 200,fCutV0sInvMassK0sMin,fCutV0sInvMassK0sMax);
          fQAV0s->Add(fhQAV0sInvMassK0s[iQA]);
          fhQAV0sInvMassLambda[iQA] = new TH1D(Form("fhQAV0sInvMassLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: InvMass; #it{m}_{inv} (GeV/#it{c}^{2});", 80,fCutV0sInvMassLambdaMin,fCutV0sInvMassLambdaMax);
          fQAV0s->Add(fhQAV0sInvMassLambda[iQA]);
          fhQAV0sMotherPt[iQA] = new TH1D(Form("fhQAV0sMotherPt_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{p}_{T}; #it{p}_{T}^{V0} (GeV/#it{c})", 200,0.,20.);
          fQAV0s->Add(fhQAV0sMotherPt[iQA]);
          fhQAV0sMotherPhi[iQA] = new TH1D(Form("fhQAV0sMotherPhi_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{#varphi}; #it{#varphi}^{V0} (GeV/#it{c})", 100,0.,TMath::TwoPi());
          fQAV0s->Add(fhQAV0sMotherPhi[iQA]);
          fhQAV0sMotherEta[iQA] = new TH1D(Form("fhQAV0sMotherEta_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{#eta}; #it{#eta}^{V0}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherEta[iQA]);
          fhQAV0sMotherCharge[iQA] = new TH1D(Form("fhQAV0sMotherCharge_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother charge; V^{0} charge", 3,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherCharge[iQA]);
          fhQAV0sMotherRapK0s[iQA] = new TH1D(Form("fhQAV0sMotherRapK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{y} (K^{0}_{S} hypo); #it{y}^{V0,K0s}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sMotherRapK0s[iQA]);
          fhQAV0sMotherRapLambda[iQA] = new TH1D(Form("fhQAV0sMotherRapLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Mother #it{y} (Lambda/#bar{#Lambda} hypo); #it{y}^{V0,#Lambda}", 301,-3.,3.);
          fQAV0s->Add(fhQAV0sMotherRapLambda[iQA]);
          fhQAV0sDaughterTPCRefit[iQA] = new TH1D(Form("fhQAV0sDaughterTPCRefit_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter TPC refit", 2,-0.5,1.5);
          fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODTrack::kTPCrefit");
          fhQAV0sDaughterTPCRefit[iQA]->GetXaxis()->SetBinLabel(2, "AliAODTrack::kTPCrefit");
          fQAV0s->Add(fhQAV0sDaughterTPCRefit[iQA]);
          fhQAV0sDaughterKinks[iQA] = new TH1D(Form("fhQAV0sDaughterKinks_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter Kinks", 2,-0.5,1.5);
          fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(1, "NOT AliAODVertex::kKink");
          fhQAV0sDaughterKinks[iQA]->GetXaxis()->SetBinLabel(2, "AliAODVertex:kKink");
          fQAV0s->Add(fhQAV0sDaughterKinks[iQA]);
          fhQAV0sDaughterNumTPCCls[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCCls_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of TPC clusters; # cls; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCCls[iQA]);
          fhQAV0sDaughterNumTPCClsPID[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCClsPID_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of TPC clusters for PID; # cls PID; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCClsPID[iQA]);
          fhQAV0sDaughterNumTPCFind[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCFind_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of findable TPC clusters; # findable; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCFind[iQA]);
          fhQAV0sDaughterNumTPCCrossRows[iQA] = new TH1D(Form("fhQAV0sDaughterNumTPCCrossRows_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter # of crossed TPC rows; # crossed; Counts;", 180,-10.,170);
          fQAV0s->Add(fhQAV0sDaughterNumTPCCrossRows[iQA]);
          fhQAV0sDaughterTPCCrossFindRatio[iQA] = new TH1D(Form("fhQAV0sDaughterTPCCrossFindRatio_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter crossed / findable TPC clusters; #crossed/find; Counts;", 50,0,5);
          fQAV0s->Add(fhQAV0sDaughterTPCCrossFindRatio[iQA]);
          fhQAV0sDaughterPt[iQA] = new TH1D(Form("fhQAV0sDaughterPt_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{p}_{T}; #it{p}_{T}^{daughter} (GeV/#it{c})", 200,0.,20.);
          fQAV0s->Add(fhQAV0sDaughterPt[iQA]);
          fhQAV0sDaughterPhi[iQA] = new TH1D(Form("fhQAV0sDaughterPhi_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{#varphi}; #it{#varphi}^{daughter} (GeV/#it{c})", 100,0.,TMath::TwoPi());
          fQAV0s->Add(fhQAV0sDaughterPhi[iQA]);
          fhQAV0sDaughterEta[iQA] = new TH1D(Form("fhQAV0sDaughterEta_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter #it{#eta}; #it{#eta}^{daugter}", 151,-1.5,1.5);
          fQAV0s->Add(fhQAV0sDaughterEta[iQA]);
          fhQAV0sDaughterCharge[iQA] = new TH1D(Form("fhQAV0sDaughterCharge_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: Daughter charge; daughter charge", 3,-1.5,1.5);
          fQAV0s->Add(fhQAV0sDaughterCharge[iQA]);
          fhQAV0sDaughterTPCstatus[iQA] = new TH1D(Form("fhQAV0sDaughterTPCstatus_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: PID status: TPC;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAV0s->Add(fhQAV0sDaughterTPCstatus[iQA]);
          fhQAV0sDaughterTOFstatus[iQA] = new TH1D(Form("fhQAV0sDaughterTOFstatus_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: PID status: TOF;", iNBinsPIDstatus,0,iNBinsPIDstatus);
          fQAV0s->Add(fhQAV0sDaughterTOFstatus[iQA]);
          fhQAV0sDaughterTPCdEdxK0s[iQA] = new TH2D(Form("fhQAV0sDaughterTPCdEdxK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: TPC dEdx daughters; #it{p}^{daughter} (GeV/#it{c}); TPC dEdx (au);", 100,0.,20, 101,-10,1000);
          fQAV0s->Add(fhQAV0sDaughterTPCdEdxK0s[iQA]);
          fhQAV0sDaughterNumSigmaPionK0s[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionK0s_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: K^{0}_{S}: Daughter PID (#pi); #it{p}_{T}^{daughter} (GeV/#it{c}); #pi PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionK0s[iQA]);
          fhQAV0sDaughterTPCdEdxLambda[iQA] = new TH2D(Form("fhQAV0sDaughterTPCdEdxLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda/#bar{#Lambda}: TPC dEdx daughters; #it{p}^{daughter} (GeV/#it{c}); TPC dEdx (au);", 100,0.,20, 101,-10,1000);
          fQAV0s->Add(fhQAV0sDaughterTPCdEdxLambda[iQA]);
          fhQAV0sDaughterNumSigmaPionLambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Daughter PID (#pi); #it{p}_{T}^{pion} (GeV/#it{c}); pion PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionLambda[iQA]);
          fhQAV0sDaughterNumSigmaProtonLambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaProtonLambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #Lambda: Daughter PID (p); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaProtonLambda[iQA]);
          fhQAV0sDaughterNumSigmaPionALambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaPionALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Daughter PID (#pi); #it{p}_{T}^{pion} (GeV/#it{c}); pion PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaPionALambda[iQA]);
          fhQAV0sDaughterNumSigmaProtonALambda[iQA] = new TH2D(Form("fhQAV0sDaughterNumSigmaProtonALambda_%s",sQAindex[iQA].Data()),"QA V^{0}_{S}: #bar{#Lambda}: Daughter PID (p); #it{p}_{T}^{proton} (GeV/#it{c}); proton PID (#sigma^{TPC});", 200,0.,20, 100,-10.,10.);
          fQAV0s->Add(fhQAV0sDaughterNumSigmaProtonALambda[iQA]);

          for(Int_t j = 0x0; j < iNBinsPIDstatus; j++)
          {
            fhQAV0sDaughterTOFstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
            fhQAV0sDaughterTPCstatus[iQA]->GetXaxis()->SetBinLabel(j+1, sPIDstatus[j].Data() );
          }
        } // endif {fProcessSpec[kK0s] || fProcessSpec[kLambda]}
      }
    }

  // posting data (mandatory)
  Int_t i = 0;
  for(Int_t iSpec(0); iSpec < kUnknown; ++iSpec) { PostData(++i, fListFlow[iSpec]); }
  PostData(++i, fQAEvents);
  PostData(++i, fQACharged);
  PostData(++i, fQAPID);
  PostData(++i, fQAV0s);
  PostData(++i, fQAPhi);
  PostData(++i, fFlowWeights);

  return;
}
