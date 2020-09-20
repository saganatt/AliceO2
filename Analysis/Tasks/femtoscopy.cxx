// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

#include "Analysis/TrackSelectionTables.h"
#include "Analysis/EventSelection.h"
#include "Analysis/Centrality.h"
#include "Analysis/Multiplicity.h"
#include "Analysis/StepTHn.h"
#include "Analysis/CorrelationContainer.h"

#include <TH1F.h>
#include <cmath>
#include <TDirectory.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FemtoscopyTask {
  using myCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Cents, aod::Mults>;
  using myCollision = soa::Join<aod::Collisions, aod::EvSels, aod::Cents, aod::Mults>::iterator;
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  const double PionMass = 0.13956995;

  bool performSharedDaughterCut = false;
  bool enablePairMonitors = true;

  // TODO: Just one kind for pions?
  // Pions: 6
  const char* chrgs[numOfChTypes] = {"PIpPIp"};

  int runqinv = 1;

  int runtype = 0; // Types 0 - global, 1 - ITS only, 2 - TPC Inner
  int owncuts = 0;
  int owndca = 0;

  int gammacut = 1; // cut na ee z gamma

  double shqmax = 2.0;
  int nbinssh = 20;

  // Configuration
  O2_DEFINE_CONFIGURABLE(cfgFilterBit, int, 96, "Select filterbit");
  // TODO: Corresponding cut in O2?
  O2_DEFINE_CONFIGURABLE(cfgMinPlpContribSPD, int, 3, "SPD pile-up cut");

  // TODO: Boolean once enabled
  O2_DEFINE_CONFIGURABLE(cfgIfGlobalTracks, int, 0, "Whether to use global tracks: 1 - yes, 0 - no");
  O2_DEFINE_CONFIGURABLE(cfgIfElectronRejection, int, 1, "Whether to reject electrons: 1 - yes, 0 - no");
  O2_DEFINE_CONFIGURABLE(cfgIfIsPileUp, int, 1, "Whether to pile up: 1 - yes, 0 - no");
  // TODO: Are there 'monitors' in O2? What is this?
  O2_DEFINE_CONFIGURABLE(cfgIfMonitors, int, 1, "Whether to use monitors: 1 - yes, 0 - no");

  O2_DEFINE_CONFIGURABLE(cfgNSigmaVal, float, 2.0f, "nSigmaVal");
  O2_DEFINE_CONFIGURABLE(cfgNSigmaVal2, float, 3.0f, "nSigmaVal2");
  O2_DEFINE_CONFIGURABLE(cfgNEtaMin, float, -0.8f, "nEtaMin");
  O2_DEFINE_CONFIGURABLE(cfgNEtaMax, float, 0.8f, "nEtaMax");
  O2_DEFINE_CONFIGURABLE(cfgMaxPt, float, 2.5f, "MaxPt");

  // TODO: AliFemtoEventReaderAOD::CopyAODToFemtoEvent() - new task and femto data model?

  // Additional codes from MultSelection reader:
  // See AliPhysics/PWGCF/FEMTOSCOPY/AliFemto/AliFemtoEventReaderAODMultSelection.cxx:25
  // Data from MultSelection task
  // femto_event->SetCentralityV0(mult_selection->GetMultiplicityPercentile("V0M"));
  // TODO: CL1 not present in O2?
  // femto_event->SetCentralityCL1(mult_selection->GetMultiplicityPercentile("CL1"));
  // TODO: What fNormalizedMult is used for?
  // femto_event->SetNormalizedMult(lrint(10 * mult_selection->GetMultiplicityPercentile("V0M")));

  Filter filterBit96 = (aod::track::isGlobalTrack == (uint8_t)1) || (aod::track::isGlobalTrackSDD == (uint8_t)1);
  // From correlations.cxx
  // Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt) && (aod::track::isGlobalTrack == (uint8_t)1);

  // TODO: Cuts

  // Basic cut and monitor
  Filter basicEventCut = (aod::mult::multV0M >= 0.001) && (aod::mult::multV0M <= 100000) && 
                         (aod::collision::posZ > -7.0) && (aod::collision::posZ < 7.0);

  struct cutMonitorEventMult : public HistogramRegistry {
    cutMonitorEventMult(char const* const name, const std::string& histSuffix, int nBins=5000, double multMax=5000.5)
      : HistogramRegistry(name, true, {//
          {"EvMult" + histSuffix, "Event Multiplicity", {HistogramType::kTH1D,//
            {nBins + 1, -0.5, multMax}//
          }},//
          {"NormEvMult" + histSuffix, "Normalized Event Multiplicity", {HistogramType::kTH1D,//
            {nBins + 1, -0.5, multMax}//
          }},//
          {"PsiEPVZERO" + histSuffix, "Event Plane Angle from VZero", {HistogramType::kTH1D,//
            {157, -1.575, 1.565}//
          }}//
        }) {}
  };
  //fEvMult->Fill(aEvent->NumberOfTracks());
  //fNormEvMult->Fill(aEvent->UncorrectedNumberOfPrimaries());
  //fPsiVZERO->Fill(aEvent->ReactionPlaneAngle());
  cutMonitorEventMult cutPassEvMetaphitpc{"cutPassEvMetaphitpc", fmt::format("cutPass%stpcM0", chrgs[ichg]), 2000, 20000.5}; 
  cutMonitorEventMult cutFailEvMetaphitpc{"cutFailEvMetaphitpc", fmt::format("cutFail%stpcM0", chrgs[ichg]), 2000, 20000.5}; 

  // Single particle track cuts
  //dtc1etaphitpc[aniter] = new AliFemtoMJTrackCut();
  //dtc1etaphitpc[aniter]->SetCharge(1.0);
  //dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
  //dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal); //w zaleznosci od pedow jest brane albo nsigma1 albo nsigma2
  //dtc1etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
  //dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);//bierzemy jednoczesnie sigma z TPC i TOF jednoczesnie
  //dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);

  // TODO: Monitors (histograms)
  // TODO: Correlation functions

  // TODO: Event mixing according to multiplicity and z-vertex
  //                          bins vtx, min vtx, max vtx, bins mult, min mult, max mult
  // 5 events to mix, collection of size min 1
  // AliFemtoVertexMultAnalysis(7, -7.0, 7.0, 20, 0, 1000);

  void process(myCollisions const& collision, myTracks const& tracks)
  {
    LOGF(info, "Tracks for collision: %d | Trigger mask: %lld | INT7: %d | V0M: %.1f | multV0M: %5.0f", tracks.size(), collision.bc().triggerMask(), collision.sel7(), collision.centV0M(), collision.multV0M());

    if (!collision.sel7())
      return;

    int bSign = 1; // TODO magnetic field from CCDB
    const float pTCut = 1.0;

    for (auto track1 = tracks.begin(); track1 != tracks.end(); ++track1) {

#ifndef MYFILTER
      if (track1.pt2() < pTCut)
        continue;
      if (track1.eta2() < -0.8 || track1.eta2() > 0.8)
        continue;
#endif

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.charge() < 0)
        continue;

      //LOGF(info, "TRACK %f %f | %f %f | %f %f", track1.eta(), track1.eta2(), track1.phi(), track1.phi2(), track1.pt(), track1.pt2());

      double eventValues[3];
      eventValues[0] = track1.pt();
      eventValues[1] = collision.centV0M();
      eventValues[2] = collision.posZ();

      same->getEventHist()->Fill(eventValues, CorrelationContainer::kCFStepReconstructed);
      //mixed->getEventHist()->Fill(eventValues, CorrelationContainer::kCFStepReconstructed);

      for (auto track2 = track1 + 1; track2 != tracks.end(); ++track2) {
#ifndef MYFILTER
        if (track2.pt2() < pTCut)
          continue;
        if (track2.eta2() < -0.8 || track2.eta2() > 0.8)
          continue;
#endif

        if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * track2.charge() < 0)
          continue;
        if (cfgPairCharge != 0 && cfgPairCharge * track1.charge() * track2.charge() < 0)
          continue;

        if (cfg.mPairCuts && conversionCuts(track1, track2))
          continue;

        if (cfgTwoTrackCut > 0 && twoTrackCut(track1, track2, bSign))
          continue;

        double values[6] = {0};

        values[0] = track1.etam() - track2.etam();
        values[1] = track1.pt();
        values[2] = track2.pt();
        values[3] = collision.centV0M();

        values[4] = track1.phim() - track2.phim();
        if (values[4] > 1.5 * TMath::Pi())
          values[4] -= TMath::TwoPi();
        if (values[4] < -0.5 * TMath::Pi())
          values[4] += TMath::TwoPi();

        values[5] = collision.posZ();

        same->getTrackHist()->Fill(values, CorrelationContainer::kCFStepReconstructed);
      }
    }
  }
};

// TODO: PID response task??

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<FemtoscopyTask>("femtoscopy-task")};
}
