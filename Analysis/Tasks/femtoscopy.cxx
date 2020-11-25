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

namespace o2::aod
{
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int);
} // namespace hash
DECLARE_SOA_TABLE(Hashes, "AOD", "HASH", hash::Bin);

using Hash = Hashes::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Event mixing according to multiplicity and z-vertex
//                          bins vtx, min vtx, max vtx, bins mult, min mult, max mult
// AliFemtoVertexMultAnalysis(7, -7.0, 7.0, 20, 0, 1000);
struct HashTask {
  std::vector<float> vtxBins{-7.0f, -5.0f, -3.0f, -1.0f, 1.0f, 3.0f, 5.0f, 7.0f};
  std::vector<float> multBins{0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f};
  Produces<aod::Hashes> hashes;

  // Calculate hash for an element based on 2 properties and their bins.
  int getHash(std::vector<float> const& vtxBins, std::vector<float> const& multBins, float vtx, float mult)
  {
    // underflow
    if (vtx < vtxBins[0]) {
      return -1;
    }
    if (mult < multBins[0]) {
      return -1;
    }

    for (int i = 1; i < vtxBins.size(); i++) {
      if (vtx < vtxBins[i]) {
        for (int j = 1; j < multBins.size(); j++) {
          if (mult < multBins[j]) {
            return i + j * (vtxBins.size() + 1);
          }
        }
      }
    }
    // overflow
    return -1;
  }

  void process(soa::Join<aod::Collisions, aod::Cents> const& collisions)
  {
    for (auto& collision : collisions) {
      int hash = getHash(vtxBins, multBins, collision.posZ(), collision.centV0M());
      LOGF(info, "Collision: %d (%f, %f) hash: %d", collision.index(), collision.posZ(), collision.centV0M(), hash);
      hashes(hash);
    }
  }
};

// Pion+pion- from pp events, deta dphi correlations
struct FemtoscopyTask {
  using myCollisions = soa::Join<aod::Collisions, aod::Hashes, aod::EvSels, aod::Cents, aod::Mults>;
  using myCollision = soa::Join<aod::Collisions, aod::Hashes, aod::EvSels, aod::Cents, aod::Mults>::iterator;
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

  // Basic cut
  Filter basicEventCut = (aod::mult::multV0M >= 0.001) && (aod::mult::multV0M <= 100000) && 
                         (aod::collision::posZ > -7.0) && (aod::collision::posZ < 7.0);

  //fEvMult->Fill(aEvent->NumberOfTracks());
  //fNormEvMult->Fill(aEvent->UncorrectedNumberOfPrimaries());
  //fPsiVZERO->Fill(aEvent->ReactionPlaneAngle());
  //cutMonitorEventMult cutPassEvMetaphitpc{"cutPassEvMetaphitpc", fmt::format("cutPass%stpcM0", chrgs[ichg]), 2000, 20000.5};
  //cutMonitorEventMult cutFailEvMetaphitpc{"cutFailEvMetaphitpc", fmt::format("cutFail%stpcM0", chrgs[ichg]), 2000, 20000.5};

  // Single particle track cuts - to be implemented as filters
  //dtc1etaphitpc[aniter] = new AliFemtoMJTrackCut();
  //dtc1etaphitpc[aniter]->SetCharge(1.0);
  //dtc1etaphitpc[aniter]->SetEta(nEtaMin,nEtaMax);
  //dtc1etaphitpc[aniter]->SetNsigma(nSigmaVal); //w zaleznosci od pedow jest brane albo nsigma1 albo nsigma2
  //dtc1etaphitpc[aniter]->SetNsigma2(nSigmaVal2);
  //dtc1etaphitpc[aniter]->SetNsigmaTPCTOF(kTRUE);//bierzemy jednoczesnie sigma z TPC i TOF jednoczesnie
  //dtc1etaphitpc[aniter]->SetElectronRejection(ifElectronRejection);

  // TODO: Monitors (histograms)
  // TODO: Correlation functions

  void process(myCollisions const& collisions, myTracks const& tracks)
  {
    // Event mixing according to multiplicity and z-vertex
    // FIXME: Hide internals when possible
    collisions.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    AnalysisDataProcessorBuilder::GroupSlicer slicer(collisions, tracksTuple);
    // Strictly upper categorised collisions, for 5 combinations per bin, skipping those in entry -1
    for (auto& [collision1, collision2] : selfCombinations("fBin", 5, -1, collisions, collisions)) {

      LOGF(INFO, "Collisions bin: %d pair: %d (%f), %d (%f)", collision1.bin(), collision1.index(), collision1.posZ(), collision2.index(), collision2.posZ());

      auto it1 = slicer.begin();
      auto it2 = slicer.begin();
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == collision2.index()) {
          it2 = slice;
          break;
        }
      }

      auto tracks1 = std::get<myTracks>(it1.associatedTables());
      tracks1.bindExternalIndices(&collisions);
      auto tracks2 = std::get<myTracks>(it2.associatedTables());
      tracks2.bindExternalIndices(&collisions);

      for (auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(INFO, "Track 1: %d 2: %d", track1.index(), track2.index());

        // TODO: Single particle cuts
        // TODO: Pair cuts
        // TODO: Histogram filling
      }
    }
  }
};

// TODO: PID response task??

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<HashTask>("produce-hashes"),
    adaptAnalysisTask<FemtoscopyTask>("femtoscopy-task")};
}
