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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// This is a very simple example showing how to iterate over tracks
// and create a new collection for them.
// FIXME: this should really inherit from AnalysisTask but
//        we need GCC 7.4+ for that
struct ATask {
  float fPI = static_cast<float>(M_PI);
  Configurable<float> ptlow{"pTlow", 0.5f, "Lowest pT"};
  Configurable<float> ptup{"pTup", 2.0f, "highest pT"};
  Filter ptFilter_a = aod::track::pt > ptlow;
  Filter ptFilter_b = aod::track::pt < ptup;

  Configurable<float> etalow{"etaLow", -1.0f, "lowest eta"};
  Configurable<float> etaup{"etaUp", 1.0f, "highest eta"};
  Filter etafilter = (aod::track::eta < etaup) && (aod::track::eta > etalow);

  Configurable<float> philow{"phiLow", 1.0f, "lowest phi"};
  Configurable<float> phiup{"phiUp", 2.0f, "highest phi"};

  using myTracks = soa::Filtered<aod::Tracks>;

  Partition<myTracks> leftPhi = aod::track::phiraw < philow;
  Partition<myTracks> midPhi = aod::track::phiraw >= philow && aod::track::phiraw < phiup;
  Partition<myTracks> rightPhi = aod::track::phiraw >= phiup;

  void process(aod::Collision const& collision, myTracks const& tracks)
  {
    LOGF(INFO, "Collision: %d [N = %d] [left phis = %d] [mid phis = %d] [right phis = %d]",
         collision.globalIndex(), tracks.size(), leftPhi.size(), midPhi.size(), rightPhi.size());

    for (auto& track : leftPhi) {
      LOGF(INFO, "id = %d, from collision: %d, collision: %d; eta:  %.3f < %.3f < %.3f; phi: %.3f < %.3f; pt: %.3f < %.3f < %.3f",
           track.collisionId(), track.collision().globalIndex(), collision.globalIndex(), (float)etalow, track.eta(), (float)etaup, track.phiraw(), (float)philow, (float)ptlow, track.pt(), (float)ptup);
    }
    for (auto& track : midPhi) {
      LOGF(INFO, "id = %d, from collision: %d, collision: %d; eta: %.3f < %.3f < %.3f; phi: %.3f <= %.3f < %.3f; pt: %.3f < %.3f < %.3f",
           track.collisionId(), track.collision().globalIndex(), collision.globalIndex(), (float)etalow, track.eta(), (float)etaup, (float)philow, track.phiraw(), (float)phiup, (float)ptlow, track.pt(), (float)ptup);
    }
    for (auto& track : rightPhi) {
      LOGF(INFO, "id = %d, from collision: %d, collision: %d; eta: %.3f < %.3f < %.3f; phi: %.3f < %.3f; pt: %.3f < %.3f < %.3f",
           track.collisionId(), track.collision().globalIndex(), collision.globalIndex(), (float)etalow, track.eta(), (float)etaup, (float)phiup, track.phiraw(), (float)ptlow, track.pt(), (float)ptup);
    }
  }
};

// Partition inside process
// Caveat: partitioned table and its grouping (parent) table cannot be passed as const& to process()
struct BTask {
  void process(aod::Collisions& collisions, aod::Tracks& tracks)
  {
    for (auto& c : collisions) {
      Partition<aod::Tracks> groupedTracks = aod::track::collisionId == c.globalIndex();
      groupedTracks.bindTables(tracks, collisions);
      for (auto& t : groupedTracks) {
        LOGF(INFO, "collision global index: %d grouped track collision id: %d, collision index from track: %d", c.globalIndex(), t.collisionId(), t.collision().globalIndex());
      }
    }
  }
};

// MC example
struct CTask {
  void process(aod::McCollision const& collision, aod::McParticles& particles, aod::McMotherDaughters& motherDaughters)
  {
    LOGF(INFO, "This collision has %d particles", particles.size());

    for (auto& particle : particles) {
      LOGF(INFO, "Particle %d (pdg %d) has: ", particle.globalIndex(), particle.pdgCode());
      // uint32_t id = static_cast<uint32_t>(particle.globalIndex());
      Partition<aod::McMotherDaughters> daughters = aod::mcparticledaughter::motherId == particle.globalIndex();
      daughters.bindTables(motherDaughters, particles);

      for (auto& daughter : daughters) {
        //LOGF(INFO, "   - daughter %d (pdg %d)", daughter.globalIndex(), daughter.daughter().pdgCode());
        LOGF(INFO, "   - daughter %d (pdg %d)", daughter.globalIndex(), daughter.daughterId());
      }
    }
  }
};

struct DTask {
  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks, aod::McParticles& particles, aod::McMotherDaughters& motherDaughters)
  {
    LOGF(INFO, "This collision has %d tracks", tracks.size());

    for (auto& track : tracks) {
      auto particle = track.label(); // TODO cannot be by reference
      LOGF(INFO, "Track %d points to the MC particle %d (pdg %d)", track.globalIndex(), particle.globalIndex(), particle.pdgCode());
      Partition<aod::McMotherDaughters> mothers = aod::mcparticledaughter::daughterId == particle.globalIndex();
      mothers.bindTables(motherDaughters, particles);
      for (auto& mother : mothers) {
        LOGF(INFO, "   - mother %d (pdg %d)", mother.globalIndex(), mother.daughter().pdgCode());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    // adaptAnalysisTask<ATask>("consume-tracks"),
    // adaptAnalysisTask<BTask>("partition-in-process")
    adaptAnalysisTask<CTask>("daughter-grouping")//,
    //adaptAnalysisTask<DTask>("mother-grouping")
  };
}
