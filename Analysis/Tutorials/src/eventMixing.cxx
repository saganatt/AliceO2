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
#include "Framework/ASoAHelpers.h"

#include <unordered_set>

namespace o2::aod
{
namespace hash
{
DECLARE_SOA_COLUMN(Bin, bin, int);
} // namespace hash
DECLARE_SOA_TABLE(Hashes, "AOD", "HASH", hash::Bin);
namespace hashedCollision
{
DECLARE_SOA_INDEX_COLUMN_FULL(Hash, hash, int32_t, Hashes, "fHashesID");
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
} // namespace hashedCollision
DECLARE_SOA_TABLE(HashedCollisions, "AOD", "HASHEDCOLLISION", o2::soa::Index<>, hashedCollision::HashId, hashedCollision::CollisionId);

using Hash = Hashes::iterator;
using HashedCollision = HashedCollisions::iterator;
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;

// This is a tentative workflow to get mixed-event tracks
// FIXME: this should really inherit from AnalysisTask but
//        we need GCC 7.4+ for that

struct HashTask {
  std::vector<float> xBins{-1.5f, -1.0f, -0.5f, 0.0f, 0.5f, 1.0f, 1.5f};
  std::vector<float> yBins{-1.5f, -1.0f, -0.5f, 0.0f, 0.5f, 1.0f, 1.5f};
  Produces<aod::Hashes> hashes;
  // Produces<aod::HashedCollisions> hashedCollisions; // For version 2

  // Calculate hash for an element based on 2 properties and their bins.
  int getHash(std::vector<float> const& xBins, std::vector<float> const& yBins, float colX, float colY)
  {
    if (colX < xBins[0] || colY < yBins[0]) {
      return -1;
    }
    for (int i = 1; i < xBins.size(); i++) {
      if (colX < xBins[i]) {
        for (int j = 1; j < yBins.size(); j++) {
          if (colY < yBins[j]) {
            return i + j * (xBins.size() + 1);
          }
        }
        return -1;
      }
    }

    return -1;
  }

  // For version 1
  void process(aod::Collisions const& collisions)
  {
    for (auto& collision : collisions) {
      int hash = getHash(xBins, yBins, collision.posX(), collision.posY());
      LOGF(info, "Collision: %d (%f, %f, %f) hash: %d", collision.index(), collision.posX(), collision.posY(), collision.posZ(), hash);
      hashes(hash);
    }
  }

  // For version 2
  //void process2(aod::Collisions const& collisions)
  //{
  //  int ind = 0;
  //  std::unordered_set<int> uniqueHashes;
  //  for (auto& collision : collisions) {
  //    int hash = getHash(xBins, yBins, collision.posX(), collision.posY());
  //    LOGF(info, "Collision: %d (%f, %f, %f) hash: %d", collision.index(), collision.posX(), collision.posY(), collision.posZ(), hash);
  //    uniqueHashes.emplace(hash);
  //    hashedCollisions(hash, collision);
  //  }
  //  for (auto& hash : uniqueHashes) {
  //    hashes(hash);
  //  }
  //}
};

// Version 1: Using categorised combinations
struct CollisionsCombinationsTask {
  o2::framework::expressions::Filter trackFilter = (aod::track::x > -0.8f) && (aod::track::x < 0.8f) && (aod::track::y > 1.0f);

  void process(aod::Hashes const& hashes, aod::Collisions& collisions, soa::Filtered<aod::Tracks>& tracks)
  {
    collisions.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    AnalysisDataProcessorBuilder::GroupSlicer slicer(collisions, tracksTuple);

    // Strictly upper categorised collisions
    for (auto& [c1, c2] : selfCombinations("fBin", 5, -1, join(hashes, collisions), join(hashes, collisions))) {
      LOGF(info, "Collisions bin: %d pair: %d (%f, %f, %f), %d (%f, %f, %f)", c1.bin(), c1.index(), c1.posX(), c1.posY(), c1.posZ(), c2.index(), c2.posX(), c2.posY(), c2.posZ());

      auto it1 = slicer.begin();
      auto it2 = slicer.begin();
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == c1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : slicer) {
        if (slice.groupingElement().index() == c2.index()) {
          it2 = slice;
          break;
        }
      }
      auto tracks1 = std::get<soa::Filtered<aod::Tracks>>(it1.associatedTables());
      tracks1.bindExternalIndices(&collisions);
      auto tracks2 = std::get<soa::Filtered<aod::Tracks>>(it2.associatedTables());
      tracks2.bindExternalIndices(&collisions);

      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.index(), t2.index(), c1.index(), c2.index());
      }
    }
  }
};

// Version 2: Using nested grouping
struct NestedGroupingTask {
  void process(aod::Hash& hash, aod::HashedCollisions& hashedCollisions, aod::Tracks& tracks)
  {
    //    hashes.bindExternalIndices(&hashedCollisions);
    //    auto hashedTuple = std::make_tuple(hashedCollisions);
    //    AnalysisDataProcessorBuilder::GroupSlicer hashSlicer(hashes, hashedTuple);
    //
    //    for (auto& slice : hashSlicer) {
    //      auto sliceCollisions = std::get<aod::HashedCollisions>(slice.associatedTables());
    //      sliceCollisions.bindExternalIndices(&tracks);
    hashedCollisions.bindExternalIndices(&tracks);
    auto tracksTuple = std::make_tuple(tracks);
    AnalysisDataProcessorBuilder::GroupSlicer colSlicer(hashedCollisions, tracksTuple);

    for (auto it1 = hashedCollisions.begin(); it1 != hashedCollisions.end(); it1++) {
      auto& c1 = *it1;
      for (auto it2 = it1 + 1; it2 != hashedCollisions.end(); it2++) {
        auto& c2 = *it2;
        LOGF(info, "Collisions bin: %d pair: %d (%f, %f, %f), %d (%f, %f, %f)", c1.hash().bin(), c1.collision().index(), c1.collision().posX(), c1.collision().posY(), c1.collision().posZ(), c2.collision().index(), c2.collision().posX(), c2.collision().posY(), c2.collision().posZ());

        auto tit1 = colSlicer.begin();
        auto tit2 = colSlicer.begin();
        for (auto& colSlice : colSlicer) {
          if (colSlice.groupingElement().collision().index() == c1.collision().index()) {
            tit1 = colSlice;
            break;
          }
        }
        for (auto& colSlice : colSlicer) {
          if (colSlice.groupingElement().collision().index() == c2.collision().index()) {
            tit2 = colSlice;
            break;
          }
        }
        auto tracks1 = std::get<aod::Tracks>(tit1.associatedTables());
        tracks1.bindExternalIndices(&hashedCollisions);
        auto tracks2 = std::get<aod::Tracks>(tit2.associatedTables());
        tracks2.bindExternalIndices(&hashedCollisions);

        for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
          LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.index(), t2.index(), c1.collision().index(), c2.collision().index());
        }
      }
    }
    //    }
  }
};

// Version 3: Filtering instead of combinations
// Possible only after filters & grouping upgrades
// struct CollisionsFilteringTask {
// Alternatively: filter/partition directly on collisions
// expressions::Filter aod::hash::bin{0} == aod::hash::bin{1};

// Currently multiple grouping and grouping by Joins is not possible
// void process(soa::Filtered<soa::Join<aod::Hashes, aod::Collisions>>::iterator const& hashedCol1, aod::Tracks const& tracks1,
//              soa::Filtered<soa::Join<aod::Hashes, aod::Collisions>>::iterator const& hashedCol2, aod::Tracks const& tracks2)
//{
//  for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
//    LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.index(), t2.index(), hashedCol1.index(), hashedCol2.index());
//  }
//}
// };

// What we would like to have
struct MixedEventsTask {
  void process(aod::Collision const& col1, aod::Tracks const& tracks1, aod::Collision const& col2, aod::Tracks const& tracks2)
  {
    for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
      LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d)", t1.index(), t2.index(), col1.index(), col2.index());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<HashTask>("collisions-hashed"),
    adaptAnalysisTask<CollisionsCombinationsTask>("mixed-event-tracks")};
  //adaptAnalysisTask<NestedGroupingTask>("mixed-event-tracks")};
}
