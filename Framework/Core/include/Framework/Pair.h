// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef FRAMEWORK_PAIR_H
#define FRAMEWORK_PAIR_H

#include "Framework/ASoAHelpers.h"
#include "Framework/GroupSlicer.h"
#include "Framework/StringHelpers.h"

namespace o2::framework
{

template <typename AssociatedPolicy, typename G, typename A, typename T1>
struct Pair : o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<G, G>, AssociatedPolicy<A, A> {
  using PairIteratorType = std::tuple<typename G::iterator, typename A::iterator, typename G::iterator, typename A::iterator>;
  using GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<T1, G, G>;

  struct PairIterator : public std::iterator<std::forward_iterator_tag, PairIteratorType>, public GroupingPolicy, public AssociatedPolicy {
   public:
    using reference = PairIteratorType&;
    using value_type = PairIteratorType;
    using pointer = PairIteratorType*;
    using iterator_category = std::forward_iterator_tag;

    PairIterator() = default;
    PairIterator(const GroupingPolicy& groupingPolicy, const AssociatedPolicy& associatedPolicy, const G& grouping, const std::shared_ptr<GroupSlicer<G, A>>&& slicer_ptr) : GroupingPolicy(groupingPolicy), AssociatedPolicy(associatedPolicy), mGrouping{std::make_shared<G>(&grouping)}, mSlicer{std::move(slicer_ptr)} {}

    PairIterator(PairIterator const&) = default;
    PairIterator& operator=(PairIterator const&) = default;
    ~PairIterator() = default;

    std::tuple<A, A> groupAssociated()
    {
      auto& [g1, g2] = GroupingPolicy::mCurrent;
      auto it1 = mSlicer->begin();
      auto it2 = mSlicer->begin();
      for (auto& slice : *mSlicer) {
        if (slice.groupingElement().index() == g1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : *mSlicer) {
        if (slice.groupingElement().index() == g2.index()) {
          it2 = slice;
          break;
        }
      }
      auto associated1 = std::get<A>(it1.associatedTables());
      associated1.bindExternalIndices(&mGrouping);
      auto associated2 = std::get<A>(it2.associatedTables());
      associated2.bindExternalIndices(&mGrouping);
      AssociatedPolicy::setTables(associated1, associated2);
    }

    std::tuple<typename G::iterator, typename::A, typename G::iterator, typename::A> getCurrentGroupedPair()
    {
      auto [a1, a2] = AssociatedPolicy::mCurrent;
      auto [g1, g2] = GroupingPolicy::mCurrent;
      return std::make_tuple(std::move(g1), std::move(a1), std::move(g2), std::move(a2));
    }

    void moveToEnd()
    {
      GroupingPolicy::moveToEnd();
      groupAssociated();
      AssociatedPolicy::moveToEnd();
    }

    bool isEnd() {
      return GroupingPolicy::mIsEnd && AssociatedPolicy::mIsEnd;
    }

    // prefix increment
    PairIterator& operator++()
    {
      if (!isEnd()) {
        if (AssociatedPolicy::mIsEnd) {
          GroupingPolicy::addOne();
          groupAssociated();
        } else {
          AssociatedPolicy::addOne();
        }
      }
      return *this;
    }
    // postfix increment
    PairIterator operator++(int /*unused*/)
    {
      PairIterator copy(*this);
      operator++();
      return copy;
    }
    // return reference
    reference operator*()
    {
      return getCurrentGroupedPair();
    }
    bool operator==(const PairIterator& rh)
    {
      return (this->isEnd() && rh.isEnd()) || (this->mCurrentGroupedPair == rh.mCurrentGroupedPair);
    }
    bool operator!=(const PairIterator& rh)
    {
      return !(*this == rh);
    }

   private:
    std::shared_ptr<GroupSlicer<G, A>> mSlicer = nullptr;
    std::shared_ptr<G> mGrouping = nullptr;
  };

  using iterator = PairIterator;
  using const_iterator = PairIterator;

  inline iterator begin()
  {
    return iterator(mBegin);
  }
  inline iterator end()
  {
    return iterator(mEnd);
  }
  inline const_iterator begin() const
  {
    return iterator(mBegin);
  }
  inline const_iterator end() const
  {
    return iterator(mEnd);
  }

  Pair(const char* category, int catNeighbours, const T1& outsider) : mBegin(), mEnd(), mCategory(category), mCatNeighbours(catNeighbours), mOutsider(outsider) {}
  Pair(G& grouping, const std::shared_ptr<GroupSlicer<G, A>>&& slicer_ptr, const char* category, int catNeighbours, const T1& outsider) :
        mCategory(category), mCatNeighbours(catNeighbours), mOutsider(outsider),
        mBegin(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr)),
        mEnd(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr)) {}
  ~Pair() = default;

  void setTables(G& grouping, const A& associated)
  {
    //grouping.bindExternalIndices(&associated);
    //auto associatedTuple = std::make_tuple(associated);
    //std::shared_ptr slicer_ptr = std::make_shared<GroupSlicer<G, A>>(grouping, associatedTuple);
    //mBegin = iterator(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr));
    //mEnd = iterator(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr));
    //mEnd.moveToEnd();
  }

 private:
  iterator mBegin;
  iterator mEnd;
  const char* mCategory;
  const int mCatNeighbours;
  const T1 mOutsider;
};

template<typename G, typename A, typename T1>
Pair<G, A, T1> definePair(G& grouping, const A& associated, const char* category, int catNeighbours, const T1& outsider) {
  grouping.bindExternalIndices(&associated);
  auto associatedTuple = std::make_tuple(associated);
  std::shared_ptr slicer_ptr = std::make_shared<GroupSlicer<G, A>>(grouping, associatedTuple);
  return Pair(grouping, std::move(slicer_ptr), category, catNeighbours, outsider);
}

} // namespace o2::framework
#endif // FRAMEWORK_PAIR_H_
