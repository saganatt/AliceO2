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

template <typename G, typename A, typename T1>
struct Pair {
  using PairIteratorType = std::tuple<typename G::iterator, A, typename G::iterator, A>;
  using PairIteratorRefType = std::tuple<typename G::iterator&, A&, typename G::iterator&, A&>;
  using GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<T1, G, G>;

  struct PairIterator : public std::iterator<std::forward_iterator_tag, PairIteratorType>, public GroupingPolicy {
   public:
    using reference = PairIteratorType&;
    using value_type = PairIteratorType;
    using pointer = PairIteratorType*;
    using iterator_category = std::forward_iterator_tag;

    PairIterator(const GroupingPolicy& groupingPolicy) : GroupingPolicy(groupingPolicy) {}
    PairIterator(const GroupingPolicy& groupingPolicy, const G& grouping, const std::shared_ptr<GroupSlicer<G, A>>&& slicer_ptr) : GroupingPolicy(groupingPolicy), mGrouping{std::make_shared<G>(std::vector{grouping.asArrowTable()})}, mSlicer{std::move(slicer_ptr)}
    {
      if (!this->mIsEnd) {
        setCurrentGroupedPair();
      }
    }

    PairIterator(PairIterator const&) = default;
    PairIterator& operator=(PairIterator const&) = default;
    ~PairIterator() = default;

    void setTables(const G& grouping, const std::shared_ptr<GroupSlicer<G, A>>&& slicer_ptr)
    {
      mGrouping = std::make_shared<G>(std::vector{grouping.asArrowTable()});
      mSlicer.reset(slicer_ptr.get());
      GroupingPolicy::setTables(grouping, grouping);
      if (!this->mIsEnd) {
        setCurrentGroupedPair();
      }
    }

    void moveToEnd()
    {
      GroupingPolicy::moveToEnd();
    }

    // prefix increment
    PairIterator& operator++()
    {
      if (!this->mIsEnd) {
        this->addOne();
        setCurrentGroupedPair();
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
      return *mCurrentGroupedPair;
    }
    bool operator==(const PairIterator& rh)
    {
      return (this->mIsEnd && rh.mIsEnd) || (this->mCurrent == rh.mCurrent);
    }
    bool operator!=(const PairIterator& rh)
    {
      return !(*this == rh);
    }

   private:
    std::tuple<A, A> getAssociatedTables()
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
      auto a1 = std::get<A>(it1.associatedTables());
      auto a2 = std::get<A>(it2.associatedTables());
      return std::make_tuple(a1, a2);
    }

    void setCurrentGroupedPair()
    {
      auto [b1, b2] = getAssociatedTables();
      bool moveForward = b1.size() == 0 || b2.size() == 0;
      while (!this->mIsEnd && moveForward) {
        GroupingPolicy::addOne();
        auto [c1, c2] = getAssociatedTables();
        moveForward = c1.size() == 0 || c2.size() == 0;
      }
      auto [a1, a2] = getAssociatedTables();

      if (!this->mIsEnd) {
        auto& [g1, g2] = GroupingPolicy::mCurrent;
        a1.bindExternalIndices(mGrouping.get());
        a2.bindExternalIndices(mGrouping.get());

        mCurrentGroupedPair.emplace(g1, a1, g2, a2);

        //mCurrentGroupedPair = PairIteratorRefType(g1, {a1.asArrowTable()}, g2, {a2.asArrowTable()});
        //mCurrentGroupedPair = PairIteratorType(g1, a1.begin(), g2, a2.begin());
      }
    }

    std::shared_ptr<GroupSlicer<G, A>> mSlicer = nullptr;
    std::shared_ptr<G> mGrouping = nullptr;
    std::optional<PairIteratorType> mCurrentGroupedPair;
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

  Pair(const char* category, int catNeighbours, const T1& outsider) : mBegin(GroupingPolicy(category, catNeighbours, outsider)), mEnd(GroupingPolicy(category, catNeighbours, outsider)), mCategory(category), mCatNeighbours(catNeighbours), mOutsider(outsider) {}
  Pair(G& grouping, const std::shared_ptr<GroupSlicer<G, A>>&& slicer_ptr, const char* category, int catNeighbours, const T1& outsider) : mCategory(category), mCatNeighbours(catNeighbours), mOutsider(outsider), mBegin(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr)), mEnd(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr))
  {
    mEnd.moveToEnd();
  }
  ~Pair() = default;

  void setTables(G& grouping, const A& associated)
  {
    grouping.bindExternalIndices(&associated);
    auto associatedTuple = std::make_tuple(associated);
    std::shared_ptr slicer_ptr = std::make_shared<GroupSlicer<G, A>>(grouping, associatedTuple);
    mBegin.setTables(grouping, std::move(slicer_ptr));
    mEnd.setTables(grouping, std::move(slicer_ptr));
    //mBegin = iterator(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr));
    //mEnd = iterator(GroupingPolicy(mCategory, mCatNeighbours, mOutsider, grouping, grouping), grouping, std::move(slicer_ptr));
    mEnd.moveToEnd();
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
