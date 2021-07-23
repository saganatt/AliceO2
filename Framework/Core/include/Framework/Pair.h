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

template <typename T, typename T2, std::size_t N, typename... REST>
struct InterleavedNTupleType {
  using type = typename InterleavedNTupleType<T, T2, N - 1, T, T2, REST...>::type;
};

template <typename T, typename T2, typename... REST>
struct InterleavedNTupleType<T, T2, 0, REST...> {
  using type = std::tuple<REST...>;
};

template <typename G, typename... As>
struct InterleavedTupleType {
  using type = typename std::tuple<(G::iterator, As)...>;
};

template<std::size_t N, typename G, typename... REST>
void execFunctionWithVariadic(void (*f)(G), const& G, REST...) {
  execFunctionWithVariadic<N-1>(f, G, G, REST...);
}

template<typename G, typename... REST>
void execFunctionWithVariadic<1>(void (*f)(G), const& G, REST...) {
  (*f)(G, REST...);
}

// TODO: with index sequence?
template<std::size_t N, typename F>
using f_return_t = typename NTupleType<F, N>::type;

template<std::size_t N, typename F>
f_return_t<N, F> functionToTuple(F (*f)()) {
  return f_return_t<N, F>((*f)(), functionToTuple<N-1, F>(f));
}

template<typename F>
f_return_t<1, F> functionToTuple(F (*f)()) {
  return f_return_t<1, F>((*f)());
}

// For correlations between particles of the same kind (track-track, v0-v0 etc.)
//template <typename G, typename A, int N = 2, typename T1 = int, typename GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<T1, G, G>>
//struct GroupedSameKindsCombinationsGenerator {
//  using grouped_iterator_t = InterleavedNTupleType<typename G::iterator, A, N>;
//};

// Most general case - if the associated tables are of different type
// E.g. track-v0 correlations

template<typename T1, typename GroupingPolicy, typename G, typename... As>
struct GroupedCombinationsGenerator {
  using grouped_iterator_t = InterleavedTupleType<G, As...>;

  struct GroupedIterator : public std::iterator<std::forward_iterator_tag, grouped_iterator_t>, public GroupingPolicy {
   public:
    using reference = grouped_iterator_t&;
    using value_type = grouped_iterator_t;
    using pointer = grouped_iterator_t*;
    using iterator_category = std::forward_iterator_tag;

    GroupedIterator(const GroupingPolicy& groupingPolicy) : GroupingPolicy(groupingPolicy) {}
    GroupedIterator(const GroupingPolicy& groupingPolicy, const G& grouping, const std::shared_ptr<GroupSlicer<G, As...>>&& slicer_ptr) : GroupingPolicy(groupingPolicy), mGrouping{std::make_shared(std::vector{grouping.asArrowTable()})}, mSlicer{std::move(slicer_ptr)}
    {
      if (!this->mIsEnd) {
        setCurrentGroupedPair();
      }
    }

    GroupedIterator(GroupedIterator const&) = default;
    GroupedIterator& operator=(GroupedIterator const&) = default;
    ~GroupedIterator() = default;

    void setTables(const G& grouping, const std::shared_ptr<GroupSlicer<G, A>>&& slicer_ptr)
    {
      mGrouping = std::make_shared<G>(std::vector{grouping.asArrowTable()});
      mSlicer.reset(slicer_ptr.get());
      execFunctionWithVariadic<sizeof...(As)>(GroupingPolicy::setTables, grouping);
      if (!this->mIsEnd) {
        setCurrentGroupedPair();
      }
    }

    void moveToEnd()
    {
      GroupingPolicy::moveToEnd();
    }

    // prefix increment
    GroupedIterator& operator++()
    {
      if (!this->mIsEnd) {
        this->addOne();
        setCurrentGroupedPair();
      }
      return *this;
    }
    // postfix increment
    GroupedIterator operator++(int /*unused*/)
    {
      GroupedIterator copy(*this);
      operator++();
      return copy;
    }
    // return reference
    reference operator*()
    {
      return *mCurrentGrouped;
    }
    bool operator==(const GroupedIterator& rh)
    {
      return (this->mIsEnd && rh.mIsEnd) || (this->mCurrent == rh.mCurrent);
    }
    bool operator!=(const GroupedIterator& rh)
    {
      return !(*this == rh);
    }

   private:
    std::tuple<As...> getAssociatedTables()
    {
      auto& currentGrouping = GroupingPolicy::mCurrent;
      constexpr auto k = sizeof...(As);
      using its_t = NTupleType<GroupSlicer::GroupSlicerIterator, k>::type;
      its_t slicerIterators = functionToTuple<k, GroupSlicer::GroupSlicerIterator>(mSlicer->mBegin);
      for (auto& slice : *mSlicer) {
        for_<k>([&, this](auto i) {
          if (slice.groupingElement().index() == std::get<i.value>(currentGrouping)->index()) {
            while(std::get<i.value>(slicerIterators) != slice) {
              std::get<i.value>(slicerIterators)++;
            }
            return;
          }
        });
      }

      // TODO: Do it properly with index_sequence
      return std::make_tuple(std::get<As>(std::get<I>(slicerIterators).associatedTables())...);
      //auto a1 = std::get<A>(it1.associatedTables());
      //auto a2 = std::get<A>(it2.associatedTables());
      //return std::make_tuple(a1, a2);
    }

    void setCurrentGroupedPair()
    {
      std::tuple<As...> associatedTables1 = getAssociatedTables();
      bool moveForward = b1.size() == 0 || b2.size() == 0;
      while (!this->mIsEnd && moveForward) {
        GroupingPolicy::addOne();
        std::tuple<As...> associatedTables2 = getAssociatedTables();
        moveForward = c1.size() == 0 || c2.size() == 0;
      }
      std::tuple<As...> associatedTables3 = getAssociatedTables();

      if (!this->mIsEnd) {
        auto& currentGrouping = GroupingPolicy::mCurrent;
        a1.bindExternalIndices(mGrouping.get());
        a2.bindExternalIndices(mGrouping.get());

        mCurrentGrouped.emplace(g1, a1, g2, a2);

        //mCurrentGrouped = GroupedIteratorRefType(g1, {a1.asArrowTable()}, g2, {a2.asArrowTable()});
        //mCurrentGrouped = grouped_iterator_t(g1, a1.begin(), g2, a2.begin());
      }
    }

    std::shared_ptr<GroupSlicer<G, A>> mSlicer = nullptr;
    std::shared_ptr<G> mGrouping = nullptr;
    std::optional<grouped_iterator_t> mCurrentGrouped;
  };

  using iterator = GroupedIterator;
  using const_iterator = GroupedIterator;

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

// Aliases for 2-particle correlations
template <typename G1, typename A1, typename G2, typename A2, typename T1 = int, typename GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<T1, G, G>>
using Pair = GroupedCombinationsGenerator<T1, GroupingPolicy, typename G1::iterator, A1, typename G2::iterator, A2>;
template <typename G, typename A, typename T1 = int, typename GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<T1, G, G>>
using Pair = GroupedCombinationsGenerator<T1, GroupingPolicy, typename G::iterator, A, typename G::iterator, A>;

// Aliases for 3-particle correlations
template <typename G1, typename A1, typename G2, typename A2, typename G3, typename A3, typename T1 = int, typename GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<T1, G, G>>
using Triple = GroupedCombinationsGenerator<T1, GroupingPolicy, typename G1::iterator, A1, typename G2::iterator, A2, typename G3::iterator, A3>;
template <typename G, typename A, typename T1 = int, typename GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<T1, G, G>>
using Triple = GroupedCombinationsGenerator<T1, GroupingPolicy, typename G::iterator, A, typename G::iterator, A, typename G::iterator, A>;

//template<typename G, typename A, typename T1>
//Pair<G, A, T1> definePair(G& grouping, const A& associated, const char* category, int catNeighbours, const T1& outsider) {
//  grouping.bindExternalIndices(&associated);
//  auto associatedTuple = std::make_tuple(associated);
//  std::shared_ptr slicer_ptr = std::make_shared<GroupSlicer<G, A>>(grouping, associatedTuple);
//  return Pair(grouping, std::move(slicer_ptr), category, catNeighbours, outsider);
//}

} // namespace o2::framework
#endif // FRAMEWORK_PAIR_H_
