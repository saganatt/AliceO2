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

#ifndef FRAMEWORK_GROUP_SLICER_H_
#define FRAMEWORK_GROUP_SLICER_H_

#include "Framework/Pack.h"
#include "Framework/Kernels.h"

#include <arrow/util/key_value_metadata.h>
#include <type_traits>
#include <string>

namespace o2::framework
{

template <typename G, typename... A>
struct GroupSlicer {
  using grouping_t = std::decay_t<G>;
  GroupSlicer(G& gt, std::tuple<A...>& at)
    : max{gt.size()},
      mBegin{GroupSlicerIterator(gt, at)}
  {
  }

  struct GroupSlicerSentinel {
    int64_t position;
  };

  struct GroupSlicerIterator {
    using associated_pack_t = framework::pack<A...>;

    GroupSlicerIterator() = default;
    GroupSlicerIterator(GroupSlicerIterator const&) = default;
    GroupSlicerIterator(GroupSlicerIterator&&) = default;
    GroupSlicerIterator& operator=(GroupSlicerIterator const&) = default;
    GroupSlicerIterator& operator=(GroupSlicerIterator&&) = default;

    auto getLabelFromType()
    {
      if constexpr (soa::is_soa_index_table_t<std::decay_t<G>>::value) {
        using T = typename std::decay_t<G>::first_t;
        if constexpr (soa::is_type_with_originals_v<std::decay_t<T>>) {
          using O = typename framework::pack_element_t<0, typename std::decay_t<G>::originals>;
          using groupingMetadata = typename aod::MetadataTrait<O>::metadata;
          return std::string("fIndex") + groupingMetadata::tableLabel();
        } else {
          using groupingMetadata = typename aod::MetadataTrait<T>::metadata;
          return std::string("fIndex") + groupingMetadata::tableLabel();
        }
      } else if constexpr (soa::is_type_with_originals_v<std::decay_t<G>>) {
        using T = typename framework::pack_element_t<0, typename std::decay_t<G>::originals>;
        using groupingMetadata = typename aod::MetadataTrait<T>::metadata;
        return std::string("fIndex") + groupingMetadata::tableLabel();
      } else {
        using groupingMetadata = typename aod::MetadataTrait<std::decay_t<G>>::metadata;
        return std::string("fIndex") + groupingMetadata::tableLabel();
      }
    }

    GroupSlicerIterator(G& gt, std::tuple<A...>& at)
      : mAt{&at},
        mGroupingElement{gt.begin()},
        position{0}
    {
      if constexpr (soa::is_soa_filtered_t<std::decay_t<G>>::value) {
        groupSelection = &gt.getSelectedRows();
      }
      auto indexColumnName = getLabelFromType();
      /// prepare slices and offsets for all associated tables that have index
      /// to grouping table
      ///
      auto splitter = [&](auto&& x) {
        using xt = std::decay_t<decltype(x)>;
        constexpr auto index = framework::has_type_at_v<std::decay_t<decltype(x)>>(associated_pack_t{});
        if (x.size() != 0 && hasIndexTo<std::decay_t<G>>(typename xt::persistent_columns_t{})) {
          auto result = o2::framework::sliceByColumn(indexColumnName.c_str(),
                                                     x.asArrowTable(),
                                                     static_cast<int32_t>(gt.tableSize()),
                                                     &groups[index],
                                                     &offsets[index],
                                                     &sizes[index]);
          if (result.ok() == false) {
            throw runtime_error("Cannot split collection");
          }
          if (groups[index].size() > gt.tableSize()) {
            throw runtime_error_f("Splitting collection resulted in a larger group number (%d) than there is rows in the grouping table (%d).", groups[index].size(), gt.tableSize());
          };
        }
      };

      std::apply(
        [&](auto&&... x) -> void {
          (splitter(x), ...);
        },
        at);
      /// extract selections from filtered associated tables
      auto extractor = [&](auto&& x) {
        using xt = std::decay_t<decltype(x)>;
        if constexpr (soa::is_soa_filtered_t<xt>::value) {
          constexpr auto index = framework::has_type_at_v<std::decay_t<decltype(x)>>(associated_pack_t{});
          selections[index] = &x.getSelectedRows();
          starts[index] = selections[index]->begin();
          offsets[index].push_back(std::get<xt>(at).tableSize());
        }
      };
      std::apply(
        [&](auto&&... x) -> void {
          (extractor(x), ...);
        },
        at);
    }

    template <typename B, typename... C>
    constexpr bool hasIndexTo(framework::pack<C...>&&)
    {
      return (isIndexTo<B, C>() || ...);
    }

    template <typename B, typename C>
    constexpr bool isIndexTo()
    {
      if constexpr (soa::is_type_with_binding_v<C>) {
        if constexpr (soa::is_soa_index_table_t<std::decay_t<B>>::value) {
          using T = typename std::decay_t<B>::first_t;
          if constexpr (soa::is_type_with_originals_v<std::decay_t<T>>) {
            using TT = typename framework::pack_element_t<0, typename std::decay_t<T>::originals>;
            return std::is_same_v<typename C::binding_t, TT>;
          } else {
            using TT = std::decay_t<T>;
            return std::is_same_v<typename C::binding_t, TT>;
          }
        } else {
          if constexpr (soa::is_type_with_originals_v<std::decay_t<B>>) {
            using TT = typename framework::pack_element_t<0, typename std::decay_t<B>::originals>;
            return std::is_same_v<typename C::binding_t, TT>;
          } else {
            using TT = std::decay_t<B>;
            return std::is_same_v<typename C::binding_t, TT>;
          }
        }
      }
      return false;
    }

    GroupSlicerIterator& operator++()
    {
      ++position;
      ++mGroupingElement;
      return *this;
    }

    bool operator==(GroupSlicerSentinel const& other)
    {
      return O2_BUILTIN_UNLIKELY(position == other.position);
    }

    bool operator!=(GroupSlicerSentinel const& other)
    {
      return O2_BUILTIN_LIKELY(position != other.position);
    }

    auto& groupingElement()
    {
      return mGroupingElement;
    }

    GroupSlicerIterator& operator*()
    {
      return *this;
    }

    auto associatedTables()
    {
      return std::make_tuple(prepareArgument<A>()...);
    }

    template <typename A1>
    auto prepareArgument()
    {
      constexpr auto index = framework::has_type_at_v<A1>(associated_pack_t{});
      if (std::get<A1>(*mAt).size() != 0 && hasIndexTo<G>(typename std::decay_t<A1>::persistent_columns_t{})) {
        uint64_t pos;
        if constexpr (soa::is_soa_filtered_t<std::decay_t<G>>::value) {
          pos = (*groupSelection)[position];
        } else {
          pos = position;
        }
        if constexpr (soa::is_soa_filtered_t<std::decay_t<A1>>::value) {
          auto groupedElementsTable = arrow::util::get<std::shared_ptr<arrow::Table>>(((groups[index])[pos]).value);

          // for each grouping element we need to slice the selection vector
          auto start_iterator = std::lower_bound(starts[index], selections[index]->end(), (offsets[index])[pos]);
          auto stop_iterator = std::lower_bound(start_iterator, selections[index]->end(), (offsets[index])[pos] + (sizes[index])[pos]);
          starts[index] = stop_iterator;
          soa::SelectionVector slicedSelection{start_iterator, stop_iterator};
          std::transform(slicedSelection.begin(), slicedSelection.end(), slicedSelection.begin(),
                         [&](int64_t idx) {
                           return idx - static_cast<int64_t>((offsets[index])[pos]);
                         });

          std::decay_t<A1> typedTable{{groupedElementsTable}, std::move(slicedSelection), (offsets[index])[pos]};
          return typedTable;
        } else {
          auto groupedElementsTable = arrow::util::get<std::shared_ptr<arrow::Table>>(((groups[index])[pos]).value);
          std::decay_t<A1> typedTable{{groupedElementsTable}, (offsets[index])[pos]};
          return typedTable;
        }
      } else {
        return std::get<A1>(*mAt);
      }
      O2_BUILTIN_UNREACHABLE();
    }

    std::tuple<A...>* mAt;
    typename grouping_t::iterator mGroupingElement;
    uint64_t position = 0;
    soa::SelectionVector const* groupSelection = nullptr;
    std::array<std::vector<arrow::Datum>, sizeof...(A)> groups;
    std::array<std::vector<uint64_t>, sizeof...(A)> offsets;
    std::array<std::vector<int>, sizeof...(A)> sizes;
    std::array<soa::SelectionVector const*, sizeof...(A)> selections;
    std::array<soa::SelectionVector::const_iterator, sizeof...(A)> starts;
  };

  GroupSlicerIterator& begin()
  {
    return mBegin;
  }

  GroupSlicerSentinel end()
  {
    return GroupSlicerSentinel{max};
  }
  int64_t max;
  GroupSlicerIterator mBegin;
};

} // namespace o2::framework
#endif // FRAMEWORK_GROUP_SLICER_H_
