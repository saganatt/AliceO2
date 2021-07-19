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
#ifndef o2_framework_AnalysisHelpers_H_DEFINED
#define o2_framework_AnalysisHelpers_H_DEFINED

#include "Framework/Traits.h"
#include "Framework/TableBuilder.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/OutputSpec.h"
#include "Framework/OutputRef.h"
#include "Framework/InputSpec.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/StringHelpers.h"
#include "Framework/Output.h"
#include <string>
#include "Framework/Logger.h"

namespace o2::framework
{
class TableConsumer;

template <typename T>
struct WritingCursor {
  static_assert(always_static_assert_v<T>, "Type must be a o2::soa::Table");
};
/// Helper class actually implementing the cursor which can write to
/// a table. The provided template arguments are if type Column and
/// therefore refer only to the persisted columns.
template <typename... PC>
struct WritingCursor<soa::Table<PC...>> {
  using persistent_table_t = soa::Table<PC...>;
  using cursor_t = decltype(std::declval<TableBuilder>().cursor<persistent_table_t>());

  template <typename... T>
  void operator()(T... args)
  {
    static_assert(sizeof...(PC) == sizeof...(T), "Argument number mismatch");
    ++mCount;
    cursor(0, extract(args)...);
  }

  /// Last index inserted in the table
  int64_t lastIndex()
  {
    return mCount;
  }

  bool resetCursor(TableBuilder& builder)
  {
    mBuilder = &builder;
    cursor = std::move(FFL(builder.cursor<persistent_table_t>()));
    mCount = -1;
    return true;
  }

  void setLabel(const char* label)
  {
    mBuilder->setLabel(label);
  }

  /// reserve @a size rows when filling, so that we do not
  /// spend time reallocating the buffers.
  void reserve(int64_t size)
  {
    mBuilder->reserve(typename persistent_table_t::column_types{}, size);
  }

  decltype(FFL(std::declval<cursor_t>())) cursor;

 private:
  template <typename T>
  static decltype(auto) extract(T const& arg)
  {
    if constexpr (soa::is_soa_iterator_t<T>::value) {
      return arg.globalIndex();
    } else {
      static_assert(!framework::has_type_v<T, framework::pack<PC...>>, "Argument type mismatch");
      return arg;
    }
  }

  /// The table builder which actually performs the
  /// construction of the table. We keep it around to be
  /// able to do all-columns methods like reserve.
  TableBuilder* mBuilder = nullptr;
  int64_t mCount = -1;
};

template <typename T>
struct Produces {
  static_assert(always_static_assert_v<T>, "Type must be a o2::soa::Table");
};

/// This helper class allows you to declare things which will be created by a
/// given analysis task. Notice how the actual cursor is implemented by the
/// means of the WritingCursor helper class, from which produces actually
/// derives.
template <typename... C>
struct Produces<soa::Table<C...>> : WritingCursor<typename soa::PackToTable<typename soa::Table<C...>::persistent_columns_t>::table> {
  using table_t = soa::Table<C...>;
  using metadata = typename aod::MetadataTrait<table_t>::metadata;

  // @return the associated OutputSpec
  OutputSpec const spec()
  {
    return OutputSpec{OutputLabel{metadata::tableLabel()}, metadata::origin(), metadata::description()};
  }

  OutputRef ref()
  {
    return OutputRef{metadata::tableLabel(), 0};
  }
};

/// Helper template for table transformations
template <typename METADATA>
struct TableTransform {
  using SOURCES = typename METADATA::sources;
  using ORIGINALS = typename METADATA::originals;

  using metadata = METADATA;
  using sources = SOURCES;

  constexpr auto sources_pack() const
  {
    return SOURCES{};
  }

  constexpr auto originals_pack() const
  {
    return ORIGINALS{};
  }

  template <typename Oi>
  constexpr auto base_spec() const
  {
    using o_metadata = typename aod::MetadataTrait<Oi>::metadata;
    return InputSpec{
      o_metadata::tableLabel(),
      header::DataOrigin{o_metadata::origin()},
      header::DataDescription{o_metadata::description()}};
  }

  template <typename... Os>
  std::vector<InputSpec> base_specs_impl(framework::pack<Os...>) const
  {
    return {base_spec<Os>()...};
  }

  std::vector<InputSpec> base_specs() const
  {
    return base_specs_impl(sources_pack());
  }

  constexpr auto spec() const
  {
    return OutputSpec{OutputLabel{METADATA::tableLabel()}, METADATA::origin(), METADATA::description()};
  }

  constexpr auto output() const
  {
    return Output{METADATA::origin(), METADATA::description()};
  }

  constexpr auto ref() const
  {
    return OutputRef{METADATA::tableLabel(), 0};
  }
};

/// This helper struct allows you to declare extended tables which should be
/// created by the task (as opposed to those pre-defined by data model)
template <typename T>
struct Spawns : TableTransform<typename aod::MetadataTrait<framework::pack_head_t<typename T::originals>>::metadata> {
  using extension_t = framework::pack_head_t<typename T::originals>;
  using base_table_t = typename aod::MetadataTrait<extension_t>::metadata::base_table_t;
  using expression_pack_t = typename aod::MetadataTrait<extension_t>::metadata::expression_pack_t;

  constexpr auto pack()
  {
    return expression_pack_t{};
  }

  T* operator->()
  {
    return table.get();
  }
  T& operator*()
  {
    return *table.get();
  }

  auto asArrowTable()
  {
    return extension->asArrowTable();
  }
  std::shared_ptr<typename T::table_t> table = nullptr;
  std::shared_ptr<extension_t> extension = nullptr;
};

/// Policy to control index building
/// Exclusive index: each entry in a row has a valid index
struct IndexExclusive {
  /// Generic builder for in index table
  template <typename... Cs, typename Key, typename T1, typename... T>
  static auto indexBuilder(const char* label, framework::pack<Cs...>, Key const&, std::tuple<T1, T...> tables)
  {
    static_assert(sizeof...(Cs) == sizeof...(T) + 1, "Number of columns does not coincide with number of supplied tables");
    using tables_t = framework::pack<T...>;
    using first_t = T1;
    auto tail = tuple_tail(tables);
    TableBuilder builder;
    auto cursor = framework::FFL(builder.cursor<o2::soa::Table<Cs...>>());

    std::array<int32_t, sizeof...(T)> values;
    iterator_tuple_t<std::decay_t<T>...> iterators = std::apply(
      [](auto&&... x) {
        return std::make_tuple(x.begin()...);
      },
      tail);

    using rest_it_t = decltype(pack_from_tuple(iterators));

    int32_t idx = -1;
    auto setValue = [&](auto& x) -> bool {
      using type = std::decay_t<decltype(x)>;
      constexpr auto position = framework::has_type_at_v<type>(rest_it_t{});

      lowerBound<Key>(idx, x);
      if (x == soa::RowViewSentinel{static_cast<uint64_t>(x.mMaxRow)}) {
        return false;
      } else if (x.template getId<Key>() != idx) {
        return false;
      } else {
        values[position] = x.globalIndex();
        ++x;
        return true;
      }
    };

    auto first = std::get<first_t>(tables);
    for (auto& row : first) {
      idx = row.template getId<Key>();

      if (std::apply(
            [](auto&... x) {
              return ((x == soa::RowViewSentinel{static_cast<uint64_t>(x.mMaxRow)}) && ...);
            },
            iterators)) {
        break;
      }

      auto result = std::apply(
        [&](auto&... x) {
          std::array<bool, sizeof...(T)> results{setValue(x)...};
          return (results[framework::has_type_at_v<std::decay_t<decltype(x)>>(rest_it_t{})] && ...);
        },
        iterators);

      if (result) {
        cursor(0, row.globalIndex(), values[framework::has_type_at_v<T>(tables_t{})]...);
      }
    }
    builder.setLabel(label);
    return builder.finalize();
  }
};
/// Sparse index: values in a row can be (-1), index table is isomorphic (joinable)
/// to T1
struct IndexSparse {
  template <typename... Cs, typename Key, typename T1, typename... T>
  static auto indexBuilder(const char* label, framework::pack<Cs...>, Key const&, std::tuple<T1, T...> tables)
  {
    static_assert(sizeof...(Cs) == sizeof...(T) + 1, "Number of columns does not coincide with number of supplied tables");
    using tables_t = framework::pack<T...>;
    using first_t = T1;
    auto tail = tuple_tail(tables);
    TableBuilder builder;
    auto cursor = framework::FFL(builder.cursor<o2::soa::Table<Cs...>>());

    std::array<int32_t, sizeof...(T)> values;

    iterator_tuple_t<std::decay_t<T>...> iterators = std::apply(
      [](auto&&... x) {
        return std::make_tuple(x.begin()...);
      },
      tail);

    using rest_it_t = decltype(pack_from_tuple(iterators));

    auto setValue = [&](auto& x, int idx) -> bool {
      using type = std::decay_t<decltype(x)>;
      constexpr auto position = framework::has_type_at_v<type>(rest_it_t{});

      if constexpr (std::is_same_v<framework::pack_element_t<position, framework::pack<std::decay_t<T>...>>, Key>) {
        values[position] = idx;
        return true;
      } else {
        lowerBound<Key>(idx, x);
        if (x == soa::RowViewSentinel{static_cast<uint64_t>(x.mMaxRow)}) {
          values[position] = -1;
          return false;
        } else if (x.template getId<Key>() != idx) {
          values[position] = -1;
          return false;
        } else {
          values[position] = x.globalIndex();
          ++x;
          return true;
        }
      }
    };

    auto first = std::get<first_t>(tables);
    for (auto& row : first) {
      auto idx = -1;
      if constexpr (std::is_same_v<first_t, Key>) {
        idx = row.globalIndex();
      } else {
        idx = row.template getId<Key>();
      }
      std::apply(
        [&](auto&... x) {
          (setValue(x, idx), ...);
        },
        iterators);

      cursor(0, row.globalIndex(), values[framework::has_type_at_v<T>(tables_t{})]...);
    }
    builder.setLabel(label);
    return builder.finalize();
  }
};

/// This helper struct allows you to declare index tables to be created in a task
template <typename T, typename IP = IndexSparse>
struct Builds : TableTransform<typename aod::MetadataTrait<T>::metadata> {
  using Key = typename T::indexing_t;
  using H = typename T::first_t;
  using Ts = typename T::rest_t;
  using index_pack_t = typename aod::MetadataTrait<T>::metadata::index_pack_t;

  T* operator->()
  {
    return table.get();
  }
  T& operator*()
  {
    return *table.get();
  }

  auto asArrowTable()
  {
    return table->asArrowTable();
  }
  std::shared_ptr<T> table = nullptr;

  constexpr auto pack()
  {
    return index_pack_t{};
  }

  template <typename... Cs, typename Key, typename T1, typename... Ts>
  auto build(framework::pack<Cs...>, Key const& key, std::tuple<T1, Ts...> tables)
  {
    this->table = std::make_shared<T>(IP::indexBuilder(aod::MetadataTrait<T>::metadata::tableLabel(), framework::pack<Cs...>{}, key, tables));
    return (this->table != nullptr);
  }
};

template <typename T>
using BuildsExclusive = Builds<T, IndexExclusive>;

/// This helper class allows you to declare things which will be created by a
/// given analysis task. Currently wrapped objects are limited to be TNamed
/// descendants. Objects will be written to a ROOT file at the end of the
/// workflow, in directories, corresponding to the task they were declared in.
/// Each object has associated handling policy, which is used by the framework
/// to determine the target file, e.g. analysis result, QA or control histogram,
/// etc.
template <typename T>
struct OutputObj {
  using obj_t = T;

  OutputObj(T&& t, OutputObjHandlingPolicy policy_ = OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType sourceType_ = OutputObjSourceType::OutputObjSource)
    : object(std::make_shared<T>(t)),
      label(t.GetName()),
      policy{policy_},
      sourceType{sourceType_},
      mTaskHash{0}
  {
  }

  OutputObj(std::string const& label_, OutputObjHandlingPolicy policy_ = OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType sourceType_ = OutputObjSourceType::OutputObjSource)
    : object(nullptr),
      label(label_),
      policy{policy_},
      sourceType{sourceType_},
      mTaskHash{0}
  {
  }

  void setObject(T const& t)
  {
    object = std::make_shared<T>(t);
    object->SetName(label.c_str());
  }

  void setObject(T&& t)
  {
    object = std::make_shared<T>(t);
    object->SetName(label.c_str());
  }

  void setObject(T* t)
  {
    object.reset(t);
    object->SetName(label.c_str());
  }

  void setObject(std::shared_ptr<T> t)
  {
    object = t;
    object->SetName(label.c_str());
  }

  void setHash(uint32_t hash)
  {
    mTaskHash = hash;
  }

  /// @return the associated OutputSpec
  OutputSpec const spec()
  {
    header::DataDescription desc{};
    auto lhash = compile_time_hash(label.c_str());
    std::memset(desc.str, '_', 16);
    std::stringstream s;
    s << std::hex << lhash;
    s << std::hex << mTaskHash;
    s << std::hex << reinterpret_cast<uint64_t>(this);
    std::memcpy(desc.str, s.str().c_str(), 12);
    return OutputSpec{OutputLabel{label}, "ATSK", desc, 0};
  }

  T* operator->()
  {
    return object.get();
  }

  T& operator*()
  {
    return *object.get();
  }

  OutputRef ref()
  {
    return OutputRef{std::string{label}, 0,
                     o2::header::Stack{OutputObjHeader{policy, sourceType, mTaskHash}}};
  }

  std::shared_ptr<T> object;
  std::string label;
  OutputObjHandlingPolicy policy;
  OutputObjSourceType sourceType;
  uint32_t mTaskHash;
};

/// This helper allows you to fetch a Sevice from the context or
/// by using some singleton. This hopefully will hide the Singleton and
/// We will be able to retrieve it in a more thread safe manner later on.
template <typename T>
struct Service {
  T* service;

  T* operator->() const
  {
    return service;
  }
};

template <typename T>
o2::soa::Filtered<T>* getTableFromFilter(const T& table, const expressions::Filter& filter)
{
  auto schema = table.asArrowTable()->schema();
  expressions::Operations ops = createOperations(filter);
  gandiva::NodePtr tree = nullptr;
  if (isSchemaCompatible(schema, ops)) {
    tree = createExpressionTree(ops, schema);
  } else {
    throw std::runtime_error("Partition filter does not match declared table type");
  }

  if constexpr (soa::is_soa_filtered_t<std::decay_t<T>>::value) {
    return new o2::soa::Filtered<T>{{table}, tree};
  } else {
    return new o2::soa::Filtered<T>{{table.asArrowTable()}, tree};
  }
}

template <typename T>
struct Partition {
  Partition(expressions::Node&& filter_) : filter{std::move(filter_)}
  {
  }

  void bindTable(T& table)
  {
    mFiltered.reset(getTableFromFilter(table, filter));
    bindExternalIndices(&table);
    getBoundToExternalIndices(table);
  }

  void setTable(const T& table)
  {
    mFiltered.reset(getTableFromFilter(table, filter));
  }

  template <typename... Ts>
  void bindExternalIndices(Ts*... tables)
  {
    if (mFiltered != nullptr) {
      mFiltered->bindExternalIndices(tables...);
    }
  }

  template <typename T2>
  void getBoundToExternalIndices(T2& table)
  {
    if (mFiltered != nullptr) {
      table.bindExternalIndices(mFiltered.get());
    }
  }

  void updatePlaceholders(InitContext& context)
  {
    expressions::updatePlaceholders(filter, context);
  }

  expressions::Filter filter;
  std::unique_ptr<o2::soa::Filtered<T>> mFiltered = nullptr;

  using iterator = typename o2::soa::Filtered<T>::iterator;
  using const_iterator = typename o2::soa::Filtered<T>::const_iterator;
  using filtered_iterator = typename o2::soa::Filtered<T>::iterator;
  using filtered_const_iterator = typename o2::soa::Filtered<T>::const_iterator;
  inline filtered_iterator begin()
  {
    return mFiltered->begin();
  }
  inline o2::soa::RowViewSentinel end()
  {
    return mFiltered->end();
  }
  inline filtered_const_iterator begin() const
  {
    return mFiltered->begin();
  }
  inline o2::soa::RowViewSentinel end() const
  {
    return mFiltered->end();
  }

  int64_t size() const
  {
    return mFiltered->size();
  }
};

template <typename P1, typename P2, typename G>
struct PairGenerator : public P1, public P2 {
  using PairIteratorType = std::tuple<typename G::iterator, typename A::iterator..., typename G::iterator, typename A::iterator...>;
  using GroupSlicer = AnalysisDataProcessorBuilder::GroupSlicer;

  struct PairIterator : public std::iterator<std::forward_iterator_tag, PairIteratorType>, public P1, public P2 {
   public:
    using reference = PairIteratorType&;
    using value_type = PairIteratorType;
    using pointer = PairIteratorType*;
    using iterator_category = std::forward_iterator_tag;

    PairIterator() = default;
    PairIterator(const P1& groupingPolicy, const G& grouping, const GroupSlicer& slicer) : P1(groupingPolicy), P2(), mGrouping(grouping), mSlicer(slicer) {}

    PairIterator(PairIterator const&) = default;
    PairIterator& operator=(PairIterator const&) = default;
    ~PairIterator() = default;

    void setAssociatedPolicy()
    {
      auto& [g1, g2] = P1::mCurrent;
      //P1::CombinationType currentGrouping = P1::mCurrent;
      auto it1 = mSlicer.begin();
      auto it2 = mSlicer.begin();
      for (auto& slice : mSlicer) {
        if (slice.groupingElement().index() == g1.index()) {
          it1 = slice;
          break;
        }
      }
      for (auto& slice : mSlicer) {
        if (slice.groupingElement().index() == g2.index()) {
          it2 = slice;
          break;
        }
      }
      auto associated1 = std::get<A...>(it1.associatedTables());
      associated1.bindExternalIndices(&mGrouping);
      auto associated2 = std::get<A...>(it2.associatedTables());
      associated2.bindExternalIndices(&mGrouping);

      P2::setTables(associated1, associated2);
    }

    // prefix increment
    PairIterator& operator++()
    {
      if (!P1::mIsEnd) {
        if (P2::mIsEnd) {
          P1::addOne();
          setAssociatedPolicy();
        } else {
          P2::addOne();
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
      auto& [g1, g2] = P1::mCurrent;
      auto& [a1, a2] = P2::mCurrent;
      return std::tuple(g1, a1, g2, a2);
    }
    bool operator==(const PairIterator& rh)
    {
      return (this->isEnd() && rh.isEnd()) || (*this == *rh);
    }
    bool operator!=(const PairIterator& rh)
    {
      return !(*this == rh);
    }
    bool isEnd()
    {
      return P1::mIsEnd && P2::mIsEnd;
    }

   private:
    GroupSlicer mSlicer;
    const G& mGrouping;
  };

  using iterator = typename PairIterator;
  using const_iterator = typename PairIterator;

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

  PairGenerator() = default;
  //PairGenerator(const P1& groupingPolicy, const P2& associatedPolicy) : mBegin(groupingPolicy, associatedPolicy), mEnd(groupingPolicy, associatedPolicy)
  //{
  //  mEnd.moveToEnd();
  //}
  ~PairGenerator() = default;

  void setPolicies(const P1& groupingPolicy, const P2& associatedPolicy, const G& grouping, const GroupSlicer& slicer)
  {
    mBegin = iterator(groupingPolicy, associatedPolicy, grouping, slicer);
    mEnd = iterator(groupingPolicy, associatedPolicy, grouping, slicer);
    mEnd.moveToEnd();
  }

 private:
  iterator mBegin;
  iterator mEnd;
};

template <std::string category, int catNeighbours, typename T1, typename G, typename... A>
struct Pair {
  using GroupingPolicy = o2::soa::CombinationsBlockStrictlyUpperSameIndexPolicy<G, G>;
  using AssociatedPolicy = o2::soa::CombinationsFullIndexPolicy<A...>;
  using GroupingGenerator = o2::soa::CombinationsGenerator<GroupingPolicy>;
  using AssociatedGenerator = o2::soa::CombinationsGenerator<>;
  using GroupSlicer = AnalysisDataProcessorBuilder::GroupSlicer;

  Pair(const T1& outsider){} : mGenerator(), mOutsider(outsider) {} //GroupingPolicy(category, catNeighbours, outsider), AssociatedPolicy()), mOutsider(outsider) {}

  void setTables(const G& grouping, const A&... associated)
  {
    grouping.bindExternalIndices(&associated...);
    auto associatedTuple = std::make_tuple(associated...);
    mGenerator.setPolicies(GroupingPolicy(category, catNeigbours, mOutsider, grouping, grouping), AssociatedPolicy(associated...), grouping, GroupSlicer(grouping, associatedTuple));
  }

 private:
  const T1 mOutsider;
  PairGenerator<GroupingPolicy, AssociatedPolicy, G> mGenerator;
};

} // namespace o2::framework

namespace o2::soa
{
/// On-the-fly adding of expression columns
template <typename T, typename... Cs>
auto Extend(T const& table)
{
  static_assert((soa::is_type_spawnable_v<Cs> && ...), "You can only extend a table with expression columns");
  using output_t = Join<T, soa::Table<Cs...>>;
  return output_t{{o2::framework::spawner(framework::pack<Cs...>{}, table.asArrowTable().get(), "dynamic extension"), table.asArrowTable()}, 0};
}

/// Template function to attach dynamic columns on-the-fly (e.g. inside
/// process() function). Dynamic columns need to be compatible with the table.
template <typename T, typename... Cs>
auto Attach(T const& table)
{
  static_assert((framework::is_base_of_template<o2::soa::DynamicColumn, Cs>::value && ...), "You can only attach dynamic columns");
  using output_t = Join<T, o2::soa::Table<Cs...>>;
  return output_t{{table.asArrowTable()}, table.offset()};
}
} // namespace o2::soa

#endif // o2_framework_AnalysisHelpers_H_DEFINED
