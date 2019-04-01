#pragma once

#include <cstddef>
#include <type_traits>
#include <tuple>
#include <variant>

using std::size_t;



// template traits

// check that elements have same attribute
template<auto T, auto... Ts>
struct same_attribute : std::bool_constant< ((T == Ts) && ...) > {};

// check that elements appear exactly once
template<auto T, auto... Ts>
struct all_unequal : std::bool_constant< all_unequal<Ts...>::value && ((T != Ts) && ...) > {};
template<auto T>
struct all_unequal<T> : std::true_type {};

// get nth element
template<size_t Idx, typename... Ts>
using nth_element = typename std::tuple_element<Idx, std::tuple<Ts...>>::type;

// check whether something is a tuple
template<typename T>
constexpr bool is_pair = false;
template<typename... Ts>
constexpr bool is_pair<std::pair<Ts...>> = true;

// compile time MSVC sum bug workaround
template<auto... Args>
constexpr auto sum = (Args + ...);


#if __has_include(<experimental/type_traits>)
#include <experimental/type_traits>
template<typename Default, template<class...> class Op, class... Args>
using detected_or_t = std::experimental::detected_or_t<Default, Op, Args...>;
#else
#include <boost/type_traits/detected_or.hpp>
template<typename Default, template<class...> class Op, class... Args>
using detected_or_t = boost::detected_or_t<Default, Op, Args...>;
#endif


// iterator proxy type

// features a begin and end function to be used in a
// range-based for loop
template <typename IteratorT, typename DoneT = typename IteratorT::DoneT>
struct StatefulIteratorProxy {

	IteratorT begin() const { return IteratorT(); }
	DoneT end() const { return DoneT(); }
};
