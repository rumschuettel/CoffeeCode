#pragma once

#include <cstddef>
#include <type_traits>
#include <tuple>

using std::size_t;



// template traits

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

// cast-and-divide
template<typename T, auto Enumerator, auto Denominator>
constexpr auto cdiv = static_cast<T>(Enumerator) / static_cast<T>(Denominator);

#include <boost/type_traits/detected_or.hpp>
template<typename Default, template<class...> class Op, class... Args>
using detected_or_t = boost::detected_or_t<Default, Op, Args...>;

// iterator proxy type

// features a begin and end function to be used in a
// range-based for loop
template <typename IteratorT, typename DoneT = typename IteratorT::DoneT>
struct StatefulIteratorProxy {

	IteratorT begin() const { return IteratorT(); }
	DoneT end() const { return DoneT(); }
};
