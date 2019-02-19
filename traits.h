#pragma once

#include <cstddef>
using std::size_t;

namespace {
	// template traits

	// check that elements have same attribute
	template <auto T, auto... Ts>
	struct same_attribute : std::bool_constant< ((T == Ts) && ...) > {};

	// check that elements appear exactly once
	template <auto T, auto... Ts>
	struct all_unequal : std::bool_constant< all_unequal<Ts...>::value && ((T != Ts) && ...) > {};
	template <auto T>
	struct all_unequal<T> : std::true_type {};

	// get nth element
	template <size_t Idx, typename... Ts>
	using nth_element = typename std::tuple_element<Idx, std::tuple<Ts...>>::type;
}