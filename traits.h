#pragma once

#include <iostream>

namespace {
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



	// compile time integer exponent
	constexpr size_t ipow(const size_t base, const int exp, const size_t result = 1) {
		return exp < 1 ? result : ipow(base*base, exp / 2, (exp % 2) ? result * base : result);
	}
	// compile time size of base k tuple
	template<const size_t base, const size_t tuple_length>
	struct BaseKSubsets {
		static constexpr auto count = ipow(base, tuple_length);
	};
	// compile time bitmask with k 1s
	template<const size_t number_of_1s>
	struct Bitmask1s {
		static constexpr auto mask = (1ull << number_of_1s) - 1;
	};


	// debug helper
	template<typename T>
	void print(const T& vec)
	{
		for (const auto e : vec)
			std::cout << (int)e << " ";
		std::cout << "\n";

	}
}