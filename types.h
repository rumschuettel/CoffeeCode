#pragma once

//#include "src/vectorclass/vectorclass.h"
#include <type_traits>
#include <array>

namespace CoffeeCode {

	// automatic selection of smallest-possible integer
	template<bool raise = true>
	struct InvalidIntegerType {
		static_assert(raise, "not enough bits");
	};

	// for storing length bits
	template<size_t Width>
	using StdStoreT = typename std::conditional <
		Width <= 8,
		uint8_t,
		typename std::conditional <
		Width <= 16,
		uint16_t,
		typename std::conditional <
		Width <= 32,
		uint32_t,
		typename std::conditional <
		Width <= 64,
		uint64_t,
		InvalidIntegerType<true>
		>::type>::type>::type>::type;

	using StdBitT = StdStoreT<1>;



	// compile-time log2
	constexpr size_t log2(const size_t n)
	{
		if (n == 0) throw "log2(0) undefined";
		return ((n < 2) ? 1 : 1 + log2(n / 2));
	}

	// tuple type
	template<size_t Width, size_t Length>
	using StdStoreExT = std::array< StdStoreT<Width>, Length >;

	template<size_t Base, size_t Length>
	using StdTupleStoreT = StdStoreExT< log2(Base - 1), Length >;
	// TODO: make this explicitly SIMD
}