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
	template<size_t length>
	using StdStoreT = typename std::conditional <
		length <= 8,
		uint8_t,
		typename std::conditional <
		length <= 16,
		uint16_t,
		typename std::conditional <
		length <= 32,
		uint32_t,
		typename std::conditional <
		length <= 64,
		uint64_t,
		InvalidIntegerType<true>
		>::type>::type>::type>::type;

	using StdBitT = StdStoreT<1>;


	// for storing length numbers in bit-blocks
	constexpr size_t log2(const size_t n)
	{
		if (n == 0) throw "log2(0) undefined";
		return ((n < 2) ? 1 : 1 + log2(n / 2));
	}

	template<size_t base, size_t length>
	using StdStoreExT = StdStoreT< log2(base-1) * length >;


	// for symmetry calculations we store pairs of bits in bytes
	template<size_t Length>
	using Bit2StoreT = std::array< StdStoreT<2>, Length >;
	// TODO: make this explicitly SIMD
}