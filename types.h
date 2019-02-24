#pragma once

#include "ctmath.h"

//#include "src/vectorclass/vectorclass.h"

#include <boost/multiprecision/cpp_int.hpp>


#include <type_traits>
#include <array>


namespace CoffeeCode {
	// TODO: better integer type that can save at least Width number of bits


	// automatic selection of smallest-possible integer
	template<bool raise = true>
	struct InvalidIntegerType {
		static_assert(raise, "not enough bits");
	};

	// for storing Width number of bits
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


	template<typename StoreT>
	inline void OrBit(StoreT& s, const bool bit, const size_t i)
	{
		s |= static_cast<StoreT>(static_cast<StoreT>(bit) << i);
	}

	// tuple type
	template<size_t Width, size_t Length>
	using StdStoreExT = std::array< StdStoreT<Width>, Length >;

	template<size_t Base, size_t Length>
	using StdTupleStoreT = StdStoreExT< ilog2(Base - 1), Length >;
	// TODO: make this explicitly SIMD
}