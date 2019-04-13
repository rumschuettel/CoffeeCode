#pragma once

#include "ctmath.h"
#include "traits.h"

#include <array>

namespace CoffeeCode {

	// automatic selection of smallest-possible integer
	template<bool raise = true>
	struct InvalidIntegerType {
		static_assert(raise, "not enough bits");
	};

	// BIT STORES
	// we use these types for bit manipulations, i.e. binary vectors, binary matrices, masks, colorings

	// for storing Width number of bits
	// TODO: add wider-than-64 bit bit store that allows constexpr operations
	template<size_t Width>
	using BitStorageType = std::conditional_t<
		Width <= 8,	uint8_t,
		std::conditional_t<
		Width <= 16, uint16_t,
		std::conditional_t<
		Width <= 32, uint32_t,
		std::conditional_t<
		Width <= 64, uint64_t,
		InvalidIntegerType<>
	>>>>;

	using BitType = BitStorageType<1>;

	// storage for multiple limbs of bits in array form
	template<size_t Width, size_t Count>
	using BitStorageTypeArray = std::array< BitStorageType<Width>, Count >;

	// storage for multiple nibbles from 0...Base-1
	template<size_t Base, size_t Count>
	using NibbleStorageTypeArray = std::array< BitStorageType<ilog2(Base-1)>, Count >;


	// SIZE STORE
	// use these for storing elements up to some specific size
	// Note: Log2Base for 256 is 8, so if we actually want to store 256 we need 9 bits.
	template<size_t Base, size_t Log2Base = ilog2(Base)>
	using SizeStorageType = std::conditional_t <
		Log2Base < 64, uint64_t,
#ifdef FLOATING_POINT_MULTIPLICITY
		double
#else
		std::conditional_t<
		Log2Base < 128, uint128_t,
		std::conditional_t<
		Log2Base < 256, uint256_t,
		std::conditional_t<
		Log2Base < 512, uint512_t,
		InvalidIntegerType<>
		>>>
#endif
	>;

	// orbit sizes can be huge
	// for K_SYS, maximum orbit size is K_SYS!
	// meaning that maximum size of orbit has 2^n = K_SYS! bits.
	// It thus suffices to have an integer of size n >= log_2(K_SYS!) + 1 bits.
	// SizeStorageType allows overriding the bit count
	constexpr static size_t LOG2_MAX_GROUP_ORBIT_SIZE = CoffeeCode::ilog2factorial(K_SYS);
	using OrbitType = SizeStorageType<0, LOG2_MAX_GROUP_ORBIT_SIZE>;
		
	// multiplicity sizes are smaller, and upper bounded by Base^K_SYS
	template<size_t Base>
	using MultiplicityType = SizeStorageType<0, ilog2(Base)*K_SYS>;


	///////////////////
	template<size_t Colors>
	using SystemColoringType = NibbleStorageTypeArray<Colors, K_SYS>;
	template<size_t Colors>
	using FullColoringType = NibbleStorageTypeArray<Colors, K_SYS + K_ENV>;


	///////////////////
	template<size_t Colors>
	using CosetGeneratorCallbackType = std::function<void(const SystemColoringType<Colors>&, const MultiplicityType<Colors>, const size_t)>;

}