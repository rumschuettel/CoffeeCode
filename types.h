#pragma once

#include "ctmath.h"
#include "traits.h"

// clone from https://github.com/calccrypto/uint256_t.git
// provides fast 128 and 256 bit integers

//#include "src/vectorclass/vectorclass.h"

#include <boost/multiprecision/cpp_int.hpp>

namespace LibPopcount {
	#include "libpopcnt.h"
}

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

	inline size_t Popcount(const uint64_t& s)
	{
		return LibPopcount::popcount64(s);
	}
	template<typename T>
	inline size_t Popcount(const boost::multiprecision::number<T>& s)
	{
		const auto* limbs = s.backend().limbs();
		const size_t limb_count = s.backend().size();
		return LibPopcount::popcnt(limbs, limb_count * sizeof(limbs[0]));
	}

	template<typename StoreT>
	inline void OrBit(StoreT& s, const bool bit, const size_t i)
	{
		s |= static_cast<StoreT>(static_cast<StoreT>(bit) << i);
	}

	// storage for multiple limbs of bits in array form
	template<size_t Width, size_t Count>
	using BitStorageTypeArray = std::array< BitStorageType<Width>, Count >;



	// SIZE STORE
	// use these for storing elements up to some specific size
	template<size_t Base>
	using SizeStorageType = std::conditional_t<
		Base <= std::numeric_limits<uint8_t>::max(), uint8_t,
		std::conditional_t<
		Base <= std::numeric_limits<uint16_t>::max(), uint16_t,
		std::conditional_t<
		Base <= std::numeric_limits<uint32_t>::max(), uint32_t,
		std::conditional_t<
		Base <= std::numeric_limits<uint64_t>::max(), uint64_t,
		InvalidIntegerType<>
		>>>>;

	template<size_t Base, size_t Count>
	using SizeStorageTypeArray = std::array< SizeStorageType<Base>, Count >;


	// orbit sizes can be huge
	// for K_SYS, maximum orbit size is K_SYS!
	// meaning that maximum size of orbit has 2^n = K_SYS! bits.
	// It thus suffices to have an integer of size n >= log_2(K_SYS!) + 1 bits.
	constexpr static size_t MAX_GROUP_ORBIT_BITS = CoffeeCode::ilog2factorial(K_SYS) + 1;
	using OrbitType = BitStorageType<MAX_GROUP_ORBIT_BITS>;
		
	// multiplicity sizes are smaller, and upper bounded by Base^K_SYS
	template<size_t Base>
	using MultiplicityType = BitStorageType< ilog2(Base)*K_SYS + 1 >;

}