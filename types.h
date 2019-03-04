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
	// TODO: better integer type that can save at least Width number of bits
	template<size_t Width>
	using WideIntegerT = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
		Width,
		Width,
		boost::multiprecision::unsigned_magnitude,
		boost::multiprecision::unchecked
	>>;

	// automatic selection of smallest-possible integer
	template<bool raise = true>
	struct InvalidIntegerType {
		static_assert(raise, "not enough bits");
	};

	// for storing Width number of bits
	template<size_t Width>
	using StdStoreT = std::conditional_t<
		Width <= 8,	uint8_t,
		std::conditional_t<
		Width <= 16, uint16_t,
		std::conditional_t<
		Width <= 32, uint32_t,
		std::conditional_t<
		Width <= 64, uint64_t,
		std::conditional_t<
		Width <= 128, boost::multiprecision::uint128_t,
		std::conditional_t<
		Width <= 256, boost::multiprecision::uint256_t,
		std::conditional_t<
		Width <= 512, boost::multiprecision::uint512_t,
		InvalidIntegerType<>
	>>>>>>>;

	using StdBitT = StdStoreT<1>;


	inline size_t Popcount(const uint64_t& s)
	{
		return LibPopcount::popcount64(s);
	}
	// TODO: wider popcounts


	template<typename StoreT>
	inline void OrBit(StoreT& s, const bool bit, const size_t i)
	{
		s |= static_cast<StoreT>(static_cast<StoreT>(bit) << i);
	}

	// tuple type
	template<size_t Width, size_t Length>
	using StdStoreExT = std::array< StdStoreT<Width>, Length >;

	template<size_t Base, size_t Length>
	using StdTupleStoreT = StdStoreExT< ilog2(Base), Length >;
	// TODO: make this explicitly SIMD


	// orbit sizes can be huge
	// for K_SYS, maximum orbit size is K_SYS!
	// meaning that maximum size of orbit has 2^n = K_SYS! bits.
	// It thus suffices to have an integer of size n >= log_2(K_SYS!) + 1 bits.
	constexpr static size_t MAX_GROUP_ORBIT_BITS = CoffeeCode::ilog2factorial(K_SYS) + 1;
	using OrbitType = StdStoreT<MAX_GROUP_ORBIT_BITS>;
		
	// multiplicity sizes are smaller, and upper bounded by Base^K_SYS
	template<size_t Base>
	using MultiplicityType = StdStoreT< ilog2(Base)*K_SYS + 1 >;

}