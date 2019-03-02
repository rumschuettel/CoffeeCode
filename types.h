#pragma once

#include "ctmath.h"

//#include "src/vectorclass/vectorclass.h"

#include <boost/multiprecision/cpp_int.hpp>

namespace LibPopcount {
	#include "libpopcnt.h"
}

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
	inline size_t Popcount(const StoreT& s)
	{
		return LibPopcount::popcount64(s);
	}


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


	// orbit sizes can be huge
	// for K_SYS, maximum orbit size is K_SYS!
	// meaning that maximum size of orbit has 2^n = K_SYS! bits.
	// It thus suffices to have an integer of size n >= log_2(K_SYS!) + 1 bits.
	template<size_t BaseSize = K_SYS>
	struct OrbitType {
		constexpr static size_t MAX_GROUP_ORBIT = CoffeeCode::ilog2factorial(BaseSize) + 1;
		using SizeT = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<
			MAX_GROUP_ORBIT,
			MAX_GROUP_ORBIT,
			boost::multiprecision::unsigned_magnitude,
			boost::multiprecision::unchecked
		>>;

		constexpr inline static SizeT Factorial(size_t n)
		{
			return n <= 1 ? 1 : (Factorial(n-1) * n);
		}

		constexpr inline static SizeT Binomial(size_t n, size_t k)
		{
			return (k == 0 || k == n) ? 1 : Binomial(n-1, k-1) + Binomial(n-1, k);
		}
	};
}