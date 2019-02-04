#pragma once

#include "types.h"
#include "traits.h"



namespace CoffeeCode {
	// compile time permutation
	// we don't check whether the indices given are a valid permutation
	template <size_t... Indices>
	struct Permutation
	{
		// length
		constexpr static const size_t Length = sizeof...(Indices);
		static_assert(Length >= 1);

	private:
		// check all indices are smaller than Length
		using are_indices_within_range = std::bool_constant<((Indices < Length) && ...)>;
		static_assert(are_indices_within_range::value);
		using are_indices_all_unequal = all_unequal<Indices...>;
		static_assert(are_indices_all_unequal::value);

	public:
		// keep around as list
		using PermutationIndicesT = std::array<size_t, Length>;
		constexpr static const PermutationIndicesT IndicesList = { Indices... };

		// acting on a tuple
		// the template folding checks that the tuple is long enough, but not whether it is too long
		template<typename T>
		constexpr inline static T permute(const T& src) noexcept
		{
			// rely on return value optimization
			return { src[Indices]... };
		}
	};
}