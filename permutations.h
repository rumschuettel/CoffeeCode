#pragma once

#include "types.h"

#include <iterator>

namespace CoffeeCode {
	// compile time integer sequence
	template <typename IndexT, IndexT... List>
	struct PermutationList
	{
		typedef IndexT IndexT;

		constexpr static const size_t N_indices = sizeof...(List);
		constexpr static const IndexT list[N_indices] = { List... };

		template <typename T>
		inline static std::array<T, N_indices> permute(const std::array<T, N_indices>& src) noexcept
		{
			// rely on return value optimization
			return { src[List]... };
		}
	};


	// compile time sgs generator
	template <typename IndexT, IndexT Index, IndexT... List>
	struct SGSGenerator
	{
		constexpr static const IndexT Index = Index;
		using Permutation = PermutationList<IndexT, List...>;
	};


	// strong generating set
	template <typename... Generators>
	class StrongGeneratingSet
	{
		// all generators have the same type (i.e. length)?
		template <class T, class... Ts>
		struct are_same : std::bool_constant<(std::is_same_v<typename T::Permutation, typename Ts::Permutation> && ...)> {};
		using generators_same_type = are_same<Generators...>;

		// number of generators
		constexpr static const size_t N_generators = sizeof...(Generators);
		static_assert(N_generators > 0);

		// generator length
		using FirstGenerator = typename std::tuple_element<0, std::tuple<Generators...>>::type;
		constexpr static const size_t N_indices = FirstGenerator::N_indices;


		// orbit iterator
		template <size_t base>
		class OrbitIterator
		{
		public:
			using iterator_category = std::forward_iterator_tag;

			// store N_indices many values of size base
			using value_type = StdStoreExT<base, N_indices>;
			using pointer = value_type*;
			using reference = value_type&;

			// note that we never need the difference_type, so we don't care whether it's representable or not
			using difference_type = std::ptrdiff_t; 

			// block length (bases are not 0-indexed, but start at 1, so subtract 1)
			constexpr static const size_t block_length = log2(base-1);
		
			OrbitIterator()
			{
				std::cout << "block length " << block_length;
			}

		};

	public:
		StrongGeneratingSet()
		{
			((std::cout << Generators::list[1] << "\n"), ...);

			OrbitIterator<4> iterator;
		}
	};

	void test() {
		using IndexT = StdStoreT<32>;
		using p1 = PermutationList<IndexT, 0, 1, 2, 3, 4, 5, 6, 7>;
		using p2 = PermutationList<IndexT, 0, 7, 1, 6, 2, 5, 3, 4>;

		Bit2StoreT<8> vec{ 10, 20, 30, 40, 50, 60, 70, 80 };
		for (const auto e : vec)
			std::cout << (int)e << " ";
		std::cout << "\n";

		auto vec2 = p2::permute(vec);
		for (const auto e : vec2)
			std::cout << (int)e << " ";
		std::cout << "\n";

		StrongGeneratingSet<p1, p2> group;
	}
}