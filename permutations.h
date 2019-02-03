#pragma once

#include "types.h"

#include <iterator>
#include <numeric>


namespace CoffeeCode {
	// compile time permutation
	// we don't check whether the indices given are a valid permutation
	template <size_t... Indices>
	struct Permutation
	{
		// length
		constexpr static const size_t Length = sizeof...(Indices);

	private:
		// check all generators have the same type (i.e. length)?
		using indices_valid = std::bool_constant<((Indices <= Length) && ...)>;

	public:
		// keep around as list
		using PermutationListT = std::array<size_t, Length>;
		constexpr static const PermutationListT IndicesList = { Indices... };

		// acting on a tuple
		// the template folding checks that the tuple is long enough, but not whether it is too long
		template<typename T>
		constexpr inline static T permute(const T& src) noexcept
		{
			// rely on return value optimization
			return { src[Indices]... };
		}
	};


	// compile time sgs generator
	template <size_t Pivot, typename _Permutation>
	struct SGSGenerator
	{
		constexpr static const size_t Pivot = Pivot;
		using Permutation = _Permutation;
	};


	// strong generating set
	template <typename... Generators>
	class StrongGeneratingSet
	{
		// check all generators have the same type (i.e. length)?
		template <class T, class... Ts>
		struct same_length_generators : std::bool_constant< ((T::Permutation::Length == Ts::Permutation::Length) && ...) > {};
		using same_length_generators_true = same_length_generators<Generators...>;

		// type of StrongGeneratingSet
		using StrongGeneratingSetT = StrongGeneratingSet<Generators...>;

		// number of generators
		constexpr static const size_t N_generators = sizeof...(Generators);
		static_assert(N_generators > 0);

		// generator length
		using FirstGenerator = typename std::tuple_element<0, std::tuple<Generators...>>::type;
		constexpr static const size_t Length = FirstGenerator::Permutation::Length;

	public:
		// tuple type that the SGS acts on
		template<size_t Base>
		using TupleT = StdStoreExT< log2(Base-1), Length >;

		// orbit iterator
		// this is a minimal implementation which satisfies the conditions
		// required for a range-based for. A proper input iterator
		// requires too much boilerplate for our purposes
		template <size_t Base>
		class OrbitIterator
		{
		public:
			using IteratorT = OrbitIterator<Base>;
			using TupleT = TupleT<Base>;

		private:
			using TupleNonBaseIndexT = size_t;

			// fast compile-time stack
			size_t currentDepth;
			std::array<
				std::pair<TupleNonBaseIndexT, TupleT>,
				Base * Length
			> currentPath;

		public:				
			OrbitIterator() :
				currentDepth{ 0 },
				currentPath{ std::make_pair(0, TupleT{0}) }
			{
				// we increment once initially;
				// this is so that the final state is the all zeros tuple
				// note that the all-zeros one will not be listed
				++(*this);
			}

			enum StepStatus { CallAgain, Done };

			// this implements a depth-first search
			StepStatus Step()
			{
				// mutable reference
				// e.g. [2, {0, 2, 1, 0, 0, 0}] or [2, {0, 2, 2, 0, 0, 0}]
				auto& [firstNonBaseIndex, currentTuple] = currentPath[currentDepth];

				// CASE 1. if index pointer points past last element then move up
				if (firstNonBaseIndex >= Length) {
					// we've traversed the whole tree
					if (currentDepth == 0) {
						return Done;
					}

					// otherwise go up
					else {
						currentDepth--;

						// get parent element and modify index
						auto& [firstNonBaseIndexParent, currentTupleParent] = currentPath[currentDepth];
						firstNonBaseIndexParent++;

						// try incrementing from there
						return CallAgain;
					}
				}


				// CASE 2. if element that index pointer points to cannot be incremented further,
				// move pointer right;
				// e.g. when {0, 2, 2, 0, 0, 0} in the second example
				if (currentTuple[firstNonBaseIndex] == Base-1) {
					firstNonBaseIndex++;
					return CallAgain;
				}

				// CASE 3. increment digit at non-base index
				// e.g. to {0, 2, 2, 0, 0, 0} in the first example
				TupleT newTuple = currentTuple;
				newTuple[firstNonBaseIndex]++;

				/*** this is where we deviate from full DFS search ***/
				// CASE 3a. element is not canonically ordered: move right immediately without exploring further down
				if (!StrongGeneratingSetT::IsCanonical<Base>(newTuple)) {
					firstNonBaseIndex++;
					return CallAgain;
				}
				/*** end ***/
				// CASE 3b. element is canonically ordered: we will look at its children

				// put on stack
				assert(currentDepth < currentPath.size() - 1);
				currentDepth++;
				currentPath[currentDepth] = { firstNonBaseIndex, newTuple };

				return Done;
			}


		public:
			OrbitIterator& operator++()
			{
				do {
				} while (Step() == CallAgain);

				return *this;
			}

			// inequality
			bool operator!=(const TupleT& end) const
			{
				return currentPath[currentDepth].second != end;
			}
			// dereference
			const TupleT& operator*() const
			{
				return currentPath[currentDepth].second;
			}
		};

		// features a begin and end function to be used in a
		// range-based for loop
		template <size_t Base>
		struct OrbitIteratorProxy {
			using IteratorT = OrbitIterator<Base>;
			using TupleT = typename IteratorT::TupleT;

			IteratorT begin() const { return IteratorT(); }
			const TupleT end() const { return TupleT{ 0 }; }
		};

		// range-based for loop iterator for tuples over Base
		template<size_t Base>
		inline static OrbitIteratorProxy<Base> TupleCosets()
		{
			return OrbitIteratorProxy<Base>();
		}

	private:
		// lexicographical order of tuples up to index n
		enum Order { LARGER, EQUAL, SMALLER };

		template<size_t N, typename T>
		inline static Order NumericalOrder(const T& lhs, const T& rhs)
		{
			assert(lhs.size() > N && rhs.size() > N);

			for (size_t i = 0; i < N; i++) {
				if (lhs[i] < rhs[i])
					return LARGER;
				if (lhs[i] > rhs[i])
					return SMALLER;
			}

			return EQUAL;				
		}

		// can't fold since we rely on early exit, so use recursive template
		template<typename GeneratorX, typename... GeneratorXS, typename T>
		inline static bool IsCanonicalImpl(const T& todo)
		{
			using Permutation = typename GeneratorX::Permutation;
			constexpr const size_t Pivot = GeneratorX::Pivot;
			constexpr const bool RecursionEnd = sizeof...(GeneratorXS) > 0;

			T new_todo;

			for (const auto& parent : todo) {

				// iterate over orbit of tuple under Permutation
				auto child = parent;
				do {
					// check numerical order
					const auto order = NumericalOrder<Pivot>(parent, child);

					// return false if new tuple
					// lexicographically larger up to the pivot point
					if (order == LARGER)
						return false;

					// equal up till pivot point? need to check
					// further generators
					if constexpr (!RecursionEnd) {
						if (order == EQUAL)
							new_todo.push_back(child);
					}

					// if smaller we don't gain information; ignore

					// permute to get new orbit element
					child = Permutation::permute(child);

				} while (child != parent);
			}

			// further generators to check?
			if constexpr (!RecursionEnd) {
				return IsCanonicalImpl<GeneratorXS>(new_todo);
			}
			return true;
		}


	public:
		template<size_t Base>
		inline static bool IsCanonical(const TupleT<Base>& tuple)
		{
			std::vector<TupleT<Base>> todo{ tuple };

			return IsCanonicalImpl<Generators...>(todo);
		}
	};

	void test() {
		using p1 = Permutation<0, 3, 2, 4, 6, 5, 1, 7>;
		using p2 = Permutation<0, 7, 1, 6, 2, 5, 3, 4>;
		using s1 = SGSGenerator<5, p1>;
		using s2 = SGSGenerator<7, p2>;
		using sgs = StrongGeneratingSet<s1, s2>;
		
		
		using TupleT = sgs::TupleT<100>;
		TupleT vec{ 1, 2, 3, 4, 5, 6, 7, 8 };
		for (const auto e : vec)
			std::cout << (int)e << " ";
		std::cout << "\n";

		TupleT vec2 = p2::permute(vec);
		for (const auto e : vec2)
			std::cout << (int)e << " ";
		std::cout << "\n\n";

		std::cout << std::boolalpha << sgs::IsCanonical<100>(vec) << "\n";

		for (const auto& tuple : sgs::TupleCosets<3>()) {
			for (const auto e : tuple)
				std::cout << (int)e << " ";
			std::cout << "\n";
		}
	}
}