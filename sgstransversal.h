#pragma once

#include <assert.h>

#include "group.h"

namespace CoffeeCode {
	// compile time sgs generator
	template <size_t _Pivot, typename _Group>
	struct SGSGenerator
	{
		constexpr static const size_t Pivot = _Pivot;
		using GeneratingGroup = _Group;
		static_assert(1 <= Pivot && Pivot <= GeneratingGroup::Length);
	};


	// strong generating set
	template <typename... Generators>
	class SGSTransversal
	{
		// check all generators have the same type (i.e. length)?
		using are_generators_same_length = same_attribute<Generators::GeneratingGroup::Length...>;
		static_assert(are_generators_same_length::value);

		// type of StrongGeneratingSet
		using StrongGeneratingSetT = SGSTransversal<Generators...>;

		// number of generators
		constexpr static const size_t N_generators = sizeof...(Generators);
		static_assert(N_generators >= 1);

		// first generator
		using FirstGenerator = nth_element<0, Generators...>;

	public:
		constexpr static const size_t Length = FirstGenerator::GeneratingGroup::Length;

		// tuple type that the SGS acts on
		template<size_t Base>
		using _TupleT = NibbleStorageTypeArray<Base, Length>;

	private:
		// orbit iterator
		// this is a minimal implementation which satisfies the conditions
		// required for a range-based for. A proper input iterator
		// requires too much boilerplate for our purposes
		template <size_t Base>
		class CosetIterator
		{
		public:
			using IteratorT = CosetIterator<Base>;
			using TupleT = _TupleT<Base>;

		private:
			using TupleNonBaseIndexT = size_t;

			// fast compile-time stack
			size_t currentDepth;
			std::array<
				std::pair<TupleNonBaseIndexT, TupleT>,
				Base * Length
			> currentPath;

		public:
			CosetIterator() :
				currentDepth{ 0 },
				currentPath{ std::make_pair(0, TupleT{0}) },
				done{ false }
			{
			}

			// we want to prevent recursion from ever happening
			enum StepStatus { Unsuccessful, Successful };
			bool done;

			// this implements a depth-first search
			StepStatus Step()
			{
				// mutable reference
				// e.g. [2, {0, 2, 1, 0, 0, 0}] or [2, {0, 2, 2, 0, 0, 0}]
				auto&[firstNonBaseIndex, currentTuple] = currentPath[currentDepth];

				// CASE 1. if index pointer points past last element then move up
				if (firstNonBaseIndex >= Length) {
					// we've traversed the whole tree
					if (currentDepth == 0) {
						done = true;
						return Successful;
					}

					// otherwise go up
					else {
						currentDepth--;

						// get parent element and modify index
						auto&[firstNonBaseIndexParent, currentTupleParent] = currentPath[currentDepth];
						firstNonBaseIndexParent++;

						// try incrementing from there
						return Unsuccessful;
					}
				}


				// CASE 2. if element that index pointer points to cannot be incremented further,
				// move pointer right;
				// e.g. when {0, 2, 2, 0, 0, 0} in the second example
				if (currentTuple[firstNonBaseIndex] == Base - 1) {
					firstNonBaseIndex++;
					return Unsuccessful;
				}

				// CASE 3. increment digit at non-base index
				// e.g. to {0, 2, 2, 0, 0, 0} in the first example
				TupleT newTuple = currentTuple;
				newTuple[firstNonBaseIndex]++;

				/*** this is where we deviate from full DFS search ***/
				// CASE 3a. element is not canonically ordered: move right immediately without exploring further down
				if (!StrongGeneratingSetT::IsCanonical<Base>(newTuple)) {
					firstNonBaseIndex++;
					return Unsuccessful;
				}
				/*** end ***/
				// CASE 3b. element is canonically ordered: we will look at its children

				// put on stack
				assert(currentDepth < currentPath.size() - 1);
				currentDepth++;
				currentPath[currentDepth] = { firstNonBaseIndex, newTuple };

				return Successful;
			}


		public:
			CosetIterator& operator++()
			{
				while (Step() == Unsuccessful);

				return *this;
			}

			// inequality for end-checking
			struct Done {};
			bool operator!=(const Done&) const
			{
				return !done;
			}

			// dereference
			const TupleT& operator*() const
			{
				return currentPath[currentDepth].second;
			}
		};

	public:
		// range-based for loop iterator for tuples over Base
		template<size_t Base>
		inline static IteratorProxy<CosetIterator<Base>> TupleCosets()
		{
			return IteratorProxy<CosetIterator<Base>>();
		}

	private:
		// lexicographical order of tuples truncated to indices {0, ..., N-1}
		enum Order { LARGER, EQUAL, SMALLER };

		template<size_t N, typename T>
		inline static Order NumericalOrder(const T& lhs, const T& rhs)
		{
			assert(lhs.size() >= N && rhs.size() >= N);

			for (size_t i = 0; i < N; i++) {
				if (lhs[i] < rhs[i])
					return LARGER;
				if (lhs[i] > rhs[i])
					return SMALLER;
			}

			return EQUAL;
		}

		// can't fold since we rely on early exit, so use recursive template
		template<typename GeneratorX, typename... GeneratorXS, typename T, typename ListT = std::vector<T>>
		inline static bool IsCanonicalImpl(const T& tuple, const ListT&& todo)
		{
			using GeneratingGroup = typename GeneratorX::GeneratingGroup;
			constexpr const size_t Pivot = GeneratorX::Pivot;
			constexpr const bool RecursionEnd = sizeof...(GeneratorXS) == 0;

			ListT new_todo;

			for (const auto& parent : todo) {
				// iterate over orbit of tuple under Permutation
				for (const auto& child : GeneratingGroup::FastOrbit(parent)) {
					// check numerical order
					const auto order = NumericalOrder<Pivot>(tuple, child);

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
				};
			}

			// further generators to check?
			if constexpr (!RecursionEnd) {
				return IsCanonicalImpl<GeneratorXS...>(tuple, std::move(new_todo));
			}
			return true;
		}


	public:
		template<size_t Base>
		inline static bool IsCanonical(const _TupleT<Base>& tuple)
		{
			return IsCanonicalImpl<Generators...>(tuple, { tuple });
		}
	};
}