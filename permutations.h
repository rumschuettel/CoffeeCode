#pragma once

#include "types.h"

#include <iterator>
#include <numeric>
#include <set>

namespace {
	// debug helper
	template<typename T>
	void print(const T& vec)
	{
		for (const auto e : vec)
			std::cout << (int)e << " ";
		std::cout << "\n";
	}
}

namespace CoffeeCode {
	// check that elements have same attribute
	template <auto T, auto... Ts>
	struct same_attribute : std::bool_constant< ((T == Ts) && ...) > {};

	// check that elements appear exactly once
	template <auto T, auto... Ts>
	struct all_unequal : std::bool_constant< all_unequal<Ts...>::value && ((T != Ts) && ...) > {};
	template <auto T>
	struct all_unequal<T> : std::true_type {};

	// get nth element
	template <size_t Idx, typename... Ts>
	using nth_element = typename std::tuple_element<Idx, std::tuple<Ts...>>::type;

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

	// compile time permutation group
	template <typename... Permutations>
	struct Group
	{
		// number of generators and copy of generators as list
		constexpr static const size_t GeneratorCount = sizeof...(Permutations);
		static_assert(GeneratorCount >= 1);

		// group length
		using FirstPermutation = nth_element<0, Permutations...>;
		constexpr static const size_t Length = FirstPermutation::Length;

	private:
		// check all generators have the same type (i.e. length)?
		using are_generators_same_length = same_attribute<Permutations::Length...>;
		static_assert(are_generators_same_length::value);


		// orbit iterator
		template<typename T>
		class OrbitIterator {
			// orbit of a tuple
			using OrbitSetT = std::set<T>; // not unordered_set, as we can actually compare elements
			using OrbitListT = std::vector<T>;

			OrbitSetT orbit;
			OrbitListT new_orbit_elements;
			size_t next_generator_to_check;
			T current_element;
			T current_parent;
			bool done;

		private:
			enum { Unsuccessful, Successful };

			template<typename Generator>
			bool CheckAndAdd()
			{
				const T element = Generator::permute(current_parent);
				// trivial check
				if (element == current_parent) return Unsuccessful;
				// if element actually changed, try to add to set
				auto[_, inserted] = orbit.insert(element);
				// if it was not in the set yet, we found a new one
				if (inserted) {
					new_orbit_elements.push_back(element);
					current_element = element;
				}
				return inserted ? Successful : Unsuccessful;
			}


			template<size_t GeneratorIdx = 0>
			inline bool Step() {
				// trivial case of one generator
				if constexpr (GeneratorCount == 1) {
					const auto new_element = FirstPermutation::permute(current_element);
					if (new_element == current_parent) {
						// done; cycled once
						done = true;
						return Successful;
					}
					
					// otherwise we found a new element
					current_element = new_element;
					return Successful;
				}

				// multiple generators
				else if constexpr (GeneratorIdx == GeneratorCount) {
					// we are completely done
					if (new_orbit_elements.size() == 0) {
						done = true;
						return Successful;
					}

					// otherwise set new parent element
					current_parent = new_orbit_elements.back();
					new_orbit_elements.pop_back();
					next_generator_to_check = 0;
					return Unsuccessful;
				}
				else {
					if (GeneratorIdx < next_generator_to_check)
						return Step<GeneratorIdx + 1>();

					// now we're at a valid generator that we can apply to the current_parent
					next_generator_to_check = GeneratorIdx + 1;
					using Generator = nth_element<GeneratorIdx, Permutations...>;
					return CheckAndAdd<Generator>();
				}
			}

		public:
			OrbitIterator(const T& start) :
				orbit{ start },
				next_generator_to_check{ 0 },
				current_element{ start },
				current_parent{ start },
				done{ false }
			{
			}

			OrbitIterator& operator++()
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
			const T& operator*() const
			{
				return current_element;
			}
		};

	public:

		template<typename T>
		struct OrbitIteratorProxy {
		private:
			const T tuple;

		public:
			OrbitIteratorProxy(const T& tuple)
				: tuple(tuple)
			{}

			using IteratorT = OrbitIterator<T>;
			IteratorT begin() const { return IteratorT(tuple); }
			const auto end() const { return IteratorT::Done(); }
		};

		template<typename T>
		inline static OrbitIteratorProxy<T> FastOrbit(const T& tuple)
		{
			return OrbitIteratorProxy(tuple);
		}

		// orbit of a tuple as a set
		template<typename T>
		using OrbitSetT = std::set<T>; // not unordered_set, as we can actually compare elements
		template<typename T>
		using OrbitListT = std::vector<T>;

		template<typename T>
		inline static OrbitSetT<T> Orbit(const T& tuple)
		{
			OrbitSetT<T> orbit{ tuple };
			OrbitListT<T> new_orbit_elements{ tuple };

			auto check_and_add = [&](const T& original, const T&& element) {
				if (original == element) return;
				auto[_, inserted] = orbit.insert(element);
				if (inserted) new_orbit_elements.push_back(element);
			};

			while (new_orbit_elements.size() > 0) {
				// we move the new elements to a separate list to iterate over,
				// and free the original one so that we can fill it up anew
				auto new_orbit_elements_copy = std::move(new_orbit_elements);
				new_orbit_elements.clear();

				// now iterate over all new orbit elements previously found,
				// and apply each generator in turn
				for (const auto& element : new_orbit_elements_copy)
					(check_and_add(element, Permutations::permute(element)), ...);
			}

			return orbit;
		}
	};


	// compile time sgs generator
	template <size_t Pivot, typename _Group>
	struct SGSGenerator
	{
		constexpr static const size_t Pivot = Pivot;
		using GeneratingGroup = _Group;
		static_assert(1 <= Pivot && Pivot <= GeneratingGroup::Length);
	};


	// strong generating set
	template <typename... Generators>
	class StrongGeneratingSet
	{
		// check all generators have the same type (i.e. length)?
		using are_generators_same_length = same_attribute<Generators::GeneratingGroup::Length...>;
		static_assert(are_generators_same_length::value);

		// type of StrongGeneratingSet
		using StrongGeneratingSetT = StrongGeneratingSet<Generators...>;

		// number of generators
		constexpr static const size_t N_generators = sizeof...(Generators);
		static_assert(N_generators >= 1);

		// generator length
		using FirstGenerator = nth_element<0, Generators...>;
		constexpr static const size_t Length = FirstGenerator::GeneratingGroup::Length;

	public:
		// tuple type that the SGS acts on
		template<size_t Base>
		using TupleT = StdStoreExT< log2(Base-1), Length >;

		// orbit iterator
		// this is a minimal implementation which satisfies the conditions
		// required for a range-based for. A proper input iterator
		// requires too much boilerplate for our purposes
		template <size_t Base>
		class CosetIterator
		{
		public:
			using IteratorT = CosetIterator<Base>;
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
				auto& [firstNonBaseIndex, currentTuple] = currentPath[currentDepth];

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
						auto& [firstNonBaseIndexParent, currentTupleParent] = currentPath[currentDepth];
						firstNonBaseIndexParent++;

						// try incrementing from there
						return Unsuccessful;
					}
				}


				// CASE 2. if element that index pointer points to cannot be incremented further,
				// move pointer right;
				// e.g. when {0, 2, 2, 0, 0, 0} in the second example
				if (currentTuple[firstNonBaseIndex] == Base-1) {
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

		// features a begin and end function to be used in a
		// range-based for loop
		template <size_t Base>
		struct OrbitIteratorProxy {
			using IteratorT = CosetIterator<Base>;
			IteratorT begin() const { return IteratorT(); }
			const auto end() const { return IteratorT::Done(); }
		};

		// range-based for loop iterator for tuples over Base
		template<size_t Base>
		inline static OrbitIteratorProxy<Base> TupleCosets()
		{
			return OrbitIteratorProxy<Base>();
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
		inline static bool IsCanonical(const TupleT<Base>& tuple)
		{
			return IsCanonicalImpl<Generators...>(tuple, { tuple });
		}
	};


	void test() {
		using sgs = StrongGeneratingSet<
			SGSGenerator<2, Group<
			Permutation<0, 4, 5, 6, 1, 2, 3>
			>>,
			SGSGenerator<3, Group<
			Permutation<0, 1, 3, 2, 4, 5, 6>
			>>,
			SGSGenerator<6, Group<
			Permutation<0, 1, 2, 3, 4, 6, 5>
			>>
			> ;

		using TupleT = sgs::TupleT<100>;
		TupleT vec{ 0, 0, 0, 1, 0, 0, 0 };
		print(vec);

		std::cout << std::boolalpha << sgs::IsCanonical<100>(vec) << "\n";

		size_t count = 0;
		for (const auto& tuple : sgs::TupleCosets<8>()) {
			if ((count++ % 10000) > 0)
				continue;
			print(tuple);
		}
	}
}