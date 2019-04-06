#pragma once

#include "permutation.h"

#include <set>
#include <vector>

namespace CoffeeCode {

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
			auto begin() const { return IteratorT(tuple); }
			const auto end() const { return typename IteratorT::Done{}; }
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
					( check_and_add(element, Permutations::permute(element)), ... );
			}

			return orbit;
		}

		// apply each generator to given tuple
		template<typename T>
		inline static OrbitSetT<T> GroupAction(const T& tuple)
		{
			OrbitSetT<T> action;
			( action.insert(Permutations::permute(tuple)), ... );
			return action;
		}
	};

}