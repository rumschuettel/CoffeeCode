#pragma once

#include "traits.h"
#include "utility.h"

namespace CoffeeCode {
	// compile time tuple indicating start and end of orbit
	template <size_t _Start, size_t _End>
	struct TrivialSGSOrbit
	{
		constexpr static size_t Start = _Start;
		constexpr static size_t End = _End;

        // we only allow non-trivial orbits
		static_assert(Start < End);

		constexpr static size_t Length = End - Start;

        // bitmask for a given type; 1s between Start and End, 0 elsewhere
        template<typename T>
        struct Mask {
            constexpr static T value = Bitmask<T, End>::mask0111 & ~(Bitmask<T, Start>::mask0111);
        };
	};


    // SGSTransversal equivalent featuring a TupleCosets function,
    // as well as a plug-in replacement for NautyLink's CanonicalImage
    // we require that the orbits are in a consecutive block at the beginning of the graph,
    // e.g. {0, 2}, {2, 5}, {5, 7}, {7, 33}
    // where the second index points to the element AFTER the orbit.
    template<size_t _Length, typename... Orbits>
    struct TrivialSGSTransversal
    {
        constexpr static size_t Length = _Length;

		// number of orbits
		constexpr static size_t N_orbits = sizeof...(Orbits);
        static_assert(N_orbits >= 1);

        // number of points not covered by any orbit
        constexpr static size_t N_remaining = Length - sum<Orbits::Length...>;


		// first and last orbit
		using FirstOrbit = nth_element<0, Orbits...>;
        static_assert(FirstOrbit::Start == 0);
        using LastOrbit = nth_element<N_orbits-1, Orbits...>;
        static_assert(LastOrbit::End <= Length);

        // check that orbits are consecutive, e.g. <0, 10>, <10, 12>, ...
        // TODO; for now we just assume that

    private:
        // store orbits' start and end
        constexpr static std::array<size_t, N_orbits> OrbitStarts = { Orbits::Start... };
        constexpr static std::array<size_t, N_orbits> OrbitEnds = { Orbits::End... };
        constexpr static std::array<size_t, N_orbits> OrbitLengths = { Orbits::Length... };


        template<size_t Base>
        struct OrbitProductIterator {
            using IteratorT = OrbitProductIterator<Base>;
			using TupleT = NibbleStorageTypeArray<Base, Length>;
            using MultiplicityT = MultiplicityType<Base>;

        private:
            // within one orbit, we need Base-1 indices to mark the boundaries
            // between the segments 0000 1111 .... Base-1...
            using IndexTupleT = std::array<size_t, Base-1>;
            
            // since we have N_orbits of these segments, we need to store an array of them
            std::array<IndexTupleT, N_orbits> indices;

            // for the remainder of the tuples we have to brute force count to base B
            // that number can be as large as Base^N_remaining
            using BaseCounterT = SizeStorageType<ipow<size_t>(Base, N_remaining)>;
            BaseCounterT remainder;
            constexpr static auto MaxRemainder = ipow<BaseCounterT>(Base, N_remaining);

            enum StepResult { NoCarry, Carry };
			StepResult Step(IndexTupleT& indices, const size_t OrbitLength)
			{
				// if this was the first item and we are saturated,
				// this means all indices are at the rhs, so we are done
                // reset and indicate that we have a carry
				if (indices[0] == OrbitLength) {
                    indices = IndexTupleT{ 0 };
					return Carry;
                }

                // try incrementing any of the indices but the last one
                for (size_t i = 0; i < indices.size()-1; i ++) {

                    // e.g. __A__B_ map to __A_B__
                    if (indices[i] < indices[i+1]) {
                        indices[i]++;
                        return NoCarry;
                    }

                    // A and B in same position: reset B to 0,
                    // then re-enter the loop which will try to increment
                    // the next index A
                    if (indices[i] == indices[i+1]) {
                        indices[i] = 0;
                        continue;
                    }
                }

                // unsuccessful attempt incrementing any but the last index,
                // so all but the last one are now at 0;
                // we only need to verify that the last index is <= OrbitSize
				auto& last = indices[indices.size()-1];
				if (last <= OrbitLength)
					last++;
                /* else does not happen by construction */

				return NoCarry;
			}

        private:
            bool done;

		public:
			OrbitProductIterator() : indices{ 0 }, remainder{ 0 }, done{ false }
			{
			}

			OrbitProductIterator& operator++()
			{
                // we treat all the individual orbit counters as
                // digits and implement a standard mixed base counter here
                for (size_t i = 0; i < N_orbits; i ++) {
                    switch (Step(indices[i], OrbitLengths[i])) {
                        case Carry:
                            // for any but the last orbit, we carry over to the next orbit
                            if (i < N_orbits-1)
                                continue;

                            // last orbit also saturated, so try incrementing remainder
                            if (++remainder == MaxRemainder)
                                done = true;
                            
                            break;

                        case NoCarry:
                            return *this;
                    }
                }

				return *this;
			}

			// comparison operator is only ever invoked for end-of-loop-checking
			// so return the state of our internal "done" variable and ignore the argument
            using DoneT = IteratorT;
			bool operator!=(const DoneT&) const
			{
				return !done;
			}


			// dereference
			const std::pair<TupleT, MultiplicityT> operator*() const
			{
				// this is always the same
				static const OrbitType orbit_enumerator{ (Factorial<OrbitType>(Orbits::Length) * ...) };
				OrbitType orbit_denominator{ 1 };

                TupleT out{0};
                using BaseT = typename TupleT::value_type;

                for (size_t ii = 0; ii < N_orbits; ii ++) {
                    const auto& idcs = indices[ii];
                    //const auto Start = OrbitStarts[ii];
                    const auto End = OrbitEnds[ii];
                    const auto Length = OrbitLengths[ii];

                    // calculate orbit size
                    // this is simply word_length! / prod_i #color_i!
					// todo: make this use MultiplicityT
					orbit_denominator *= Factorial<OrbitType>(idcs[0]);
                    for (size_t i = 0; i < idcs.size()-1; i ++)
						orbit_denominator *= Factorial<OrbitType>(idcs[i+1] - idcs[i]);
					orbit_denominator *= Factorial<OrbitType>(Length - idcs[idcs.size()-1]);
					

                    // fill vector in with numbers 0, 1, 2, ..., Base-1
                    // within sections determined by indices
                    for (size_t i = 0; i < idcs.size()-1; i ++)
                        for (size_t j = idcs[i]; j < idcs[i+1]; j ++)
                            out[End - j - 1] = checked_cast<BaseT>(i+1);
                    for (size_t j = idcs[idcs.size()-1]; j < Length; j ++)
                        out[End - j - 1] = checked_cast<BaseT>(Base-1);
                }

                // fill remainder
                BaseCounterT rem = remainder;
                for (size_t j = 0; j < N_remaining; j ++) {
                    const BaseT digit = rem % Base;
                    rem /= Base;

                    out[j + Length - N_remaining] = digit;
                }

                return { out, checked_cast<MultiplicityT>(orbit_enumerator / orbit_denominator) };
			}
        };

    public:

		// range-based for loop iterator for tuples over Base
		template<size_t Base>
		inline static StatefulIteratorProxy<OrbitProductIterator<Base>> TupleCosets()
		{
			return StatefulIteratorProxy<OrbitProductIterator<Base>>();
		}


        // canonicalization struct
        template<typename MatrixT>
        struct SymmetryProvider {
            using ColoringRawT = typename MatrixT::RowVectorT::StoreT;
            using CanonicalImageT = ColoringRawT;
            using CanonicalImageHashT = std::hash<CanonicalImageT>;
            using MultiplicityT = MultiplicityType<2>;

            SymmetryProvider(const MatrixT&) {};

            inline static auto CanonicalColoring(const ColoringRawT coloring)
            {
                ColoringRawT out = coloring;
				MultiplicityT mult = 1;

                // for every orbit segment, we sort the bit subsets
                // and accumulate the multiplicities
                ( SortBitSubset<Orbits>(out, mult), ... );

                return std::make_pair(out, mult);
            }


        private:
            // sort and return multiplicity, which is simply
            // Orbit::Length choose #1s
            template<typename Orbit>
            inline static void SortBitSubset(ColoringRawT& src, MultiplicityT& mult)
            {
                // bitmask with 1s between orbit start and end
                constexpr ColoringRawT mask = Orbit::template Mask<ColoringRawT>::value;

                // count 1s
                const size_t count1s = Popcount(src & mask);

                // if all were 1s, nothing to do
                if (count1s == Orbit::Length)
                    return;

                // otherwise overwrite section with padded 1s
                const ColoringRawT sorted = ((static_cast<ColoringRawT>(0b01) << count1s) - 1) << Orbit::Start;

                src &= ~mask;
                src ^= sorted;

                // update multiplicity
                mult *= Binomial<MultiplicityT>(Orbit::Length, count1s);
            }

        };
	};
}