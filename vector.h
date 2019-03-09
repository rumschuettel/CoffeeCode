#pragma once

#include "types.h"

#include <cctype>
#include <iostream>
#include <algorithm>
#include <tuple>

namespace CoffeeCode {
	// binary vector class
	template <size_t Width>
	struct Vector {
		using StoreT = BitStorageType<Width>;
		using BitT = BitType;

		// enough bits to store vector?
		static_assert(Width <= 8 * sizeof(StoreT), "StoreT too short");

		StoreT vec;

		// default constructors
		constexpr Vector() : vec{ 0 } {};
		constexpr Vector(const StoreT vec) : vec{ vec } {};
		constexpr Vector(const Vector& c) = default;

	private:
		template<size_t... Idx>
		constexpr Vector(const std::array<BitT, Width>& entries, std::index_sequence<Idx...>)
			: vec{static_cast<StoreT>( 
					(Bitmask<StoreT, Idx>::mask1000 ^ ... ^ 0)
				)}
		{
		}
	public:
		constexpr explicit Vector(const std::array<BitT, Width>& entries)
			: Vector(entries, std::make_index_sequence<Width>{})
		{
		}

		// implicit conversion to StoreT
		inline explicit operator StoreT() const {
			return vec;
		}

		// get number of 1s
		inline size_t popcount() const {
			return Popcount(vec);
		}

		// addition
		// pass by value
		inline Vector operator+(const Vector rhs) const {
			return vec ^ rhs.vec;
		}
		inline Vector& operator+=(const Vector rhs) {
			vec ^= rhs.vec;
			return *this;
		}

		// dot product
		// pass by value
		inline StoreT operator*(const Vector rhs) const {
			return Popcount(vec & rhs.vec) & 0b01;
		}

		// bit access
		constexpr inline BitT operator[](size_t index) const {
			return !!(vec & (StoreT{ 0b01 } << index));
		}

		// to string
		friend std::ostream& operator<< (std::ostream& stream, const Vector& vector) {
			for (size_t i = 0; i < Width - 1; i++)
				// cast to normal int, otherwise it'll print chars
				stream << static_cast<uint32_t>(vector[i]) << " ";
			stream << static_cast<uint32_t>(vector[Width - 1]);
			return stream;
		}
	};
}