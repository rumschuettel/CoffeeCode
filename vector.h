#pragma once

#include "types.h"
#include "libpopcnt.h"

#include <cctype>
#include <iostream>
#include <algorithm>

namespace CoffeeCode {

	// fast popcount implementation
	constexpr auto __popcount = popcnt64;


	// binary vector class
	template <size_t Width>
	struct Vector {
		using StoreT = StdStoreT<Width>;
		using BitT = StdBitT;

		// enough bits to store vector?
		static_assert(Width <= 8 * sizeof(StoreT), "StoreT too short");

		StoreT vec;

		// default constructors
		constexpr Vector() : vec{ 0 } {};
		constexpr Vector(const StoreT vec) : vec{ vec } {};
		constexpr Vector(const Vector& c) = default;

		// get number of 1s
		inline size_t popcount() const {
			return static_cast<size_t>(__popcount(vec));
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
		inline BitT operator*(const Vector rhs) const {
			return __popcount(vec & rhs.vec) & BitT { 1 };
		}

		// bit access
		inline BitT operator[](size_t index) const {
			return !!(vec & (BitT{ 1 } << index));
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