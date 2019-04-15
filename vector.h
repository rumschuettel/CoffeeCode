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
		constexpr Vector() noexcept : vec{ 0 } {};
		constexpr Vector(const StoreT vec) noexcept : vec{ vec } {};
		constexpr Vector(const Vector& c) noexcept = default;

	private:
		template<size_t... Idx>
		constexpr Vector(const std::array<BitT, Width>& entries, std::index_sequence<Idx...>) noexcept
			: vec{static_cast<StoreT>( 
					((static_cast<StoreT>(entries[Idx]) << Idx) ^ ... ^ 0)
				)}
		{
		}
	public:
		constexpr explicit Vector(const std::array<BitT, Width>& entries) noexcept
			: Vector(entries, std::make_index_sequence<Width>{})
		{
		}

		// implicit conversion to StoreT
		inline explicit operator StoreT() const noexcept
		{
			return vec;
		}

		// get number of 1s
		inline size_t popcount() const noexcept
		{
			return Popcount(vec);
		}

		// addition
		// pass by value
		inline Vector operator+(const Vector rhs) const noexcept
		{
			return vec ^ rhs.vec;
		}
		inline Vector& operator+=(const Vector rhs) noexcept
		{
			vec ^= rhs.vec;
			return *this;
		}

		// dot product
		// pass by value
		inline StoreT operator*(const Vector rhs) const noexcept
		{
			return Popcount(vec & rhs.vec) & 0b01;
		}

		// bit access
		constexpr inline BitT operator[](size_t index) const noexcept
		{
			return !!(vec & (StoreT{ 0b01 } << index));
		}

		// to string
		friend std::ostream& operator<< (std::ostream& stream, const Vector& vector)
		{
			for (size_t i = 0; i < Width - 1; i++)
				stream << +vector[i] << " ";
			stream << +vector[Width - 1];
			return stream;
		}
	};
}