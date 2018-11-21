// Vector is small and used a lot.
// For the operations we expect to use a lot we thus pass them by value, not reference.

#pragma once

#include "libpopcnt.h"
#include <iostream>
#include <algorithm>
#include <cctype>
#include <assert.h>
#include <type_traits>

namespace CoffeeCode {
	// automatic selection of smallest-possible integer
	template<bool raise = true>
	struct InvalidIntegerType {
		static_assert(raise, "not enough bits");
	};

	template<size_t length>
	using StdStoreT = typename std::conditional <
		length <= 8,
		uint8_t,
		typename std::conditional <
		length <= 16,
		uint16_t,
		typename std::conditional <
		length <= 32,
		uint32_t,
		typename std::conditional <
		length <= 64,
		uint64_t,
		InvalidIntegerType<true>
		>::type>::type>::type>::type;

	using StdBitT = StdStoreT<1>;
		
	// fast popcount implementation
	constexpr auto __popcount = popcnt64;


	// binary vector class
	template <size_t length>
	struct Vector {
		using StoreT = StdStoreT<length>;
		using BitT = StdBitT;

		// enough bits to store vector?
		static_assert(length <= 8 * sizeof(StoreT), "StoreT too short");

		StoreT vec;

		// default constructors
		Vector() : vec{0} {};
		Vector(const StoreT vec) : vec{ vec } {};
		Vector(const Vector& c) = default;
		Vector(const std::string& bitstring) : vec{ 0 }
		{
			assert(bitstring.length() == length);
			for (size_t i = 0; i < length; i++) {
				assert(bitstring[i] == '0' || bitstring[i] == '1');
				BitT bit = bitstring[i] == '0' ? 0 : 1;
				vec |= StoreT{ bit } << i;
			}
		}

		// get number of 1s
		inline size_t popcount() const {
			return __popcount(vec);
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
			return __popcount(vec & rhs.vec) & BitT{1};
		}

		// bit access
		inline BitT operator[](size_t index) const {
			return !!(vec & (BitT{1} << index));
		}

		// to string
		friend std::ostream& operator<< (std::ostream& stream, const Vector& vector) {
			for (size_t i = 0; i < length-1; i++)
				// cast to normal int, otherwise it'll print chars
				stream << static_cast<uint32_t>(vector[i]) << " ";
			stream << static_cast<uint32_t>(vector[length - 1]);
			return stream;
		}
	};


	// binary matrix class
	template<size_t row_count, size_t column_count>
	struct Matrix {
		using RowVectorT = Vector<column_count>;
		using ColumnVectorT = Vector<row_count>;

		RowVectorT rows[row_count];

		Matrix() = default;
		Matrix(const std::string bitstring) {
			// trim whitespace
			std::string cleaned{ bitstring };
			const auto predicate = [](unsigned char const c) { return std::isspace(c); };
			cleaned.erase(std::remove_if(cleaned.begin(), cleaned.end(), predicate), cleaned.end());
			assert(cleaned.length() == row_count * column_count);

			// load into rows
			for (size_t i = 0; i < row_count; i++)
				rows[i] = RowVectorT(cleaned.substr(column_count*i, column_count));
		}

		// dot product
		// pass by value
		inline ColumnVectorT operator*(const RowVectorT rhs) const {
			ColumnVectorT out;
			size_t i = 0;
			for (const auto row : rows)
				out += typename ColumnVectorT::StoreT{ row * rhs } << i++;
			return out;
		}

		// to string
		friend std::ostream& operator<< (std::ostream& stream, const Matrix& matrix) {
			for (size_t i = 0; i < row_count - 1; i++)
				stream << matrix.rows[i] << "\n";
			stream << matrix.rows[row_count - 1];
			return stream;
		}
	};


	// adjacency matrix class
	template<size_t size_A, size_t size_B>
	struct AdjacencyMatrix : public Matrix<size_A + size_B, size_A + size_B> {
		using Matrix<size_A + size_B, size_A + size_B>::Matrix;
		using ABBlockT = Matrix<size_A, size_B>;

		// return AB block
		ABBlockT AB() const {
			ABBlockT blockAB;
			for (size_t i = 0; i < size_A; i++)
				blockAB.rows[i] = static_cast<typename ABBlockT::RowVectorT>( (*this).rows[i].vec >> size_A );

			return blockAB;
		}
	};
}