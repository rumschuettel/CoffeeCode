#pragma once

#include "vector.h"

#include <assert.h>
#include <utility>

namespace CoffeeCode {
	// binary matrix class
	template<size_t row_count, size_t column_count>
	struct Matrix {
		using RowVectorT = Vector<column_count>;
		using ColumnVectorT = Vector<row_count>;

		const std::array<RowVectorT, row_count> rows;

		constexpr static size_t row_count = row_count;
		constexpr static size_t column_count = column_count;

		// uninitialized matrix
		Matrix() = delete;

		// for compile-time construction we want explicit numbers
		template<typename... T>
		constexpr Matrix(const T... rows) : rows{ static_cast<typename RowVectorT::StoreT>(rows)... }
		{}
		using MatrixInitializerRowT = std::array<RowVectorT, row_count>;
		constexpr Matrix(const MatrixInitializerRowT& rows) : rows{ rows }
		{}

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
		using BaseT = Matrix<size_A + size_B, size_A + size_B>;
		using BaseT::Matrix;
		using AdjacencyMatrixT = AdjacencyMatrix<size_A, size_B>;
		using ABBlockT = Matrix<size_A, size_B>;

		using BaseT::row_count;
		using BaseT::column_count;
		using typename BaseT::RowVectorT;

		// initialize from string; clean out whitespace
		static AdjacencyMatrixT FromString(const std::string& bitstring) {
			// trim whitespace
			std::string cleaned{ bitstring };
			const auto predicate = [](unsigned char const c) { return std::isspace(c); };
			cleaned.erase(std::remove_if(cleaned.begin(), cleaned.end(), predicate), cleaned.end());
			assert(cleaned.length() == row_count * column_count);

			// load into rows
			std::array<RowVectorT, row_count> rows;
			for (size_t i = 0; i < row_count; i++) {
				typename RowVectorT::StoreT row{ 0 };
				for (size_t j = 0; j < column_count; j++)
					row |= static_cast<typename RowVectorT::StoreT>((bitstring[column_count*i + j] == '0' ? 0 : 1) << j);
				rows[i] = row;
			}

			// return new matrix
			return AdjacencyMatrixT{ rows };
		}

		// return AB block
		constexpr ABBlockT AB() const
		{
			return ABImpl(std::make_index_sequence<size_A>{});
		}

	private:
		template<size_t... Idx>
		constexpr ABBlockT ABImpl(const std::index_sequence<Idx...>) const
		{
			std::array<typename ABBlockT::RowVectorT, size_A> rows;
			((
				rows[Idx] = static_cast<typename ABBlockT::RowVectorT>((*this).rows[Idx].vec >> size_A)
			), ...);
			return ABBlockT{ rows };
		}
	};
}