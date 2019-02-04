#pragma once

#include "vector.h"

#include <assert.h>

namespace CoffeeCode {
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
				blockAB.rows[i] = static_cast<typename ABBlockT::RowVectorT>((*this).rows[i].vec >> size_A);

			return blockAB;
		}
	};
}