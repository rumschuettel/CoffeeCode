#pragma once

#include "vector.h"
#include "utility.h"

#include <assert.h>
#include <utility>
#include <tuple>

namespace CoffeeCode {
	// binary matrix class
	template<size_t _row_count, size_t _column_count>
	struct Matrix {
		constexpr static size_t row_count = _row_count;
		constexpr static size_t column_count = _column_count;

		using RowVectorT = Vector<column_count>;
		using ColumnVectorT = Vector<row_count>;
		using BitT = typename RowVectorT::BitT;

		std::array<RowVectorT, row_count> rows;


		// uninitialized matrix
		Matrix() = delete;

		// for compile-time construction we want explicit numbers
	private:
		template<size_t... Idx>
		constexpr Matrix(const std::array<std::array<BitT, column_count>, row_count>& entries, const std::index_sequence<Idx...>) noexcept
			: rows{ RowVectorT(entries[Idx])... }
		{
		}
	public:
		constexpr Matrix(const std::array<std::array<BitT, column_count>, row_count>& entries) noexcept
			: Matrix(entries, std::make_index_sequence<row_count>{})
		{
		}
		constexpr Matrix(const std::array<RowVectorT, row_count>& rows) noexcept
			: rows{ rows }
		{
		}

		// dot product
		// pass by value
		inline ColumnVectorT operator*(const RowVectorT rhs) const noexcept
		{
			ColumnVectorT out;
			using ColumnVectorStoreT = typename ColumnVectorT::StoreT;
			size_t i = 0;
			for (const auto row : rows)
				out += checked_cast<ColumnVectorStoreT>( 
					checked_cast<ColumnVectorStoreT>(row * rhs) << i++
				);
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
	template<size_t _k_sys, size_t _k_env>
	struct AdjacencyMatrix : public Matrix<_k_sys + _k_env, _k_sys + _k_env> {
		constexpr static auto k_sys = _k_sys;
		constexpr static auto k_env = _k_env;

		using BaseT = Matrix<k_sys + k_env, k_sys + k_env>;
		using BitT = typename BaseT::BitT;
		using AdjacencyMatrixT = AdjacencyMatrix<k_sys, k_env>;
		using ABBlockT = Matrix<k_sys, k_env>;

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
			std::array<std::array<BitT, column_count>, row_count> entries;
			for (size_t i = 0; i < row_count; i++) {
				for (size_t j = 0; j < column_count; j++)
					entries[i][j] = bitstring[column_count*i + j] == '0' ? 0 : 1;
			}

			// return new matrix
			return AdjacencyMatrixT(entries);
		}

		// return AB block
		constexpr ABBlockT AB() const
		{
			return ABImpl(std::make_index_sequence<k_sys>{});
		}

	private:
		template<size_t... IdxRow>
		constexpr ABBlockT ABImpl(const std::index_sequence<IdxRow...>) const noexcept
		{
			static_assert(sizeof...(IdxRow) == k_sys);

			using ABRowVectorT = typename ABBlockT::RowVectorT;

			// extract subblock from (k_sys, 0) to but not including (k_sys+k_env, k_sys)
			// this works since IdxRow runs from 0...k_sys-1, and IdxCol runs from 0...k_env
			const auto rows = this->rows;
			std::array<ABRowVectorT, k_sys> new_rows{
				ABRowVectorT(static_cast<typename ABRowVectorT::StoreT>(rows[IdxRow].vec >> k_sys)) ...
			};
			return ABBlockT(new_rows);
		}
	};
}