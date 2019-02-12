#pragma once

#include <algorithm>

namespace CoffeeCode::NautyLink {
	#if MAXN != (K_SYS+K_ENV)
	#error "MAXN should be equal to K_SYS+K_ENV"
	#endif

	// maximum number of vertices equals k_tot
	constexpr size_t K_TOT = (K_SYS + K_ENV);

	// conflict with libpopcount
	#undef HAVE_POPCNT

	// nauty header
	// note that it should be customized to the target hardware
	// on windows, open nauty.h and set HAVE_UNISTD_H to 0
	#include <nauty.h>

	namespace {
		// implementation details due to nauty interface

		// group order
		// this is calling legacy C code so we allow a global variable
		// TODO: exchange for a multiprecision integer since this number can be huge
		size_t __grouporder;
		void UserLevelProc(int*, int*, int, int*, statsblk*, int, int index, int, int, int, int)
		{
			__grouporder *= index;
		}
	}


	// make sure we picked right word size for nauty sets
	#if MAXN <= 64 && WORDSIZE < MAXN
	#error "consider increasing WORDSIZE to " MAXN " for a performance boost"
	#endif

	// number of setwords needed to store K_TOT bits
	constexpr static auto K_TOT_SETWORDS = SETWORDSNEEDED(K_TOT);

	// nauty graph instance
	template<typename MatrixT>
	class NautyLink {
		graph G[K_TOT*K_TOT_SETWORDS];
		int lab[K_TOT];
		int ptn[K_TOT];
		int orbits[K_TOT];
		statsblk stats;

		enum Color {
			// 4 colors at most since our Us have at most 4 different states per vertex
			C0 = 0,
			C1 = 1,
			C2 = 2,
			C3 = 3,

			ENVIRONMENT_COLOR = 4,

			COLOR_COUNT = 5
		};
		using ColoringT = std::array<Color, K_TOT>;
		template<typename T>
		using PartialColoringT = std::array<T, K_SYS>;

	public:
		NautyLink(const MatrixT& M) :
			lab{ 0 },
			ptn{ 0 },
			orbits{ 0 }
		{
			EMPTYGRAPH(G, K_TOT_SETWORDS, K_TOT);

			// verify matrix given has right shape
			static_assert(MatrixT::row_count == MatrixT::column_count && MatrixT::row_count == K_TOT);
			static_assert(MatrixT::k_sys == K_SYS && MatrixT::k_env == K_ENV);

			// verify linked nauty has correct parameters
			nauty_check(WORDSIZE, static_cast<int>(K_TOT_SETWORDS), static_cast<int>(K_TOT), NAUTYVERSIONID);

			// copy CoffeeCode matrix to nauty matrix
			// since we anticipate only doing this once, we go for portability instead of memcpy-ing the bits
			using BitT = typename MatrixT::RowVectorT::BitT;
			for (size_t i = 0; i < MatrixT::row_count; i++)
				for (size_t j = 0; j < MatrixT::column_count; j++)
					if (M.rows[i][j] == BitT{ 1 })
						ADDONEEDGE(G, i, j, K_TOT_SETWORDS);

			// set default coloring such that environment qubits have ID 2
			ColoringT coloring{ 0 };
			for (size_t i = K_SYS; i < K_TOT; i++)
				coloring[i] = ENVIRONMENT_COLOR;
			SetColoring(coloring);
		}

		auto GroupOrder()
		{
			// default options
			static DEFAULTOPTIONS_GRAPH(options);
			__grouporder = 1;
			options.userlevelproc = UserLevelProc;
			options.defaultptn = false;
			options.writeautoms = false;

			densenauty(G, lab, ptn, orbits, &options, &stats, K_TOT_SETWORDS, K_TOT, NULL);

			return __grouporder;
		}

		template<typename T>
		void SetColoring(const PartialColoringT<T>& partial)
		{
			// set default coloring such that environment qubits have ID 2
			ColoringT coloring;
			for (size_t i = 0; i < K_SYS; i++)
				coloring[i] = static_cast<Color>(partial[i]);
			for (size_t i = K_SYS; i < K_TOT; i++)
				coloring[i] = ENVIRONMENT_COLOR;
			SetColoring(coloring);
		}

		void SetColoring(const ColoringT& coloring)
		{
			using LabEntryT = int;

			std::array<size_t, COLOR_COUNT> bucket_idx{ 0 };
			std::array<std::array<LabEntryT, K_TOT>, COLOR_COUNT> lab_bucket;

			// first we sort the vertices in subsets corresponding to their color
			for (size_t i = 0; i < coloring.size(); i++) {
				const auto color = coloring[i];

				// put vertex in right bucket and increment index for that color
				lab_bucket[color][bucket_idx[color]] = static_cast<LabEntryT>(i);
				bucket_idx[color]++;
			}

			// then we merge the buckets into the lab vector
			size_t lab_idx = 0;
			for (size_t i = 0; i < lab_bucket.size(); i++) {
				std::copy_n(lab_bucket[i].begin(), bucket_idx[i], lab + lab_idx);
				lab_idx += bucket_idx[i];

				// accumulate bucket indices
				bucket_idx[i] = lab_idx;
			}


			// mark boundaries on ptn vector
			for (size_t i = 0; i < coloring.size(); i++)
				ptn[i] = 1;
			for (size_t i = 0; i < COLOR_COUNT; i++)
				ptn[bucket_idx[i]-1] = 0;
		}
	};
}