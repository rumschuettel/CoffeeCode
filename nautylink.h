#pragma once

#include "nauty-color-generator-vcolg.h"

#include "ctmath.h"
#include "types.h"

#include "utility.h"

#include <algorithm>
#include <boost/container_hash/hash.hpp>

#include <tuple>

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



	// implementation details for nauty interface
	namespace ImplDetails {
		// group order
		// this is calling legacy C code so we allow a global variable
		using OrbitSizeT = CoffeeCode::OrbitType;
		static thread_local OrbitSizeT __grouporder;


		void UserLevelProc_GroupOrder(int*, int*, int, int*, statsblk*, int, int index, int, int, int, int)
		{
			assert(index > 0);
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
		decltype(G) G_canon;
		int lab[K_TOT];
		int ptn[K_TOT];
		int orbits[K_TOT];
		statsblk stats;

		// not static since we need to first call nauty to fill it with a value
		ImplDetails::OrbitSizeT FULL_GROUP_ORDER;

		enum Color {
			// 4 colors at most since our Us have at most 4 different states per vertex
			C0 = 0b00,
			C1 = 0b01,
			C2 = 0b10,
			C3 = 0b11,

			ENVIRONMENT_COLOR = 0b100,

			// we xor the environment color on top
			E0 = C0 | ENVIRONMENT_COLOR,
			E1 = C1 | ENVIRONMENT_COLOR,
			E2 = C2 | ENVIRONMENT_COLOR,
			E3 = C3 | ENVIRONMENT_COLOR,

			COLOR_COUNT = 8
		};

		// map integer values from CoffeeCode to graph colors in Nauty
		template<typename NumberT>
		static inline Color ToColor(NumberT nc) noexcept
		{
			// 4 colors
			assert(0 <= nc && nc <= 3);

			// cast
			return static_cast<Color>(nc);
		}

		using ColoringT = std::array<Color, K_TOT>;
		// colorings where environment is not yet specified
		using PartialColoringT = NibbleStorageTypeArray<4, K_SYS>;
		using ColoringRawT = typename MatrixT::RowVectorT::StoreT;

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
			nauty_check(WORDSIZE, checked_cast<int>(K_TOT_SETWORDS), checked_cast<int>(K_TOT), NAUTYVERSIONID);

			// copy CoffeeCode matrix to nauty matrix
			// since we anticipate only doing this once, we go for portability instead of memcpy-ing the bits
			using BitT = typename MatrixT::RowVectorT::BitT;
			for (size_t i = 0; i < MatrixT::row_count; i++)
				for (size_t j = 0; j < MatrixT::column_count; j++)
					if (M.rows[i][j] == BitT{ 1 })
						ADDONEEDGE(G, i, j, K_TOT_SETWORDS);

			// set default coloring such that environment qubits have ID 2
			ColoringT coloring{ C0 };
			for (size_t i = K_SYS; i < K_TOT; i++)
				coloring[i] = ENVIRONMENT_COLOR;
			SetColoring(coloring);

			// calculate full group order with default coloring
			FULL_GROUP_ORDER = GroupOrder();
		}


		template<size_t Colors>
		void Colorings(const size_t CallbackStride, const size_t CallbackOffset, const CosetGeneratorCallbackType<Colors>& callback) noexcept
		{
			vcolg([&](const int* cols, size_t counter) -> void {
				// get multiplicity of coloring
				PartialColoringT system_coloring;
				for (size_t i = 0; i < K_SYS; i++)
					system_coloring[i] = checked_cast<Color>(cols[i]);

				SetColoring(system_coloring);
				callback(system_coloring, ColoringMultiplicity<Colors>(), counter);
			}, G, Colors, CallbackStride, CallbackOffset);
		}


		auto inline GroupOrder() noexcept
		{
			// default options
			static DEFAULTOPTIONS_GRAPH(options);
			ImplDetails::__grouporder = 1;
			options.userlevelproc = ImplDetails::UserLevelProc_GroupOrder;
			options.defaultptn = false;

			densenauty(G, lab, ptn, orbits, &options, &stats, K_TOT_SETWORDS, K_TOT, nullptr);

			return ImplDetails::__grouporder;
		}

		template<size_t ColorCount>
		auto inline ColoringMultiplicity() noexcept
		{
			return static_cast<MultiplicityType<ColorCount>>(FULL_GROUP_ORDER / GroupOrder());
		}

		struct CanonicalImage {
			CanonicalImage() = default;
			CanonicalImage(const CanonicalImage&) = default;
			CanonicalImage(const decltype(G_canon)& in, const ColoringT& coloring) noexcept
				: vertexColorCounts{ 0 }
			{
				std::copy(std::begin(in), std::end(in), std::begin(G_canon));
				for (const auto c : coloring)
					vertexColorCounts[c]++;
			}
			// stores canonical image
			decltype(G) G_canon;
			// stores number of colors in coloring
			NibbleStorageTypeArray<K_TOT, COLOR_COUNT> vertexColorCounts;

			// hash for this type
			struct Hash {
				inline std::size_t operator()(CanonicalImage const& image) const noexcept
				{
					std::size_t h = boost::hash_range(std::begin(image.G_canon), std::end(image.G_canon));
					std::size_t h2 = boost::hash_range(std::begin(image.vertexColorCounts), std::end(image.vertexColorCounts));
					boost::hash_combine(h, h2);
					return h;
				}
			};
			// comparison for this type
			bool operator==(const CanonicalImage& rhs) const noexcept
			{
				return (vertexColorCounts == rhs.vertexColorCounts) && std::equal(std::begin(G_canon), std::end(G_canon), std::begin(rhs.G_canon));
			}
		};
		

		using CanonicalImageT = CanonicalImage;
		using CanonicalImageHashT = typename CanonicalImage::Hash;


		// reorders to get canonical image of the partial coloring given
		// also returns the group order
		using CanonicalColoringReturnT = std::pair<CanonicalImage, MultiplicityType<2>>;
		CanonicalColoringReturnT CanonicalColoring(const ColoringRawT raw_coloring) noexcept
		{
			// memoization
			thread_local unordered_map<ColoringRawT, CanonicalColoringReturnT> memo;
			typename decltype(memo)::const_iterator p = memo.find(raw_coloring);
			if (p != memo.end()) return p->second;

			// get canonical coloring
			const auto coloring = SetColoring(raw_coloring);

			static DEFAULTOPTIONS_GRAPH(options);
			ImplDetails::__grouporder = 1;
			options.userlevelproc = ImplDetails::UserLevelProc_GroupOrder;
			options.defaultptn = false;
			options.getcanon = true;

			densenauty(G, lab, ptn, orbits, &options, &stats, K_TOT_SETWORDS, K_TOT, G_canon);

			// build canonical image type, memoize
			CanonicalImage out(G_canon, coloring);

			const CanonicalColoringReturnT ret = std::make_pair(
				out,
				static_cast<MultiplicityType<2>>(FULL_GROUP_ORDER / ImplDetails::__grouporder)
			);

			// cache; if memory full, flush
			try {
				memo[raw_coloring] = ret;
			} catch (std::bad_alloc&) {
				memo.swap(decltype(memo){});
				memo[raw_coloring] = ret;
			}
			return ret;
		}

		auto inline SetColoring(const ColoringRawT partial) noexcept
		{
			ColoringT coloring;
			// careful: coloring[0] is the first vertex, which corresponds to the rightmost bit of the partial
			for (size_t i = 0; i < K_TOT; i++) 
				coloring[i] = ToColor((partial & (static_cast<ColoringRawT>(0b01) << i)) >> i);

			// mark environment colors
			for (size_t i = K_SYS; i < K_TOT; i++)
				coloring[i] = static_cast<Color>(coloring[i] | ENVIRONMENT_COLOR);
			return SetColoring(coloring);
		}

		// set coloring given as an array for all system vertices, by converting to internal color type
		auto inline SetColoring(const PartialColoringT partial) noexcept
		{
            // set default coloring such that environment qubits have an extra flagged color
			ColoringT coloring;
			for (size_t i = 0; i < K_SYS; i++)
				coloring[i] = ToColor(partial[i]);
			for (size_t i = K_SYS; i < K_TOT; i++)
				coloring[i] = ENVIRONMENT_COLOR;
			return SetColoring(coloring);
		}

		// set coloring given as a full array
		auto inline SetColoring(const ColoringT& coloring) noexcept
		{
			using LabEntryT = int;

			std::array<size_t, COLOR_COUNT> bucket_idx{ 0 };
			std::array<std::array<LabEntryT, K_TOT>, COLOR_COUNT> lab_bucket;

			// first we sort the vertices in subsets corresponding to their color
			for (size_t i = 0; i < coloring.size(); i++) {
				const auto color = coloring[i];

				// put vertex in right bucket and increment index for that color
				lab_bucket[color][bucket_idx[color]] = checked_cast<LabEntryT>(i);
				bucket_idx[color]++;
			}

			// then we merge the buckets into the lab vector
			size_t lab_idx = 0;
			for (size_t i = 0; i < lab_bucket.size(); i++) {
				std::copy_n(lab_bucket[i].begin(), bucket_idx[i], lab + lab_idx);
				lab_idx += bucket_idx[i];
			}

			// mark boundaries on ptn vector
			for (size_t i = 0; i < coloring.size(); i++)
				ptn[i] = 1;
			size_t bucket_idx_acc{ 0 };
			for (size_t i = 0; i < COLOR_COUNT; i++) {
				if (bucket_idx[i] == 0) continue;

				bucket_idx_acc += bucket_idx[i];
				ptn[bucket_idx_acc - 1] = 0;
			}

			// return for future use
			return coloring;
		}
	};
}
