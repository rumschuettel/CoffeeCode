// this is a compile-time problem instance as exported from Mathematica
struct graphstate_instance {
	using sgs = SGSTransversal<
		SGSGenerator<1, Group<
		Permutation<2, 1, 0>
		>>
		>;
	constexpr static AdjacencyMatrixT<4> adjacency_matrix{ {{0, 1, 0, 0}, {1, 0, 1, 1}, {0, 1, 0, 0}, {0, 1, 0, 0}} };
};