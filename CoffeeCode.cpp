#ifndef K_SYS
#error "need to specify -DK_SYS=number of system qubits"
#endif

#ifndef K_ENV
#error "need to specify -DK_ENV=number of environment qubits"
#endif


#ifdef SYMMETRIC_SOLVER
int SymmetricSolver();
#else
int FullSolver();
#endif

int main()
{
#ifdef SYMMETRIC_SOLVER
	return SymmetricSolver();
#else
	return FullSolver();
#endif
}
