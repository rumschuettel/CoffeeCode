int FullSolver();
int SymmetricSolver();

#define SYMMETRIC_SOLVER

int main()
{
#ifdef SYMMETRIC_SOLVER
	return SymmetricSolver();
#else
	return FullSolver();
#endif
}
