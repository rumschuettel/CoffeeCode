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
