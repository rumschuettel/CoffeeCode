﻿#ifndef K_SYS
#error "need to specify -DK_SYS=number of system qubits"
#endif

#ifndef K_ENV
#error "need to specify -DK_ENV=number of environment qubits"
#endif

#ifdef FLOATING_POINT_MULTIPLICITY
#ifdef REDUCE_LAMBDA_IF_POSSIBLE
#if K_SYS > 16  // means MultiplicityType<4> = CoefficientType for polynomials will be a floating point number
#error "REDUCE_LAMBDA_IF_POSSIBLE does not play well with FLOATING_POINT_MULTIPLICITY, as coefficients need to be compared exactly."
#endif
#endif
#endif

#ifdef PARALLELIZE
#ifndef SYMMETRIC_SOLVER
#error "PARALLELIZE only makes sense for SYMMETRIC_SOLVER"
#endif
#endif

#ifdef SYMMETRIC_SOLVER
#ifdef OPTIMIZE_FOR_DEPOLARIZING
#error "SYMMETRIC_SOLVER and OPTIMIZE_FOR_DEPOLARIZING are incompatible. Polynomial type is specified in cc-instance-custom.h."
#endif
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
