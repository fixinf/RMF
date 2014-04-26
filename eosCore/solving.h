#include "set_const.h"


double solve(set_const * C, double n0, double f0,
	double E0, double K0, double A0);
int Jacobian(const gsl_vector * x, void * params, gsl_matrix * J);