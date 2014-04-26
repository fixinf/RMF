/*
 * eq_volkov.cpp
 * Contains Oppenheimer-Volkov's equation in GSL-compatible form. Hard-coded dimensionless constants are evaluated to
 * correspond [E] = [P] = km^-2, [M] = [r] = km, in equation, but EoS must provide them in MeV^4.
 *  Created on: 31.01.2013
 *      Author: fixinf
 */
#include "set_const.h"
#include "EoS.h"
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

using namespace EoS;


int eq_volkov(double r, const double y[], double f[], void *params) {
	bool debug = false;
	double A = 1.4766;
	double B = 1.0;
	double C = 1.47934585e-12;
	double D = 2.9532;
	double E1 = 1.479345850e-12;
	if (debug)
		printf("r :: %f P :: %f M :: %.e E :: %f \n", r, y[0], y[1], y[2]);
	f[0] = -A * y[1] * y[2] / (r * r);
	if (debug)
		printf("f[0]:1 :: %f \n", f[0]);
	f[0] *= (1.0 + B * y[0] / y[2]);
	if (debug)
		printf("f[0]:2 :: %f \n", f[0]);
	f[0] *= (1.0 + C * y[0] * pow(r, 3.0) / y[1]);
	if (debug)
		printf("f[0]:3 :: %f \n", f[0]);
	f[0] /= (1.0 - D * y[1] / r);
	if (debug)
		printf("f[0]:4 :: %f \n", f[0]);
	f[1] = r * r * E1 * y[2];
	if (debug)
		printf("f[1] :: %f \n", f[1]);
	f[2] = 0.0;
	return GSL_SUCCESS;
}

struct vofv_params{
	double V0;
	set_const* C;
};


double func_eofp(double x, void * params) {
	vofv_params *p =  (vofv_params * ) params;
	double P0 = p->V0;
	set_const* C = p->C;
	return P(x, C) - P0;
}

const gsl_root_fsolver_type *T_eofp = gsl_root_fsolver_brent;
gsl_root_fsolver *s_eofp = gsl_root_fsolver_alloc(T_eofp);
double EofP(double P, set_const* C) {
	//this stupid root-finder again
	int status;
	int iter = 0, max_iter = 800;
	
	double x = 0.0, x_expect = 10.0;
	double xmin = 0.3e-2, xmax = 20.0;
	gsl_function F;
	F.function = &func_eofp;
	vofv_params params =  {P, C};
	F.params = &params;

	gsl_root_fsolver_set(s_eofp, &F, xmin, xmax);

	do {
		iter++;
		status = gsl_root_fsolver_iterate(s_eofp);
		x = gsl_root_fsolver_root(s_eofp);
		xmin = gsl_root_fsolver_x_lower(s_eofp);
		xmax = gsl_root_fsolver_x_upper(s_eofp);
		status = gsl_root_test_interval(xmin, xmax, 0, 1e-15);
		//if (status == GSL_SUCCESS)
		//printf("Hi, i'm done");
	} while (status == GSL_CONTINUE && iter < max_iter);
	//printf("x::%f \n", x);
	return E(x, C);
}

double func_pofe(double x, void * params) {
	vofv_params *p  = (vofv_params *) params;
	double E0 = p->V0;
	set_const* C = p->C;
	return E(x, C) - E0;
}

const gsl_root_fsolver_type *T_pofe = gsl_root_fsolver_brent;
gsl_root_fsolver *s_pofe = gsl_root_fsolver_alloc(T_pofe);

double PofE(double E, set_const *C) {
	//this stupid root-finder again
	int status;
	int iter = 0, max_iter = 800;
	double x = 1e-3, x_expect = 10.0;
	double xmin = 0.0, xmax = 20.0;
	gsl_function F;
	F.function = &func_pofe;
	vofv_params params = {E, C};
	F.params = &params;

	gsl_root_fsolver_set(s_pofe, &F, xmin, xmax);

	do {
		iter++;
		status = gsl_root_fsolver_iterate(s_pofe);
		x = gsl_root_fsolver_root(s_pofe);
		xmin = gsl_root_fsolver_x_lower(s_pofe);
		xmax = gsl_root_fsolver_x_upper(s_pofe);
		status = gsl_root_test_interval(xmin, xmax, 0, 1e-15);
		//if (status == GSL_SUCCESS)
		//printf("Hi, i'm done");
	} while (status == GSL_CONTINUE && iter < max_iter);

	return P(x, C);
}
