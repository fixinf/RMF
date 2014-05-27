

#include "set_const.h"
#include "gsl/gsl_vector.h"
#include "EoS.h"
#include "constants.h"
#include "Auxillary.h"
#include "gsl/gsl_errno.h"
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_deriv.h>
#include <iostream>

struct func_params{
	set_const * C;
	double f0;
	double E0;
	double K0;
	double A0;
	double n0;
};

double map_l = 0.0;
double map_u = 0.1;
double map_range = 10.0;
int runs = 0;
void map(double & what, double x_l, double x_u, double range){
	what = x_l + (1.0 / pi)*(x_u - x_l)*(atan(what / range) + pi / 2);
}

int func(const gsl_vector * x, void * params, gsl_vector * f){
	func_params * p = (func_params *)params;
	set_const * C = p->C;
	double f0 = p->f0;
	double E0 = p->E0;
	double K0 = p->K0;
	double A0 = p->A0;
	double n0 = p->n0;
// 
// 	printf(" Cs = %.15f, Co = %.15f, Cr = %.15f, b = %.15f, c = %.15f \n",
// 		gsl_vector_get(x, 0),
// 		gsl_vector_get(x, 1), 
// 		gsl_vector_get(x, 2), 
// 		gsl_vector_get(x, 3), 
// 		gsl_vector_get(x, 4));
	C->C_s = fabs(gsl_vector_get(x, 0));
	C->C_o = fabs(gsl_vector_get(x, 1));
	C->C_r = fabs(gsl_vector_get(x, 2));
	C->b = gsl_vector_get(x, 3);
	//C->c = abs(gsl_vector_get(x, 4));
	C->c = gsl_vector_get(x, 4);

	//map(C->C_s, 0.0, 20.0, map_range);
	//map(C->c, map_l, map_u, map_range);
	//map(C->b, map_l, map_u, map_range);
	if (C->c < 0){
		int i = 0;
	}
	//printf(" Cs = %f, Co = %f, Cr = %f, b = %f, c = %f \n",
	//	C->C_s, C->C_o, C->C_r, C->b, C->c);
	double temp = 0.0;
	temp = Auxillary::SNM_BED(n0, C);// - (E0);
	//printf("Cs= %f, temp1=%f\n", C->C_s, temp);
	temp = temp - E0;
	gsl_vector_set(f, 0, temp);
	double dn = 1e-6;
	double Ep = Auxillary::SNM_BED(n0 + dn, C);
	double Em = Auxillary::SNM_BED(n0 - dn, C);

	temp = (Ep - Em) / (2.0*dn);
	gsl_vector_set(f, 1, temp);
	temp = Auxillary::Compr(C) - K0;
	gsl_vector_set(f, 2, temp);
	temp = Auxillary::Asymm(n0, C) - A0;
	gsl_vector_set(f, 3, temp);
	temp = EoS::f_eq(n0/2, n0/2, C) - f0;
	gsl_vector_set(f, 4, temp);

	return GSL_SUCCESS;
}

struct f_params{
	set_const * C;
	int index_x;
	int index_y;
	const gsl_vector * x;
};

gsl_vector * z = gsl_vector_alloc(5);
gsl_vector * f = gsl_vector_alloc(5);

double fun(double x, void * params){
	f_params * p = (f_params *)params;
	set_const * C = p->C;

	for (int i = 0; i < 5; i++){
		gsl_vector_set(z, i, gsl_vector_get(p->x, i));
	}

	gsl_vector_set(z, p->index_x, x);
	func(z, C, f);
	double y = gsl_vector_get(f, p->index_y);
	return y;
}

int Jacobian(const gsl_vector * x, void * params, gsl_matrix * J){
	set_const * C = (set_const *)params;
	double dz[5] = { 1e-5, 1e-5, 1e-5, 1e-5, 1e-5};
	gsl_function ff;
	ff.function = &fun;
	f_params fpar;
	fpar.C = C;
	fpar.x = x;
	double res, err;
	for (int i = 0; i < 5; i++){
		for (int j = 0; j < 5; j++){
			fpar.index_y = i;
			fpar.index_x = j;
			ff.params = &fpar;
			gsl_deriv_central(&ff, gsl_vector_get(x, j), dz[j], &res, &err);
		//	printf("%f , %f \n", gsl_vector_get(x, j), res);
			gsl_matrix_set(J, i, j, res);
		}
	}
	return GSL_SUCCESS;
}


int func_fdf(const gsl_vector * x, void *params,
	gsl_vector * f, gsl_matrix * J){

	func(x, params, f);
	Jacobian(x, params, J);

	return GSL_SUCCESS;
}

// int
// print_state(size_t iter, gsl_multiroot_fsolver * s)
// {
// 	printf("iter = %3u x = % .7f % .7f % .7f % .7f % .7f \n "
// 		"f(x) = % .7f % .7f % .7f % .7f % .7f \n",
// 		iter,
// 		gsl_vector_get(s->x, 0),
// 		gsl_vector_get(s->x, 1),
// 		gsl_vector_get(s->x, 2),
// 		gsl_vector_get(s->x, 3),
// 		gsl_vector_get(s->x, 4),
// 		gsl_vector_get(s->f, 0),
// 		gsl_vector_get(s->f, 1),
// 		gsl_vector_get(s->f, 2),
// 		gsl_vector_get(s->f, 3),
// 		gsl_vector_get(s->f, 4)
// 		);
// 	return 0;
// }
// 
// 
// using namespace std;
// double solve(set_const * C){
// 	const gsl_multiroot_fsolver_type *T;
// 	gsl_multiroot_fsolver *s;
// 
// 	int status;
// 	size_t i, iter = 0;
// 
// 	const size_t n = 5;
// 	//struct rparams p = {1.0, 10.0};
// 	gsl_multiroot_function m_f = {
// 		&func,
// 		n, C
// 	};
// 
// 	//double x_init[5] = { C->C_s, C->C_o, C->C_r, C->b, C->c };
// 	double x_init[5] = { 13.0, 9.0, 10.0, 0.004, 0.002 };
// 	gsl_vector * x = gsl_vector_alloc(n);
// 
// 	for (int i = 0; i < n; i++){
// 		gsl_vector_set(x, i, x_init[i]);
// 	}
// 
// 
// 	f_params fpar = { C, 0, 0, x };
// 	cout << "F = " << f(x_init[0], &fpar) << endl;
// 
// 	gsl_matrix * Jac = gsl_matrix_alloc(5, 5);
// 	Jacobian(x, C, Jac);
// 	cout << "Jacobian = " << gsl_matrix_get(Jac, 0, 0) << endl;
// 
// 	for (int i = 0; i < 5; i++){
// 		for (int j = 0; j < 5; j++){
// 			printf("%.15f  ", gsl_matrix_get(Jac, i, j));
// 		}
// 		printf("\n");
// 	}
// 
// 	T = gsl_multiroot_fsolver_hybrid;
// 	s = gsl_multiroot_fsolver_alloc(T, n);
// 	gsl_multiroot_fsolver_set(s, &m_f, x);
// 
// 	print_state (iter, s);
// 
// 	do
// 	{
// 		iter++;
// 		printf("DX = %.27f ", s->dx);
// 		status = gsl_multiroot_fsolver_iterate(s);
// 
// 		print_state(iter, s);
// 
// 		if (status)   /* check if solver is stuck */
// 			break;
// 
// 		status =
// 			gsl_multiroot_test_residual(s->f, 1e-7);
// 	} while (status == GSL_CONTINUE && iter < 1000);
// 
// 	printf("status = %s\n", gsl_strerror(status));
// 
// 	gsl_multiroot_fsolver_free(s);
// 	gsl_vector_free(x);
// 	return 0;
// }




int
	print_state (size_t iter, gsl_multiroot_fdfsolver * s)
{
		double cc = gsl_vector_get(s->x, 4);
		//map(cc, map_l, map_u, map_range);
	printf ("iter = %3u x = % .7f % .7f % .7f % .7f % .7f \n "
		"f(x) = % .7f % .7f % .7f % .7f % .7f \n",
		iter,
		gsl_vector_get (s->x, 0), 
		gsl_vector_get (s->x, 1),
		gsl_vector_get (s->x, 2),
		gsl_vector_get (s->x, 3),
		cc,
		gsl_vector_get (s->f, 0), 
		gsl_vector_get (s->f, 1),
		gsl_vector_get (s->f, 2),
		gsl_vector_get (s->f, 3),
		gsl_vector_get (s->f, 4)
		);
	return 0;
}


using namespace std;
double solve(set_const * C, double n0, double f0, 
	double E0, double K0, double A0){
	const gsl_multiroot_fdfsolver_type *T;
	gsl_multiroot_fdfsolver *s;

	int status;
	size_t i, iter = 0;

	const size_t n = 5;
	//struct rparams p = {1.0, 10.0};
	func_params params = { C, f0, E0, K0, A0, n0 };
	gsl_multiroot_function_fdf m_f = {
		&func,
		&Jacobian,
		&func_fdf,
		n, &params
	};
		
	double x_init[5] = { C->C_s, C->C_o, C->C_r, C->b, C->c };
	//double x_init[5] = { 13.0, 9.0, 10.0, 0.004, 0.001 };
	gsl_vector * x = gsl_vector_alloc (n);

	for (int i = 0; i<n; i++){
		gsl_vector_set(x, i, x_init[i]);
	}


// 	f_params fpar = { C, 0, 0, x };
// 	cout << "F = " << f(x_init[0], &fpar) << endl;
// 
// 	gsl_matrix * Jac = gsl_matrix_alloc(5, 5);
// 	Jacobian(x, C, Jac);
// 	cout << "Jacobian = " << gsl_matrix_get(Jac, 0, 0) << endl;
// 	
// 	for (int i = 0; i < 5; i++){
// 		for (int j = 0; j < 5; j++){
// 			printf("%.15f  ", gsl_matrix_get(Jac, i, j));
// 		}
// 		printf("\n");
// 	}

	T = gsl_multiroot_fdfsolver_hybridsj;
	s = gsl_multiroot_fdfsolver_alloc(T, n);
	gsl_multiroot_fdfsolver_set (s, &m_f, x);

//	print_state (iter, s);

	do
	{
		iter++;
		//printf("DX = %.20f", s->dx);
		status = gsl_multiroot_fdfsolver_iterate (s);

//		print_state (iter, s);

		if (status)   /* check if solver is stuck */
			break;

		status = 
			gsl_multiroot_test_residual (s->f, 1e-6);
	}
	while (status == GSL_CONTINUE && iter < 1000);
	if (status != 0){
		printf ("status = %s\n", gsl_strerror (status));
	}

	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free (x);
	return 0;	
}
