/*
 * EoS.h
 *
 *  Created on: 16.05.2013
 *      Author: fixinf
 */

#ifndef EOS_H_
#define EOS_H_
#include "set_const.h"
#include <boost/python.hpp>
#include <vector>
#include <gsl/gsl_vector_double.h>

typedef boost::python::list pylist;

namespace EoS{

struct func_f_eq_multi_params{
	std::vector<double> n;
	set_const * C;
};
double E(double, set_const*);
double t_E(double, double, double, set_const*);
double t_E(std::vector<double>, double, set_const*);
double t_E(double, double, set_const*);

double p_f(double);
double f_eq(double, double, set_const*);
double f_eq(vec, set_const*);
double np_eq(double, vec, set_const*);
//int np_eq(double ntot, vec & result, set_const * C);
double m_eff(double);
double P(double, set_const*);
int func_np(gsl_vector *x, void * params, gsl_vector *f);
double mu_e(double, double, double, set_const *);
double mu(int i, vec n, set_const * C);

struct func_np_params {
	double ntot;
	set_const * C;
	vec n_ext;
};

struct func_f_eq_params{
	double nn, np;
	set_const * C;
};
double func_f_eq(double f, void * params);
//double mu_e(double n, double np, double f, set_const * C);
double func_f_eq_multi(double f, void * params);



}

#endif /* EOS_H_ */
