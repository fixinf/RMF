/*
 * fitting.h
 *
 *  Created on: 20.05.2013
 *      Author: fixinf
 */



#include <gsl/gsl_vector.h>
#include "set_const.h"
#ifndef FITTING_H_
#define FITTING_H_
struct fit_eq_params{
	double E0, n0, K0, Asym, f0, z;
};


int fit(fit_eq_params, set_const *);
int func_fit(const gsl_vector *, void * , gsl_vector *);




#endif /* FITTING_H_ */
