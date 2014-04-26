/*
 * fitting.cpp
 *
 *  Created on: 20.05.2013
 *      Author: fixinf
 */
#include <stdlib.h>
#include <stdio.h>
#include "fitting.h"
#include "set_const.h"
#include "EoS.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "constants.h"
#include <iostream>
#include <math.h>


using namespace std;
using namespace EoS;

double Ener(double n, set_const* C){
	return t_E(0.5*n, 0.5*n, C)/(D*n) - m_n;
}

double Ener(double n, double f, set_const* C){
	return t_E(0.5*n,0.5*n, f, C)/(D*n) - m_n;
}

double K(double n0, double f, set_const* C){
	double meff = m_n*C->phi_n(f);
	double x = p_f(n0);
	double res = 0.5* x *sqrt(meff*meff+x*x)- 0.5* meff*meff*log(x+sqrt(meff*meff+x*x));
	return res * (-m_n*meff)/(pi*pi);
}




//Возвращает уравнение на параметры для nn = np = n/2
int func_fit(const gsl_vector * x, void * params, gsl_vector * f){
	fit_eq_params * p = (fit_eq_params *) params;

	bool debug = 1;

	double E0 = p->E0;
	double n0 = p->n0;
	double K0 = p->K0;
	double f0 = p->f0;
	double Asym0 = p->Asym;
//			  C_s                  C_o                  C_r
	set_const* C = new set_const(gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),
			gsl_vector_get(x,3),gsl_vector_get(x,4),p->z);
//                        b,			c
	double dn = 1e-10;
	
	//TODO Производная от самой энергии

	double dEner_n0 = (Ener(n0+dn, f0, C) - Ener(n0-dn,f0,C))/(2*dn);

	//Энергия в равновесной плотности

	double E_n0 = Ener(n0, f0, C);
	cout << "E_n0_f0 = " << E_n0 <<endl;
	double d2n = 1e-5;
	//TODO вторая производная от самой энергии

//	double d2Ener_n0 = d2E_n0/(D*n0) - 2.0*dE_n0/(D*n0*n0) +
//			2.0*t_E(0.5*n0, 0.5*n0, C)/(D*n0*n0*n0);

	double d2Ener_n0 = (Ener(n0+d2n, f0, C) - 2.0*Ener(n0, f0, C) + Ener(n0-d2n, f0, C))/(d2n*d2n);
	//Выражение для К

	double K = 9.0*n0*n0*d2Ener_n0;


	//double f0 = f_eq(0.5*n0, 0.5*n0, C);

	

	//Каждый раз нужно прописывать эфф. массу
	double m_eff_n0 = C->phi_n(f0);

	// Энергия симметрии

	double Asym_n0 = pow(C->C_r/m_n, 2.0)*D*n0/(8.0*C->eta_r(f0)) +
			pi*pi*D*n0/(4.0*p_f(n0)*sqrt(m_eff_n0*m_eff_n0 + p_f(n0)*p_f(n0)));


	if (debug){
		cout <<"\n ---------------------in func_fit-----------------------" << endl;
		cout << "n0 = " << n0 << endl;
		cout << "E(n0) = " << E_n0 << endl;
		cout << "K0 = " << K << endl;
		cout << "f_eq = " << f0 << endl;
		cout << "phi_n(f_eq) = " << C->phi_n(f0) << endl;
		cout << "Asym = " << Asym_n0 << endl;
		cout << "\n --------------------out func_fit----------------------" << endl;
	}
	double y0 = dEner_n0;
	double y1 = E_n0 - E0;
	double y2 = m_n*C->phi_n(f0) - m_eff_n0;
	double y3 = K - K0;
	double y4 = Asym_n0 - Asym0;

	gsl_vector_set(f, 0, y0);
	gsl_vector_set(f, 1, y1);
	gsl_vector_set(f, 2, y2);
	gsl_vector_set(f, 3, y3);
	gsl_vector_set(f, 4, y4);

	if(debug){
		printf("Equation: %f, %f, %f, %f, %f \n", y0,y1,y2,y3,y4);
	}
	return GSL_SUCCESS;

}



int fit(fit_eq_params p, set_const * C){
	const gsl_multiroot_fsolver_type * T;
	gsl_multiroot_fsolver *s;
	int status;
	size_t i, iter = 0;
	const size_t n = 5;
	fit_eq_params * params = &p;
	gsl_multiroot_function F = {&func_fit, n, params};

	double x_init[5] = {11.0, 5.0, 10.0, 0.0,0.0}; // MW(u) working params
	//double x_init[5] = {13.0, 7.0, 10.0, 1e-2, 1e-2}; // Phi(f) = 1/(1+f) working params
	//double x_init[5] = {13.4, 9.0, 10.0, 5e-3, 3e-4};
	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector_set(x, 0, x_init[0]);
	gsl_vector_set(x, 1, x_init[1]);
	gsl_vector_set(x, 2, x_init[2]);
	gsl_vector_set(x, 3, x_init[3]);
	gsl_vector_set(x, 4, x_init[4]);
	
	T = gsl_multiroot_fsolver_hybrid;
	s = gsl_multiroot_fsolver_alloc(T, n);
	gsl_multiroot_fsolver_set(s, &F, x);
//	print_state(iter, s);
	do{
		iter++;
		
		printf(" const    %f %f %f %f %f \n", gsl_vector_get(s->x, 0),
			gsl_vector_get(s->x, 1),
			gsl_vector_get(s->x, 2),
			gsl_vector_get(s->x, 3),
			gsl_vector_get(s->x, 4));

		status = gsl_multiroot_test_residual(s->f, 1e-3);
		status = gsl_multiroot_fsolver_iterate(s);
//		print_state(iter, s);

		if (status){
			break;
		}

	}
	while (status == GSL_CONTINUE && iter < 10000);

	cout << "status = " << gsl_strerror(status) << endl;

	double Cs = gsl_vector_get(s->x, 0);
	double Co = gsl_vector_get(s->x, 1);
	double Cr = gsl_vector_get(s->x, 2);
	double b = gsl_vector_get(s->x, 3);
	double c = gsl_vector_get(s->x, 4);


	C->C_s = Cs;
	C->C_o = Co;
	C->C_r = Cr;
	C->b = b;
	C->c = c;
	C->z = p.z;



	gsl_multiroot_fsolver_free(s);
	gsl_vector_free(x);
	return 0;
}
