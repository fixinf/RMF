
#include "APR_fit.h"

#include "EoS.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include "set_const.h"
#include <functional>
#include "constants.h"
#include <complex>
#include <math.h>
#include "MixedFit.h"

namespace APR_fit{

	int func_fit_mix(const gsl_vector * x, void * params, gsl_vector * f){
		size_t n_snm = ((struct data_mix *)params)->n_snm;
		size_t n_pnm = ((struct data_mix *)params)->n_pnm;
		size_t n = n_snm + n_pnm;
		double *y = ((struct data_mix *)params)->y;
		double *t = ((struct data_mix *) params)->t;
		set_const* C_init = ((struct data_mix *) params)-> Init;
		double Cs = gsl_vector_get (x, 0);
		double Co = gsl_vector_get (x, 1);
		double Cr = gsl_vector_get(x,4);
		double b = gsl_vector_get (x, 2);
		double c = gsl_vector_get (x, 3);
		//double z = gsl_vector_get(x,5);
		//double a = gsl_vector_get(x,6);
		double z = 0.65;
		set_const *C = C_init;
		C->name = "Equation dummy constant set";
		cond(C, Cs, Co, b, c, Cr);
		C->z = z;
		C->a = 1.0;
	
		/*set_const C("Equation dummy constant set",abs(Cs),Co,sqrt(100.0),b,abs(c),z,
		[](double f){return 1-f;},
		[](double f){return 1.0;},
		[=](double f){return eta_o(f);},
		[](double f){return 1.0;});*/

	/*	set_const C("Equation dummy constant set",abs(Cs),abs(Co),sqrt(100.0),abs(b),abs(c),0.65,
			[](double f){return 1-f;},
			[](double f){return 1.0;},
			[](double f){return 1.0;},
			[](double f){return 1.0;});
*/
		size_t i;

		for (i = 0; i < n_snm; i++)
		{
			/* Model Yi = A * exp(-lambda * i) + b */
			double Yi = EoS::t_E(t[i]/2,t[i]/2, C)/(D*t[i]) - m_n;
			//std::cout << i << "   " <<  Yi - y[i]<< std::endl;
			gsl_vector_set (f, i, abs(Yi - y[i]));//Try abs(Yi-y[i])
		}

		for (i = n_snm; i < n; i++)
		{
			/* Model Yi = A * exp(-lambda * i) + b */
			double Yi = EoS::t_E(t[i],0.0, C)/(D*t[i]) - m_n;
			//std::cout << i << "   " <<  Yi - y[i]<< std::endl;
			gsl_vector_set (f, i, abs(Yi - y[i]));//Try abs(Yi-y[i])
		}



		//gsl_vector_set(f, n-1, abs(asymm(0.5, C) - 32.0));
		//std::cout << "ASYMM = " << asymm(0.5,C) << std::endl;

		return GSL_SUCCESS;
	}




	double fit_mix(set_const* Init, int nn){
		const gsl_multifit_fdfsolver_type *T;
		gsl_multifit_fdfsolver *s;
		int status;
		unsigned int i, iter = 0;
		const size_t n_snm = 7;
		const size_t n_pnm = nn;
		const size_t n = n_snm + n_pnm;
		const size_t p = 5;
		double k = 0.5/0.16;
		gsl_matrix *covar = gsl_matrix_alloc (p, p);
		double y[15] =
		{/*Symm -> */ -6.48,-12.13,-15.04, -16.0, -15.09, -12.88, -5.03,
		/*Neutr -> */ 75.13, 99.75, 127.58, 205.34, 305.87}; 
		//double t[11] = {k*0.04, k*0.08,k*0.12,k*0.16,k*0.2,k*0.24, k*0.32, k*0.4,k*0.48, k*0.56, k*0.64};
		double t[15] = 
		{/*Symm ->*/k*0.04, k*0.08,k*0.12,k*0.16, k*0.20,k*0.24, k*0.32,
		/*Neutr -> */k*0.48,k*0.56,k*0.64, k*0.8, k*0.96};

		//^_^ Good, but not enough

		//double y[15] =
		//{/*Symm -> */ -6.48,-12.13,-15.04, -16.0, -15.09, -12.88, -5.03,
		///*Neutr -> */ 127.58,205.34,305.87};
		////double t[11] = {k*0.04, k*0.08,k*0.12,k*0.16,k*0.2,k*0.24, k*0.32, k*0.4,k*0.48, k*0.56, k*0.64};
		//double t[15] = 
		//{/*Symm ->*/k*0.04, k*0.08,k*0.12,k*0.16, k*0.20,k*0.24, k*0.32,
		///*Neutr -> */k*0.64,k*0.8,k*0.96};


		struct data_mix d = { n_snm, n_pnm, y, t, Init};
		gsl_multifit_function_fdf f;
		double x_init[5] = {Init->C_s,Init->C_o, Init->b,Init->c, Init->C_r};

		//double x_init[6]  = {11.56279437,7.49931859,0.00871711,0.00267620,0.86859184,0.5};
		//double x_init[4] = { sqrt(130.746),sqrt(120.7244),1.0,10.0};
		gsl_vector_view x = gsl_vector_view_array (x_init, p);
		const gsl_rng_type * type;
		gsl_rng * r;

		gsl_rng_env_setup();

		type = gsl_rng_default;
		r = gsl_rng_alloc (type);

		f.f = &func_fit_mix;
		f.df = NULL;
		f.fdf = NULL;
		f.n = n;
		f.p = p;
		f.params = &d;

		/* This is the data to be fitted */

		/*for (i = 0; i < n; i++)
		{
			double t = i;
			y[i] = 1.0 + 5 * exp (-0.1 * t) 
				+ gsl_ran_gaussian (r, 0.1);
			sigma[i] = 0.1;
			printf ("data: %u %g %g\n", i, y[i], sigma[i]);
		};*/

		T = gsl_multifit_fdfsolver_lmsder;
		
		s = gsl_multifit_fdfsolver_alloc (T, n, p);

		gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	
		do
		{
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);

			//printf ("status = %s\n", gsl_strerror (status));

//			print_state (iter, s);

			if (status)
				break;

			status = gsl_multifit_test_delta (s->dx, s->x,
				1e-8, 1e-8);
		}
		while (status == GSL_CONTINUE && iter < 2000);

		gsl_multifit_covar (s->J, 0.0, covar);

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))


		cond(Init, FIT(0), FIT(1), FIT(2), FIT(3), FIT(4));

		{ 
			double chi = gsl_blas_dnrm2(s->f);
			double dof = n - p;
			double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
			//double c = 1.0;
		/*	printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);

			printf ("Cs      = %.5f +/- %.5f\n", Init->C_s, c*ERR(0));
			printf ("Co = %.5f +/- %.5f\n", Init->C_o, c*ERR(1));
			printf ("b      = %.5f +/- %.5f\n", Init->c, c*ERR(2));
			printf ("c      = %.5f +/- %.5f\n", Init->b, c*ERR(3));
			printf ("Cr      = %.5f +/- %.5f\n", Init->C_r, c*ERR(4));*/
		}
		
		//printf ("status = %s\n", gsl_strerror (status));
		double z = 0.65;

		
		gsl_matrix_free (covar);
		gsl_rng_free (r);

		double yi = 0;
		/*for (int i = 0; i < n_snm; i++){
		double yi = EoS::t_E(t[i]/2,t[i]/2, Init)/(D*t[i]) - m_n ;
		printf("n = %.3f, %.3f  %.3f  %.3f \n",
		t[i],
		yi,
		y[i],
		yi-y[i]);

		}

		for (int i = n_snm; i < n; i++){
		double yi = EoS::t_E(t[i],0.0, Init)/(D*t[i]) - m_n ;
		printf("n = %.3f, %.3f  %.3f  %.3f \n",
		t[i],
		yi,
		y[i],
		yi-y[i]);

		}*/
		/*return *(new set_const("APR_fit return constant set",FIT(0), FIT(1), 10.0, FIT(2),abs(FIT(3)), z, 
			[](double f){return (1-f);},
			[](double f){return 1.0;},
			[=](double f){return eta_o(f);},
			[](double f){return 1.0;}));*/
		double rr = gsl_blas_dnrm2(s->x);
		gsl_multifit_fdfsolver_free (s);
		return rr;
	}

	

}
