
#include "Auxillary.h"
#include "set_const.h"
#include "EoS.h"
#include "constants.h"
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <iostream>

namespace Auxillary{

	double SD_func(double x, void * params){
		set_const * C = (set_const * ) params;
		return SNM_BED(x, C);
	}



	double SaturationDensity(set_const *C){
		const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
		gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);
		gsl_function F;
		bool debug = 0;
		double m, a, b, m_expected = 0.5;
		int status;
		int iter = 0, max_iter = 100;

		F.function = &SD_func;
		F.params = C;
		gsl_error_handler_t * def = gsl_set_error_handler_off();
		gsl_min_fminimizer_set(s, &F, 0.5, 0.1, 4.0);
		
		do {
			iter++;
			status = gsl_min_fminimizer_iterate(s);
			m = gsl_min_fminimizer_x_minimum(s);
			a = gsl_min_fminimizer_x_lower(s);
			b = gsl_min_fminimizer_x_upper(s);

			status = gsl_min_test_interval(a, b, 0.001, 0.0);

			if (status == GSL_SUCCESS)
				if (debug)
					printf("Converged:\n");

			if (debug)
				printf("%5d [%.7f, %.7f] "
				"%.7f %+.7f %.7f\n", iter, a, b, m, m - m_expected, b - a);
		} while (status == GSL_CONTINUE && iter < max_iter);
		if (debug)
		{
			std::cout << "DONE" << std::endl;
		}

		gsl_min_fminimizer_free(s);
		gsl_set_error_handler(def);
		if (status != 0){
			m = 2.0;
		}
		//cout << m << endl;
		//	cout << func_enermin(m, NULL) - m_n << endl;
		return m;
	}

	double SNM_BED(double n, set_const * C){
		return EoS::t_E(0.5*n, 0.5*n, C)/(D*n) - m_n;
	}
	double PNM_BED(double n, set_const * C){
		return EoS::t_E(n, 0.0, C)/(D*n) - m_n;
	}
	double NSM_BED(double n, set_const * C){
		return EoS::E(n, C)/(D*n) - m_n;
	}
	double Asymm(double rho, set_const * C){
		double drp = 1e-4;
		double f = EoS::f_eq(rho / 2, rho / 2, C);
		double d2Ep = (EoS::t_E(rho/2 - drp, rho/2 + drp, f, C) - 2*EoS::t_E(rho/2,rho/2, f, C) +
			EoS::t_E(rho/2 + drp, rho/2 - drp,f,  C))/(D*drp*drp);
		return rho*d2Ep/8.0;
	}
	double L(double n, set_const* C){
		double dn = 1e-4;
		double d_a = (Asymm(n-2*dn, C) - 8*Asymm(n-dn,C) + 8*Asymm(n+dn, C) - Asymm(n + 2*dn, C)) / (12*dn);
		return 3.0*n*d_a;
	}
	double eps_second(double x, set_const *C){
		double h = 1e-2;
		//return (E(x+h) - 2.0*E(x) + E(x-h))/(pow(h,2));
		return (SNM_BED(x+h, C) - 2.0*SNM_BED(x, C) + SNM_BED(x-h, C))/(pow(h,2));
	}
	double Compr(set_const* C){
		double rho_eq = n0;
		return 9.0*pow(rho_eq,2.0) * eps_second(rho_eq, C);
	}

}
