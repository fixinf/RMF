
#include "star.h"
#include "set_const.h"
#include "eq_volkov.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include "EoS.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include "constants.h"


int star(double rho_init, double result[], set_const* C) {
	gsl_odeiv2_system sys = { eq_volkov, NULL, 3, NULL };
	double delta = 5e-3;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
		gsl_odeiv2_step_rkf45, 1e-3, 1e-3, 1e-3);
	int i;
	double r_init = 1e-6;
	double t = r_init, t1 = 15.0 + r_init;
	double P_init = EoS::P(rho_init, C);
	double E_init = EoS::E(rho_init, C);
	double y[3] = { EoS::P(rho_init, C), 1.333 * M_PI * pow(r_init, 3.0) * 1.4793,
		E_init };
	double f[3];
	int status = eq_volkov(r_init, y, f, NULL);
	for (i = 1; i <= 10000; i++) {
		printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
		double ti = i * t1 / 1000.0;
		if (y[0] > delta * P_init) {
			status = gsl_odeiv2_driver_apply(d, &t, ti, y);

			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n", status);
				break;
			}
		} else {
			std::cout << "RES2 " << y[1] << "      " << t << std::endl;
			result[0] = y[1];
			result[1] = t;
			gsl_odeiv2_driver_free(d);
			return 0;
		}
		y[2] = EofP(y[0], C);
	}
	return 0;
}

int star(double rho_init, double result[], EoS_tab EoS) {
	gsl_odeiv2_system sys = { eq_volkov, NULL, 3, NULL };
	double delta = 5e-3;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys,
		gsl_odeiv2_step_rkf45, 1e-2, 1e-2, 1e-2);
	int i;
	double bMass = 0.0;
	double r_init = 1e-6;
	double t = r_init, t1 = 15.0 + r_init;
	double P_init = EoS.P(rho_init);
	double E_init = EoS.E(rho_init);
	double y[3] = { EoS.P(rho_init), 1.333 * M_PI * pow(r_init, 3.0) * 1.4793,
		E_init };
	double f[3];
	int status = eq_volkov(r_init, y, f, NULL);
	for (i = 1; i <= 10000; i++) {
		//printf("%f %f %f %f \n \r", t, y[0], y[1], y[2]);
		double ti = i * t1 / 1000.0;
		if (y[0] > delta * P_init) {
			status = gsl_odeiv2_driver_apply(d, &t, ti, y);

			if (status != GSL_SUCCESS) {
				printf("error, return value=%d\n", status);
				break;
			}
		} else {
			//cout << "RES2 " << y[1] << "      " << t << endl;
			result[0] = y[1];
			result[1] = t;
			result[2] = bMass;
			gsl_odeiv2_driver_free(d);
			return 0;
		}
		y[2] = EoS.EofP(y[0]);
		bMass += (1.479345850e-12)*D*m_n*EoS.NofE(y[2])*pow((1.0 - (2.0*y[1])/t),-0.5)*t*t*(t1/1000.0);
	}
	return 0;
}
