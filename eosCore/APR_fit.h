#pragma once

#ifndef _APR_FIT_H
#define _APR_FIT_H

#include "set_const.h"
#include <gsl/gsl_multifit_nlin.h>
#include "EoS.h"
#include <iostream>



namespace APR_fit{


	struct data {
		size_t n;
		double * y;
		double * t;
		set_const* Init;
	};

	/*double k_fit = 0.5/0.16;
	double t_snm[13] = {k_fit*0.04, k_fit*0.08, k_fit*0.12, k_fit*0.16,
	k_fit * 0.2, k_fit*0.24, k_fit*0.32, k_fit * 0.4,
	k_fit*0.48, k_fit*0.56, k_fit*0.64, k_fit*0.8,
	k_fit*0.96};
	double t_pnm[14] = {k_fit*0.02,k_fit*0.04, k_fit*0.08, k_fit*0.12, k_fit*0.16,
		k_fit * 0.2, k_fit*0.24, k_fit*0.32, k_fit * 0.4,
		k_fit*0.48, k_fit*0.56, k_fit*0.64, k_fit*0.8,
		k_fit*0.96};
	
	double y_snm[13] = {-6.48, -12.13, -15.04, -16.00,
	-15.09, -12.88, -5.03, -2.13, 
	-15.46, -34.39, -58.35, -121.25, -204.02};

	double y_pnm[14] = {4.45, 6.45, 9.65, 13.29, 17.94, 
	22.92, 27.49, 38.82, 54.95, 75.13, 99.75,
	127.58, 205.34, 305.87};

	void CompareSNM(set_const * C){
		for (int i = 0; i < 13; i++){
			std::cout << t_snm[i] << "   " << EoS::t_E(t_snm[i]/2,t_snm[i]/2,C) << "   " 
				<< y_snm[i] << "    " << EoS::t_E(t_snm[i]/2,t_snm[i]/2,C) - y_snm[i] << std::endl;
		}
	}

	void ComparePNM(set_const * C){
		for (int i = 0; i < 14; i++){
			std::cout << t_pnm[i] << "   " << EoS::t_E(t_pnm[i],0.0,C) << "   " 
				<< y_pnm[i] << "    " << EoS::t_E(t_pnm[i],0.0,C) - y_pnm[i] << std::endl;
		}
	}*/

	double fit(set_const*, double);
	double fit_n(set_const*, double);
	void cond(set_const* C, double Cs, double Co, double b, double c);
	void uncond(set_const* C, double Cs, double Co, double b, double c);
	void cond(set_const* C, double Cs, double Co, double b, double c, double Cr);
	void cond(set_const* C, double Cs, double Co, double b, double c, double Cr);
	void uncond(set_const* C);
	void print_state (size_t iter, gsl_multifit_fdfsolver * s);
}
#endif