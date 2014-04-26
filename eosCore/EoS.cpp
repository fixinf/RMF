/*
* EoS.cpp
*	Provides EoS described by arXiv:nucl-th/041003 v1 Pt.2
*
*
*
*  Created on: 16.05.2013
*      Author: const.maslov@gmail.com
*/

#include <math.h>
#include "constants.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "set_const.h"
#include <iostream>
#include "EoS.h"
#include <gsl/gsl_integration.h>
#include <vector>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_multiroots.h>
#include <fstream>

using namespace std;

namespace EoS{
	double p_f(double n) {
		return pow(3.0 * pi * pi * D * n, 1.0 / 3.0);
	}

	double intf_func_f_eq(double p, void * params) {
		double m_eff = *(double *)params;
		return p * p / sqrt(m_eff * m_eff + p * p);
	}

	gsl_integration_workspace * f_eq_w = gsl_integration_workspace_alloc(1000);

	double func_f_eq(double f, void * params) {
		struct func_f_eq_params *p = (struct func_f_eq_params *) params;
		double nn = p->nn;
		double np = p->np;
		set_const * C = p->C;
		double dx = 1e-6;
		double res = (t_E(nn, np, f + dx, C) - t_E(nn, np, f, C)) / (dx);
		//printf("\n Cs = %f, Co = %f, b = %f, c = %f \n", C->C_s, C->C_o, C->b, C->c);
		//printf("f = %f, nn = %f, np = %f, res = %f \n", f, nn, np, res);
		return res;
	}



	double func_f_eq_multi(double f, void * params){
		struct func_f_eq_multi_params *p = (struct func_f_eq_multi_params *) params;
		set_const * C = p->C;
		std::vector<double> n = p->n;
		double dx = 1e-6;
		//return (t_E(n, f + dx, C) - t_E(n, f - dx, C)) / (2 * dx);
		double res = (t_E(n, f + dx, C) - t_E(n, f, C)) / (dx);
		//printf("res=%f \n", res);
		return res;
	}


	const gsl_root_fsolver_type *T_f_eq = gsl_root_fsolver_brent;
	gsl_root_fsolver *s_f_eq = gsl_root_fsolver_alloc(T_f_eq);

	int f_eq_solve(vec n, set_const *C, double lo, double hi, double * res){
		int status;
		int iter = 0, max_iter = 3000;

		double x = 0.0, x_expect = 0.5;
		double xmin = lo, xmax = hi;

		

		gsl_function F;
		F.function = &func_f_eq_multi;
		func_f_eq_multi_params params = { n, C };
		F.params = &params;

		//Like there is no root if function values have the same sign:
		double fmin = F.function(xmin, &params);
		double fmax = F.function(xmax, &params);

		if (fmin*fmax > 0){
			return -10;
		}

		gsl_root_fsolver_set(s_f_eq, &F, xmin, xmax);
		

		do {
			iter++;
			status = gsl_root_fsolver_iterate(s_f_eq);
			//printf("status = %i\n", status);
			x = gsl_root_fsolver_root(s_f_eq);
			xmin = gsl_root_fsolver_x_lower(s_f_eq);
			xmax = gsl_root_fsolver_x_upper(s_f_eq);
			status = gsl_root_test_interval(xmin, xmax, 0, 1e-9);
			//printf("status = %i\n", status);

		} while (status == GSL_CONTINUE && iter < max_iter);

		if (status != 0) {
			cout << "error, status = " << status << endl;
		}
		//printf("func_feq_vec: nn = %f, np = %f, res = %f \n", n[0], n[1], x);

		*res = x;
		return status;
	}

	double f_eq(vec n, set_const * C){
		int status = 1;
		double lo = 1e-5;
		double hi = 0.4;
		double trystep = 0.01;
		double res;
		gsl_error_handler_t * def = gsl_set_error_handler_off();

		while ((status != 0) && (hi < 1.0)){
		//	printf("hi = %f\n", hi);
			status = f_eq_solve(n, C, lo, hi, &res);
			if (status != 0){
				hi += trystep;
			}
		}

		gsl_set_error_handler(def);
		return res;
	}

	double f_eq(double nn, double np, set_const * C) {
		int status;
		int iter = 0, max_iter = 3000;

		double x = 0.0, x_expect = 0.5;
		double xmin = 1e-5, xmax = 1.0;

		gsl_function F;
		F.function = &func_f_eq;
		func_f_eq_params params = { nn, np, C };
		F.params = &params;

		//for (int i = 0; i < 100; i++){
		//	if (func_f_eq(0.01*i, &params) < 0){
		//		xmin = 0.01*i;
		//	}
		//}

		gsl_root_fsolver_set(s_f_eq, &F, xmin, xmax);

		do {
			iter++;
			status = gsl_root_fsolver_iterate(s_f_eq);
			x = gsl_root_fsolver_root(s_f_eq);
			xmin = gsl_root_fsolver_x_lower(s_f_eq);
			xmax = gsl_root_fsolver_x_upper(s_f_eq);
			status = gsl_root_test_interval(xmin, xmax, 0, 1e-9);

		} while (status == GSL_CONTINUE && iter < max_iter);

		if (status != 0) {
			cout << "error, status = " << status << endl;
		}
		//printf("func_feq_orig: nn = %f, np = %f, res = %f \n", nn, np, x);
		return x;
	}

	double mu_e(double n, double np, double f, set_const * C) {
		double res = pow(C->C_r, 2.0) * D *(n - 2.0 * np) / (2.0 * m_n *m_n* C->eta_r(f));
		res -= sqrt(pow(m_n * C->phi_n(f), 2.0) + pow(p_f(np), 2.0));
		res += sqrt(pow(m_n * C->phi_n(f), 2.0) + pow(p_f(n - np), 2.0));
		double dn = 1e-6;
		if (dn > np){
			dn = np / 2.0;
		}
		double mu_n = (t_E(n - np - 2 * dn, np, f, C) - 8 * t_E(n - np - dn, np, f, C) +
			8 * t_E(n - np + dn, np, f, C) - t_E(n - np + 2 * dn, np, f, C)) / (12 * D * dn);
		double mu_p = (t_E(n - np, np - 2 * dn, f, C) - 8 * t_E(n - np, np - dn, f, C) +
			8 * t_E(n - np, np + dn, f, C) - t_E(n - np, np + 2 * dn, f, C)) / (12 * D * dn);
		double res2 = mu_n - mu_p;
		return res;
	}

	double mu(int i, vec n, set_const * C){
		double f = f_eq(n, C);
		double dn = 1e-7;
		double E_p, E_m, E_s;
		//mu_n
		double mu;
		n[i] = n[i] + dn;
		E_p = t_E(n, f, C);
		n[i] = n[i] - dn;
		E_m = t_E(n, f, C);
		mu = (E_p - E_m) / (D * dn);
		return mu;
	}

	double func_f_eq_hyper(double f, void * params){
		struct func_f_eq_multi_params *p = (struct func_f_eq_multi_params *) params;
		set_const * C = p->C;
		std::vector<double> n_small = p->n;
		vec n;
		
		double dx = 1e-6;
		return (t_E(n, f + dx, C) - t_E(n, f - dx, C)) / (2 * dx);
	}

	double f_eq_hyper(vec n_small, set_const *C){
		int status;
		int iter = 0, max_iter = 3000;

		double x = 0.0, x_expect = 0.5;
		double xmin = -0.001, xmax = 1.0;

		gsl_function F;
		F.function = &func_f_eq_multi;
		func_f_eq_multi_params params = { n_small, C };
		F.params = &params;

		gsl_root_fsolver_set(s_f_eq, &F, xmin, xmax);

		do {
			iter++;
			status = gsl_root_fsolver_iterate(s_f_eq);
			x = gsl_root_fsolver_root(s_f_eq);
			xmin = gsl_root_fsolver_x_lower(s_f_eq);
			xmax = gsl_root_fsolver_x_upper(s_f_eq);
			status = gsl_root_test_interval(xmin, xmax, 0, 1e-9);

		} while (status == GSL_CONTINUE && iter < max_iter);

		if (status != 0) {
			cout << "error, status = " << status << endl;
		}
		return x;
	}

	const gsl_root_fsolver_type *T_np_eq = gsl_root_fsolver_brent;
	gsl_root_fsolver *s_np_eq = gsl_root_fsolver_alloc(T_np_eq);

	double func_np(double np, void * params) {
		struct func_np_params *p = (struct func_np_params *) params;
		double n = p->ntot;


		set_const * C = p->C;
		vec n_ext = p->n_ext;
		double sum = 0.0;
		for (int i = 0; i < n_ext.size(); i++){
			sum += n_ext[i];
		}

		double sum_ch = 0.0;
		for (int i = 0; i < n_ext.size(); i++){
			sum_ch += n_ext[i]*C->Q[i+2];
		}

		//n = n - sum;

		vec _n;
		_n.push_back(n - np);
		_n.push_back(np);
	
		for (int i = 0; i < n_ext.size(); i++){
			_n.push_back(n_ext[i]);
		}


		double f = f_eq(_n, C);
		//double mue = mu_e(n, np, f, C);
		double mu_n = mu(0, _n, C);
		double mu_p = mu(1, _n, C);
		double mue = mu_n - mu_p;
	//	printf("%f \n", mue);

		double mue2 = mu_e(n, np, f, C);
		//mue = mue2;
		
		//cout << "MUE = " << mue << " N = "<< n << " NP = " << np << " F = " << f << endl;
		double result = np;
		//cout << "RESULT = " << result << endl;
		double n_l = 0.0;
		if (mue*mue - m_e*m_e >= 0){
			n_l += pow(mue*mue - m_e*m_e, 3.0 / 2.0) / (3.0*D*pi*pi);
		}
		//cout << "RESULT = " << result << endl;
		if (mue*mue - m_mu*m_mu >= 0){
			n_l += pow(mue*mue - m_mu*m_mu, 3.0 / 2.0) / (3.0*D*pi*pi);
		}
		result -= n_l - sum_ch;
		/*printf("np_eq: mue = %f, n_ch = %f \n",
			mue, sum_ch);*/
		//cout << "RESULT = " << result << endl;
	//	printf("ntot = %f, result = %f \n", n, result);
		return result;
	}


	int np_eq_solve(double n, vec n_ext, double xmin, double xmax, set_const * C, double & res){
		int status;
		int iter = 0, max_iter = 1000;
		double x = 0.0, x_expect = (xmax - xmin)/2.0;
		gsl_function F;
		F.function = &func_np;
		func_np_params par = { n, C, n_ext };
		F.params = &par;

		status = gsl_root_fsolver_set(s_np_eq, &F, xmin, xmax);
		if (status != 0){
			return status;
		}

		do {
			iter++;
			status = gsl_root_fsolver_iterate(s_np_eq);
			x = gsl_root_fsolver_root(s_np_eq);
			xmin = gsl_root_fsolver_x_lower(s_np_eq);
			xmax = gsl_root_fsolver_x_upper(s_np_eq);
			status = gsl_root_test_interval(xmin, xmax, 0, 1e-10);
		} while (status == GSL_CONTINUE && iter < max_iter);
		;
		if (status != 0) {
			cout << "error, status = " << status << endl;
		}
		res = x;
		return status;
	}
	double np_eq(double n, vec n_ext, set_const * C) {
		double res;
		int status;
		gsl_error_handler_t * def = gsl_set_error_handler_off();
		double sum = 0.0;
		for (int i = 0; i < n_ext.size(); i++){
			sum += n_ext[i];
		}
		double _n = n;
		double xmin = 0.0, xmax = 0.2 * _n;
		status = 1;
		while ((status != 0) && (xmax < n)){
			status = np_eq_solve(n, n_ext, xmin, xmax, C, res);
			xmax += 0.1*_n;
		}
		if (status != 0){
			printf("! \n");
// 			std::ofstream ofs("np.dat");
// 			func_np_params params = {n, C, n_ext};
// 			for (int i = 0; i < 100; i++){
// 				ofs << n*i/100 << "," << func_np(n*i/100, &params) << endl;
// 			}
// 			ofs.close();

			res = 0.95*_n;
		}
		//printf("np_eq(%f, %f), ntot = %f, status = %i, res = %f, \n",
		//	n_ext[0], n_ext[4], n, status, res);
		//printf("np_eq status = %i", status);
		gsl_set_error_handler(def);
		if (res < 0){
			res = 0.0;
		}
		return res;
	}
//
//	int np_eq(double ntot, vec & result, set_const * C){
//		const gsl_multiroot_fsolver_type *T;
//		gsl_multiroot_fsolver *s;
//
//		int status;
//		size_t i, iter = 0;
//
//		const size_t num = 7;
//		struct func_np_params p;
//		gsl_multiroot_function f = {&func_np, num , &p };
//
//		p.C = C;
//		p.ntot = ntot;
//
//		double x_init[8] = {/*93 * ntot / 100,*/ 1.0*ntot / 100.0, ntot / 100, ntot / 100, ntot / 100,
//			ntot / 100, ntot / 100, ntot / 100};
//		gsl_vector *x = gsl_vector_alloc(num);
//
//		for (int i = 0; i < num; i++){
//			gsl_vector_set(x, i, x_init[i]);
//		}
//
//		T = gsl_multiroot_fsolver_hybrid;
//		s = gsl_multiroot_fsolver_alloc(T, num);
//		gsl_multiroot_fsolver_set(s, &f, x);
//
////		print_state(iter, s);
//
//		do
//		{
//			iter++;
//			status = gsl_multiroot_fsolver_iterate(s);
//
//	//		print_state(iter, s);
//		//	printf("DX = %.46f", s->dx);
//
//			if (status)   /* check if solver is stuck */
//				break;
//			
//			status =
//				gsl_multiroot_test_residual(s->f, 1e-7);
//		} while (status == GSL_CONTINUE && iter < 1000);
//
//		printf("status = %s\n", gsl_strerror(status));
//
//		result.resize(num);
//		for (i = 0; i < num; i++){
//			result[i] = abs(gsl_vector_get(s->x, i));
//		}
//
//		gsl_multiroot_fsolver_free(s);
//		gsl_vector_free(x);
//		return 0;
//	}

	double func_e(double p, void * params) {
		double m_eff = *(double *)params;
		return p * p * sqrt(p * p + m_eff * m_eff);
	}


	gsl_integration_workspace * w_t_E = gsl_integration_workspace_alloc(1000);

	double KineticIntegral(double p_f, double m, set_const *C){
		double m_eff = m;

		double result, error;
		gsl_function F;
	//	cout << "PF = " << p_f << endl;
		F.function = &func_e;

		F.params = &m_eff;
		gsl_integration_qags(&F, 0.0, p_f, 0, 1e-10, 1000, w_t_E, &result, &error);
		result = result / (pi*pi);
		return result;
// 		double result2 = 2*p_f*sqrt(m_eff*m_eff + p_f*p_f)*(m_eff*m_eff + 2*p_f*p_f);
// 		if (m_eff != 0.0){
// 			result2 += pow(m_eff,4)*log(m_eff*m_eff);
// 			result2 -= 2*pow(m_eff,4)*log(p_f + sqrt(m_eff*m_eff + p_f*p_f));
// 		}
// 		result2 = result2/(16*pi*pi);
//		printf("p_f = %f, m = %f, 1 = %f, 2 = %f \n", p_f, m_eff, result, result2);
//		return result2;
	}

	double mu_e_2(double n, double np, double f, set_const * C) {
		double res = pow(C->C_r, 2.0) * D *(n - 2.0 * np) / (2.0 * m_n *m_n* C->eta_r(f));
		//cout << res << endl;
		res -= sqrt(pow(m_n * C->phi_n(f), 2.0) + pow(p_f(np), 2.0));
		//cout << res << endl;
		res += sqrt(pow(m_n * C->phi_n(f), 2.0) + pow(p_f(n - np), 2.0));
		//cout << res << endl;
		return res;
	}

	double t_E(double nn, double np, set_const * C) {
		double f = f_eq(nn, np, C);
		double res = t_E(nn, np, f, C);

		double me = m_e;

		double pf_e = 0;
		if (pow(mu_e(nn + np, np, f, C), 2.0) - me*me >= 0){
			pf_e = sqrt(pow(mu_e(nn + np, np, f, C), 2.0) - me*me);
		}

		if (np != 0.0){
			res += KineticIntegral(pf_e, m_e, C);
		}

		double mmu = m_mu;
		double pf_mu = 0;
		if (pow(mu_e(nn + np, np, f, C), 2.0) - mmu*mmu >= 0){
			pf_mu = sqrt(pow(mu_e(nn + np, np, f, C), 2.0) - mmu*mmu);
		}

		if (np != 0.0){
			res += KineticIntegral(pf_mu, m_mu, C);
		}

		return res;

	}

	double t_E(double nn, double np, double f, set_const * C) {
		//double res = 0.5 * pow(m_n * m_n * f / C->C_s, 2.0);

		//res += C->U(f);

		//res += KineticIntegral(p_f(np), m_n*C->phi_n(f), C);

		//res += KineticIntegral(p_f(nn), m_n*C->phi_n(f), C);

		//res += 0.5 * pow(C->C_o * D * (nn + np) / m_n, 2.0) / (C->eta_o(f));

		//res += pow(C->C_r * D * (nn - np) / m_n, 2.0) / (8.0 * C->eta_r(f));
		//
		//return res;
		vec n;
		n.push_back(nn);
		n.push_back(np);
		double res2 = t_E(n, f, C);

		return res2;
	}

	double t_E(std::vector<double> n, double f, set_const * C) {
		double res = 0.5 * pow(m_n * m_n * f / C->C_s, 2.0)*C->eta_s(f);

		res += C->U(f);

		int num = n.size();
		/*std::vector<double> nvec;
		for (int i = 0; i < num; i++){
			nvec.push_back(boost::python::extract<double>(n[i]));
		}*/

		double o_sum = 0;
		double r_sum = 0;

		for (int i = 0; i < num; i++){
//			cout << i << ", n[i] = " << n[i] << ", pf = " <<  p_f(n[i]) << endl;
			//res += KineticIntegral(p_f(n[i]), C->M[i] * C->phi_n(C->X_s[i] * f), C);
			res += KineticIntegral(p_f(n[i]), C->M[i] * C->phi_n(C->X_s[i] * (C->M[0]/C->M[i]) * f), C);
			o_sum += C->X_o[i] * n[i];
			r_sum += C->X_r[i] * C->T[i] * n[i];
		}


		res += 0.5 * pow(C->C_o * D * (o_sum) / m_n, 2.0) / (C->eta_o(f));

		res += 0.5 * pow(C->C_r * D * (r_sum) / m_n, 2.0) / (C->eta_r(f));


		//double me = m_e;

		//double pf_e = 0;
		//if (pow(mu_e(n[0] + n[1], n[1], f, C), 2.0) - me*me >= 0){
		//	pf_e = sqrt(pow(mu_e(n[0] + n[1], n[1], f, C), 2.0) - me*me);
		//}

		//if (n[1] != 0.0){
		//	res += KineticIntegral(pf_e, m_e, C);
		//}

		//double mmu = m_mu;
		//double pf_mu = 0;
		//if (pow(mu_e(n[0] + n[1], n[1], f, C), 2.0) - mmu*mmu >= 0){
		//	pf_mu = sqrt(pow(mu_e(n[0] + n[1], n[1], f, C), 2.0) - mmu*mmu);
		//}

		//if (n[1] != 0.0){
		//	res += KineticIntegral(pf_mu, m_mu, C);
		//}

		//gsl_function F;
		//F.function = &func_e;
		//double result = 0.0;
		//double error;
		//double me = m_e;
		//F.params = &me;
		//double pf_e = 0;
		//if (pow(mu_e_2(n[0] + n[1], n[1], f, C), 2.0) - me*me >= 0){
		//	pf_e = sqrt(pow(mu_e_2(n[0] + n[1], n[1], f, C), 2.0) - me*me);
		//}
		////	cout << "MU_E = " <<  mu_e(nn+np, np, f, C) << endl;
		//gsl_integration_qags(&F, 0.0, pf_e, 0, 1e-10, 1000, w_t_E, &result, &error);

		//if (n[1] != 0.0){
		//	//cout << "hey! " << np << endl;
		//	res += result / (pi*pi);
		//}
		////cout << "RESULT   " << result << endl;

		//double mmu = m_mu;
		//F.params = &mmu;
		//double pf_mu = 0;
		//if (pow(mu_e_2(n[0] + n[1], n[1], f, C), 2.0) - mmu*mmu >= 0){
		//	pf_mu = sqrt(pow(mu_e_2(n[0] + n[1], n[1], f, C), 2.0) - mmu*mmu);
		//}
		//gsl_integration_qags(&F, 0.0, pf_mu, 0, 1e-10, 1000, w_t_E, &result, &error);


		//if (n[1] != 0.0){
		//	//cout << "hey! " << np << endl;
		//	res += result / (pi*pi);
		//}

		return res;
	}



	double t_P(double nn, double np, set_const * C) {
		return 0.0;
	}

	double E(double n, set_const * C){
		vec n_ext;
 		double np = np_eq(n, n_ext, C);
		vec _n;
		_n.push_back(n - np);
		_n.push_back(np);
		double f = f_eq(_n, C);
		//printf("n = %.18f, np = %.18f, feq = %.18f \n",n , np, f);
 		//return t_E(n - np, np, f, C);
		return t_E(n - np, np, f, C);
		//return 0;
	}


	double P(double n, set_const * C) {
		double dn = 1e-5;
		double d_E = (E(n + dn, C) - E(n - dn, C)) / (2.0*dn);
		return n*d_E - E(n, C);
	}
};
