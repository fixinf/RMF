
#include "EoS.h"
#include <gsl/gsl_multiroots.h>
#include "constants.h"
#include "hyper.h"
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_roots.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <fstream>


using namespace EoS;
namespace Hyper_EoS{
	struct func_l0_params{
		double ntot;
		set_const * C;
		vec n_ext;
	};


	double func_l0(double x, void * params){
		func_l0_params * p = (func_l0_params *)params;
		set_const * C = p->C;
		double ntot = p->ntot;
		vec n_ext = p->n_ext;
		vec n_in;
		n_in.push_back(x);
		for (int i = 0; i < n_ext.size(); i++){
			n_in.push_back(n_ext[i]);
		}
		//��������� �������� ��� �������� ��������� �����

		double np = np_eq(ntot - x, n_in, C);
		double sum = np;
		for (int i = 0; i < n_in.size(); i++){
			sum += n_in[i];
		}
		double nn = ntot - np - x;
		if (nn < 0.0){
			nn = 1e-6;
		}

		//Chemical potentials

		vec n;
		n.push_back(nn);
		n.push_back(np);
		n.push_back(x);
		for (int i = 0; i < n_ext.size(); i++){
			n.push_back(n_ext[i]);
		}
		double f = f_eq(n, C);

		//mu_n
		double mu_n = mu(0, n, C);

		//mu_p
		double mu_p = mu(1, n, C);

		//mu_p
		double mu_l0 = mu(2, n, C);

		//printf("nn = %f, np = %f, nL0 = %f \n", n[0], n[1], n[2]);
		return mu_l0 - mu_n;
	}


	const gsl_root_fsolver_type *T_l0_eq = gsl_root_fsolver_brent;
	gsl_root_fsolver *s_l0_eq = gsl_root_fsolver_alloc(T_l0_eq);

	int n_L0_solve(double ntot, vec n_ext, double xmin, double xmax, set_const * C, vec & out){
		int status;
		int iter = 0, max_iter = 300;
		double x = 0.0, x_expect = (xmax - xmin) / 2.0;
		//double xmin = 0.0, xmax = 0.3 * ntot;
		gsl_function F;
		F.function = &func_l0;
		func_l0_params par = { ntot, C, n_ext };
		F.params = &par;

		status = gsl_root_fsolver_set(s_l0_eq, &F, xmin, xmax);
		if (status != 0){
			return status;
		}
		do {
			iter++;
			status = gsl_root_fsolver_iterate(s_l0_eq);
			x = gsl_root_fsolver_root(s_l0_eq);
			xmin = gsl_root_fsolver_x_lower(s_l0_eq);
			xmax = gsl_root_fsolver_x_upper(s_l0_eq);
			status = gsl_root_test_interval(xmin, xmax, 1e-6, 1e-6);
		} while (status == GSL_CONTINUE && iter < max_iter);
		;
		if (status != 0) {
			std::cout << "error, status = " << status << std::endl;
		}
		//vec out;
		vec _in;
		_in.push_back(x);
		for (int i = 0; i < n_ext.size(); i++){
			_in.push_back(n_ext[i]);
		}
		double np = EoS::np_eq(ntot - x, _in, C);
		double nn = ntot - x - np;
		out.push_back(nn);
		out.push_back(np);
		out.push_back(x);
		return status;
	}

	vec n_L0(double ntot, vec n_ext, set_const * C){
		gsl_error_handler_t * def = gsl_set_error_handler_off();
		double sum = 0.0;
		for (int i = 0; i < n_ext.size(); i++){
			sum += n_ext[i];
		}
		ntot = ntot - sum;
		//Check if is possible to make hyperons
		double np = np_eq(ntot, n_ext, C);
		double nn = ntot - np;
		vec in_temp;
		in_temp.push_back(nn);
		in_temp.push_back(np);
		vec in2;
		in2.push_back(ntot);
		double f = f_eq(in2, C);
		double f2 = f_eq(in_temp, C);
// 		printf("%f %f \n", mu(0, in2, C), C->M[2]);
// 		printf("%f %f \n", mu(0, in_temp, C), C->M[2]);
		if (mu(0, in_temp, C)  < C->M[2]){
			vec res;
			res.push_back(nn);
			res.push_back(np);
			res.push_back(0.0);
			return res;
		}

	
		double xmax = 0.1*ntot;
		double xmin = 0.0;
		vec res;
		int status = 1;
		while ((status != 0) && xmax < 0.95*ntot){
			status = n_L0_solve(ntot, n_ext, xmin, xmax, C, res);
			xmax += 0.1*ntot;
		}
		if (status != 0){
			vec in;
			in.push_back(0.0);
			res.clear();
			for (int i = 0; i < n_ext.size(); i++){
				in.push_back(n_ext[i]);
			}
			double np = np_eq(ntot, in, C);
			res.push_back(ntot - np);
			res.push_back(np);
			res.push_back(0.0);
		}
		//printf("%i \n", status);
		gsl_set_error_handler(def);
		/*printf("mu_n = %f, mu_L0 = %f \n", mu(0, res, C), mu(2, res, C));*/
		return res;
	}


	double func_Sm(double x, void * params){
		func_l0_params * p = (func_l0_params *)params;
		set_const * C = p->C;
		double ntot = p->ntot;
		vec n_in;
		n_in.push_back(x);
		for (int i = 0; i < p->n_ext.size(); i++){
			n_in.push_back(p->n_ext[i]);
		}

		//��������� �������� ��� �������� ��������� �����

		vec res_L0 = n_L0(ntot, n_in, C);
		double nn = res_L0[0];
		double np = res_L0[1];
		//Chemical potentials

		vec n;
		for (int i = 0; i < res_L0.size(); i++){
			n.push_back(res_L0[i]);
		}
		n.push_back(x);
		for (int i = 0; i < p->n_ext.size(); i++){
			n.push_back(p->n_ext[i]);
		}

		double f = f_eq(n, C);

		double dn = 1e-7;
		double E_p, E_m;
		//mu_n
		double mu_n = mu(0, n, C);
		//mu_p
		double mu_p = mu(1, n, C);

		//mu_p
		double mu_Sm = mu(3, n, C);

		//printf("I'm working, don't kill me! \n");

		return mu_Sm - 2*mu_n + mu_p;
	}


	const gsl_root_fsolver_type *T_sm_eq = gsl_root_fsolver_brent;
	gsl_root_fsolver *s_sm_eq = gsl_root_fsolver_alloc(T_sm_eq);

	int n_Sm_solve(double ntot, vec n_ext, double xmin, double xmax, set_const * C, vec & out){
		int status;
		int iter = 0, max_iter = 300;
		double x = 0.0, x_expect = (xmax - xmin) / 2.0;
		//double xmin = 0.0, xmax = 0.3 * ntot;
		gsl_function F;
		F.function = &func_Sm;
		func_l0_params par = { ntot, C, n_ext };
		F.params = &par;

		//test output
		// 		std::ofstream ofs("func_sm.dat");
		// 		double step = (xmax-xmin)/100;
		// 		for (int i = 0; i < 100; i++){
		// 			ofs << step*i << "   " << func_Sm(step*i, &par) << std::endl;
		// 		}
		// 		ofs.close();
		//end test


		status = gsl_root_fsolver_set(s_sm_eq, &F, xmin, xmax);
		if (status != 0){
			return status;
		}
		do {
			iter++;
			status = gsl_root_fsolver_iterate(s_sm_eq);
			x = gsl_root_fsolver_root(s_sm_eq);
			xmin = gsl_root_fsolver_x_lower(s_sm_eq);
			xmax = gsl_root_fsolver_x_upper(s_sm_eq);
			status = gsl_root_test_interval(xmin, xmax, 1e-6, 1e-6);
		} while (status == GSL_CONTINUE && iter < max_iter);

		if (status != 0) {
			std::cout << "error, status = " << status << std::endl;
		}
		//vec out;
		vec in;
		in.push_back(0.0);
		in.push_back(0.0);
		in.push_back(0.0);
		in.push_back(x); 
		vec L0;
		L0 = n_L0(ntot, in, C);
		out.push_back(L0[0]);
		out.push_back(L0[1]);
		out.push_back(L0[2]);
		out.push_back(x);
		return status;
	}

	vec n_Sm(double ntot, vec n_ext, set_const * C){
		gsl_error_handler_t * def = gsl_set_error_handler_off();
		double sum = 0;
		for (int i = 0; i < n_ext.size(); i++){
			sum += n_ext[i];
		}
		ntot = ntot - sum;
		double xmax = 0.05*ntot;
		double xmin = 0.0;
		vec res;
		int status = 1;
		vec in_temp;
		vec in1;
		in1.push_back(0.0);
		for (int i = 0; i < n_ext.size(); i++){
			in1.push_back(n_ext[i]);
		}
		vec l0 = n_L0(ntot, in1, C);
		vec l1;
		for (int i = 0; i < l0.size(); i++){
			l1.push_back(l0[i]);
		}
		//Inserting higher hyperons for mu calculation
		for (int i = 0; i < n_ext.size(); i++){
			l1.push_back(n_ext[i]);
		}
		printf("Sm: %f %f \n", 2 * mu(0, l1, C) - mu(1, l1, C), C->M[6]);
		if (2 * mu(0, l1, C) - mu(1, l1, C) < C->M[6]){
			l0.push_back(0.0);
			return l0;
		}

		while ((status != 0) && xmax < 0.5*ntot){
			status = n_Sm_solve(ntot, n_ext, xmin, xmax, C, res);
			xmax += 0.1*ntot;
		}

		if (status != 0){
			vec in;
			in.push_back(0.0);
			for (int i = 0; i < n_ext.size(); i++){
				in.push_back(n_ext[i]);
			}
			res.clear();
			vec L0;
			L0 = n_L0(ntot, in, C);
			res.push_back(L0[0]);
			res.push_back(L0[1]);
			res.push_back(L0[2]);
			res.push_back(0.0);
		}
		//printf("%i \n", status);
		//printf("%f %f \n", 2 * mu(0, res, C) - mu(1, res, C), mu(3, res, C));
		/*printf("mu_n = %f, mu_p = %f, mu_L0 = %f, mu_Xm = %f, 2*mu_n - mu_p = %f \n",
			mu(0, res, C), mu(1, res, C), mu(2, res, C), mu(6, res, C), 2*mu(0, res, C) - mu(1, res, C));*/
		gsl_set_error_handler(def);
		return res;
	}

	const gsl_root_fsolver_type *T_xm_eq = gsl_root_fsolver_brent;
	gsl_root_fsolver *s_xm_eq = gsl_root_fsolver_alloc(T_xm_eq);

	double func_Xm(double x, void * params){
		func_l0_params * p = (func_l0_params *)params;
		set_const * C = p->C;
		double ntot = p->ntot;
		vec n_in;
// 		for (int i = 0; i < p->n_ext.size(); i++){
// 			n_in.push_back(p->n_ext[i]);
// 		}
		n_in.push_back(0.0);
		n_in.push_back(0.0);
		n_in.push_back(x);
		//��������� �������� ��� �������� ��������� �����

		//vec res_L0 = n_Sm(ntot, n_in, C);
		vec res_L0 = n_L0(ntot, n_in, C);
		double nn = res_L0[0];
		double np = res_L0[1];
		//Chemical potentials

		vec n;
		for (int i = 0; i < res_L0.size(); i++){
			n.push_back(res_L0[i]);
		}
		n.push_back(0.0);
		n.push_back(0.0);
		n.push_back(x);
	
		double f = f_eq(n, C);

		double dn = 1e-7;
		double E_p, E_m;
		//mu_n
		double mu_n = mu(0, n, C);
		//mu_p
		double mu_p = mu(1, n, C);

		//mu_p
		double mu_Xm = mu(6, n, C);

		printf("func_Xm: I'm working, don't kill me! x = %f, value = %f\n", x, mu_Xm - 2*mu_n + mu_p);

		return mu_Xm - 2*mu_n + mu_p;
	}

	int n_Xm_solve(double ntot, vec n_ext, double xmin, double xmax, set_const * C, vec & out){
		int status;
		int iter = 0, max_iter = 300;
		double x = 0.0, x_expect = (xmax - xmin) / 2.0;
		//double xmin = 0.0, xmax = 0.3 * ntot;
		gsl_function F;
		F.function = &func_Xm;
		func_l0_params par = { ntot, C, n_ext };
		F.params = &par;

		//test output
//  		std::ofstream ofs("func_sm.dat");
//  		double step = (xmax-xmin)/30;
//  		for (int i = 0; i < 30; i++){
//  			ofs << step*i << "   " << func_Xm(step*i, &par) << std::endl;
// 			std::cout << i << std::endl;
//  		}
//  		ofs.close();
		//end test


		status = gsl_root_fsolver_set(s_xm_eq, &F, xmin, xmax);
		if (status != 0){
			return status;
		}
		do {
			iter++;
			status = gsl_root_fsolver_iterate(s_xm_eq);
			x = gsl_root_fsolver_root(s_xm_eq);
			xmin = gsl_root_fsolver_x_lower(s_xm_eq);
			xmax = gsl_root_fsolver_x_upper(s_xm_eq);
			status = gsl_root_test_interval(xmin, xmax, 1e-6, 1e-6);
		} while (status == GSL_CONTINUE && iter < max_iter);

		if (status != 0) {
			std::cout << "error, status = " << status << std::endl;
		}
		//vec out;
		vec in;
		in.push_back(0.0);
		in.push_back(0.0);
		in.push_back(x); 
		vec L0;
		L0 = n_L0(ntot, in, C);
		out.push_back(L0[0]);
		out.push_back(L0[1]);
		out.push_back(L0[2]);
		//out.push_back(L0[3]);
		out.push_back(0.0);
		out.push_back(0.0);
		out.push_back(0.0);
		out.push_back(x);
		return status;
	}
	
	vec n_Xm(double ntot, vec n_ext, set_const * C){
		gsl_error_handler_t * def = gsl_set_error_handler_off();
		double xmax = 0.1*ntot;
		double xmin = 0.0;
		vec res;
		int status = 1;
		vec in_temp;
		in_temp.push_back(0.0);
		in_temp.push_back(0.0);
		in_temp.push_back(0.0);
		vec l0 = n_L0(ntot, in_temp, C);//TODO higher hyperons insert
		printf("Xm: %f %f \n", 2 * mu(0, l0, C) - mu(1, l0, C), C->M[6]);
		if (2 * mu(0, l0, C) - mu(1, l0, C) < C->M[6]){
			l0.push_back(0.0);
			l0.push_back(0.0);
			l0.push_back(0.0);
			l0.push_back(0.0);
			return l0;
		}

		while ((status != 0) && xmax < 0.5*ntot){
			status = n_Xm_solve(ntot, n_ext, xmin, xmax, C, res);
			xmax += 0.1*ntot;
		}

		if (status != 0){
			vec in;
			in.push_back(0.0);
			in.push_back(0.0);
			in.push_back(0.0);
			in.push_back(0.0);
			res.clear();
			vec L0;
			L0 = n_Sm(ntot, in, C);
			res.push_back(L0[0]);
			res.push_back(L0[1]);
			res.push_back(L0[2]);
			//res.push_back(L0[3]);
			res.push_back(0.0);
			res.push_back(0.0);
			res.push_back(0.0);
			res.push_back(0.0);
		}
		//printf("%i \n", status);
		printf("%f %f \n", 2 * mu(0, res, C) - mu(1, res, C), mu(6, res, C));
		/*printf("mu_n = %f, mu_p = %f, mu_L0 = %f, mu_Xm = %f, 2*mu_n - mu_p = %f \n",
			mu(0, res, C), mu(1, res, C), mu(2, res, C), mu(6, res, C), 2*mu(0, res, C) - mu(1, res, C));*/
		gsl_set_error_handler(def);
		return res;
	}
}
