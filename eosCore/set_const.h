#include <string>
#include <functional>
#include <vector>
/*
 * set_const.h
 *
 *  Created on: 16.05.2013
 *      Author: fixinf
 */

#ifndef SET_CONST_H_
#define SET_CONST_H_

typedef std::vector<double> vec;



class set_const {
public:
	
	//			Cs		Co		Cr		b		c		z
	set_const(){

	}
	set_const(double, double, double, double, double, double);
	//			Cs		Co		Cr		b		c		z

	set_const(std::string name, double Cs, double Co , double Cr, double b, double c, double z);
	void init(double, double, double, double, double, double);
	double diff_phi_n(double);
	virtual double U(double);
	double dU(double);
	void set_name(std::string name);
	virtual double phi_n(double);
	virtual double eta_s(double);
	virtual double eta_o(double);
	virtual double eta_r(double);
	double C_s;
	double C_o;
	double C_r;
	std::string name;
	double b;
	double c;
	double z;
	double a;
	std::vector<double> X_s;
	std::vector<double> X_o;
	vec X_r;
	std::vector<double> Q;
	std::vector<double> T;
	vec M;
	bool Hyper;
	int SetHyperConstants();
	std::string repr();
	//std::function<double(double)> phi_n;
	//std::function<double(double)> eta_s;
	//std::function<double(double)> eta_o;
	//std::function<double(double)> eta_r;
};

#endif /* SET_CONST_H_ */
