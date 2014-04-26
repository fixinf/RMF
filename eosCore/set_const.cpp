/*
 * set_const.cpp
 *
 *  Created on: 16.05.2013
 *      Author: fixinf
 */
#include "set_const.h"
#include <math.h>
#include "constants.h"
#include <iostream>
#include <string>

double set_const::phi_n(double f){
	return 1.0 - f;
	//return 1.0/(1.0 + f);
}

double set_const::diff_phi_n(double f){
	double df = 0.0001;
	return (this->phi_n(f+df) - this->phi_n(f))/(df);
}

double set_const::eta_r(double f){
	return 1.0;
	//return this->eta_o(f)/(this->eta_o(f) + 
		//4*pow(this->C_o/this->C_r,2.0)*(this->eta_o(f)-1.0));
}

double set_const::dU(double f){
	return pow(m_n,4.0)*(b*pow(f,2.0) + c*pow(f,3.0));
}


double set_const::eta_o(double f){
	return 1;
	//return (1 + this->z*0.195)/( 1 + this->z*f);
}


double set_const::eta_s(double f){
	return 1.0;
}

double set_const::U(double f){

	//std::cout << "Cs = " << C_s << "B = " << this->b << " C = " << this->c << std::endl;
	return pow(m_n,4.0)*(b * pow(f,3.0)/3.0 + c*pow(f,4.0)/4.0);
}

std::string set_const::repr(){
	return this->name;
}


set_const::set_const(double C_s, double C_o, double C_r, double b, double c, double z) {
	this->init(C_s,C_o,C_r,b,c,z);
}

set_const::set_const(std::string name, double C_s, double C_o , double C_r, double b, double c, double z){
	this->init(C_s,C_o,C_r,b,c,z);
	this->name = name;
}

void set_const::set_name(std::string name){
	this->name = name;
}

void set_const::init(double C_s, double C_o, double C_r, double b, double c, double z){
	this->Hyper = true;
	this->SetHyperConstants();
	this->C_o = C_o;
	this->C_r = C_r;
	this->C_s = C_s;
	this->b = b;
	this->c = c;
	this->z = z;
}

int set_const::SetHyperConstants(){
	//Порядок следования барионов:
	// n p L0 S- S0 S+ X- X0

	this->X_o.clear();
	this->X_s.clear();
	this->X_r.clear();
	this->T.clear();
	this->Q.clear();
	this->M.clear();
	if (this->Hyper){
		double xo[8] = { 1.0, 1.0, 2.0 / 3.0, 2.0/3.0, 2.0/3.0, 2.0/3.0, 1.0 / 3.0, 1.0 / 3.0 };
		double xr[8] = { 1.0, 1.0, 0, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
		//double xr[8] = { 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 1.0};
		double ebind[8] = {0, 0, -30, 30, 30, 30, -18, -18};
		double m[8] = { 938, 938, 1116, 1195, 1195, 1195 , 1317, 1317};
		//	double xs[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		//	double xo[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		//	double xr[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		double xs[8];
		xs[0] = 1.0;
		xs[1] = 1.0;
		for (int i = 2; i < 8; i++){
			xs[i] = (80.73*xo[i] - ebind[i]) / 140.70;
		}
		double t[8] = { 0.5, -0.5, -1.0, 0.0, 1.0, -0.5, 0.5 };
		double q[8] = { 0, 1, 0, -1, 0, 1, -1, 0 };
		for (int i = 0; i < 8; i++){
				this->X_s.push_back(xs[i]);
				this->X_o.push_back(xo[i]);
				this->X_r.push_back(xr[i]);
				this->T.push_back(t[i]);
				this->Q.push_back(q[i]);
				this->M.push_back(m[i]);
		}
	} 
	else{
		double xo[8] = { 1, 1, 0,0,0,0,0,0 };
		double xr[8] = { 1, 1, 0,0,0,0,0,0 };
		double ebind[8] = { 0, 0, 0,0,0,0,0,0 };
		double m[8] = { 938, 938, 1116, 1195, 1195, 1195, 1317, 1317 };
		//	double xs[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		//	double xo[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		//	double xr[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		double xs[8] = { 1, 1, 0, 0, 0, 0, 0, 0 };
		double t[8] = { 0.5, -0.5, -1.0, 0.0, 1.0, -0.5, 0.5 };
		double q[8] = { 0, 1, 0, -1, 0, 1, -1, 0 };
		for (int i = 0; i < 8; i++){
			if (this->Hyper){
				this->X_s.push_back(xs[i]);
				this->X_o.push_back(xo[i]);
				this->X_r.push_back(xr[i]);
				this->T.push_back(t[i]);
				this->Q.push_back(q[i]);
				this->M.push_back(m[i]);
			}
		}
	}
	
	return 0;
}


//set_const::set_const(double C_s, double C_o, double C_r, double b, double c, 
//					 double z,std::function<double(double)> f_eta_o) {
//	// TODO Auto-generated constructor stub
//	set_const::C_o = C_o;
//	set_const::C_r = C_r;
//	set_const::C_s = C_s;
//	set_const::b = b;
//	set_const::c = c;
//	set_const::z = z;
//	this->eta_o = f_eta_o;
//}

//set_const::set_const(std::string name, double C_s, double C_o, double C_r, 
//					 double b, double c, double z,
//					 scaling phi_n = NULL,scaling eta_s = NULL, scaling eta_o = NULL, scaling eta_r= NULL) {
//	// TODO Auto-generated constructor stub
//	set_const::C_o = C_o;
//	set_const::C_r = C_r;
//	set_const::C_s = C_s;
//	set_const::b = b;
//	set_const::c = c;
//	set_const::z = z;
//	set_const::name = name;
//	this->phi_n = phi_n;
//	this->eta_s = eta_s;
//	this->eta_o = eta_o;
//	this->eta_r = eta_r;
//
//	//std::cout << "BBBBBB" << this->b <<"CCCC" << this->c << std::endl;
//	// TODO использование переданной функции как метода или же нафиг
//}
