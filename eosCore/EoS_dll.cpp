// EoS_dll.cpp: ���������� ���������������� ������� ��� ���������� DLL.
//

#include "set_const.h"
#include <boost/python.hpp>
#include <string>
#include "EoS.h"
#include "MixedFit.h"
#include "APR_fit.h"
#include "constants.h"
#include "star.h"
#include "Auxillary.h"
#include "solving.h"
#include "hyper.h"
#include <math.h>
#include <boost/math/special_functions/math_fwd.hpp>

using namespace boost::python;

class Walecka: public set_const{
public: Walecka(double Cs, double Co, double Cr,
			double b, double c, double z) : 
		set_const(Cs,Co,Cr,b,c,z){};

		Walecka(std::string name, double Cs, double Co, double Cr,
			double b, double c, double z): 
		set_const(name, Cs, Co, Cr, b, c, z){}

		double eta_o(double f){
			return 1.0;
		}

		double f0;

		void set_f0(double f){
			this->f0 = f;
		}

		double eta_r(double f){
			return 1.0;
		}
		double eta_s(double f){
			return 1.0;
		}

		double phi_n(double f){
			return 1.0 - f;
		}
		double U(double f){
			return pow(m_n,4.0)*(this->b * pow(f,3.0)/3.0 + this->c*pow(f,4.0)/4.0);
		}
};

class gc_test: public set_const{
public: gc_test(double Cs, double Co, double Cr,
			double b, double c, double z) : 
		set_const(Cs,Co,Cr,b,c,z){
			
		}
		gc_test(){
			this->name = "Test Class";
			//Default constants:
			this->C_s = sqrt(164.5);
			this->C_o = sqrt(54.6);
			this->C_r = sqrt(121.7);
			this->b = 0.02028;
			this->c = 0.04716;
			this->z = 0;
			this->SetHyperConstants();
		}

		gc_test(std::string name, double Cs, double Co, double Cr,
			double b, double c, double z): 
		set_const(name, Cs, Co, Cr, b, c, z){

		}

		double eta_o(double f){
			return 1.0;
		}

		double f0;

		double gp_alpha;

		void set_gp_alpha(double a){
			this->gp_alpha = a;
		}

		double gp_beta;

		void set_gp_beta(double a){
			this->gp_beta = a;
		}

		void set_f0(double f){
			this->f0 = f;
		}

		double eta_r(double f){
			return 1.0;
		}
		double eta_s(double f){
			return 1.0;
		}

		double phi_n(double f){
			return 1.0 - f;
		}
		double U(double f){
			return pow(m_n,4.0)*(this->b * pow(f,3.0)/3.0 + this->c*pow(f,4.0)/4.0);
		}
};

class KVOR: public set_const{
public: KVOR(double Cs, double Co, double Cr,
			double b, double c, double z) : 
		set_const(Cs,Co,Cr,b,c,z){
			this->f0 = 0.195;
		};

		KVOR(std::string name, double Cs, double Co, double Cr,
			double b, double c, double z): 
		set_const(name, Cs, Co, Cr, b, c, z){
			this->f0 = 0.195;
		}

		void set_f0(double f){
			this->f0 = f;
		}

		double eta_o(double f){
			return (1.0 + this->z*this->f0)/(1.0 + this->z*f);
		}
		double eta_r(double f){
			return this->eta_o(f)/(this->eta_o(f) + 
				4*pow(this->C_o/this->C_r,2.0)*(this->eta_o(f)-1.0));
		}
		double eta_s(double f){
			return 1.0;
		}

		double phi_n(double f){
			return (1.0 - f);//*(1 + pow(f-f0, 2)*exp(f-f0));
		}
		double U(double f){
			return pow(m_n,4.0)*(this->b * pow(f,3.0)/3.0 + this->c*pow(f,4.0)/4.0);
		}
		double dU(double f){
			return pow(m_n, 4.0)*(this->b * pow(f, 2.0) + this->c*pow(f, 3.0));
		}
		double f0;
};

class gc_M : public KVOR{
public:
	gc_M(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "KVOR MOD";
		this->gp_alpha = 1.0;
		this->gp_beta = 1.0;
	}
	double gp_beta;
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}

	void set_gp_beta(double x){
		this->gp_beta = x;
	}
	double eta_o(double f){
		return this->gp_alpha/(1 + f*f);
	}
};

class gc_M4 : public KVOR{
public:
	gc_M4(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = a/(1 + bf^4)";
		this->gp_alpha = 1.0;
		this->gp_beta = 1.0;
	}
	double gp_beta;
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}

	void set_gp_beta(double x){
		this->gp_beta = x;
	}
	double eta_o(double f){
		return this->gp_alpha/(1 + this->gp_beta*f*f*f*f);
	}
};

class gc_M5 : public KVOR{
public:
	gc_M5(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = (1+z*f0/1+zf)^alpha";
		this->gp_alpha = 1.0;
		this->gp_beta = 1.0;
	}
	double gp_beta;
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}

	void set_gp_beta(double x){
		this->gp_beta = x;
	}
	double eta_o(double f){
		return pow( (1.0 + z*f0)/(1.0 + z*f), gp_alpha);
	}
};

class gc_M6 : public KVOR{
public:
	gc_M6(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = (1+z*f0/1+zf)^alpha, eta_r = 1";
		this->gp_alpha = 1.0;
		this->gp_beta = 1.0;
	}
	double gp_beta;
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}

	void set_gp_beta(double x){
		this->gp_beta = x;
	}
	double eta_o(double f){
		return pow( (1.0 + z*f0)/(1.0 + z*f), gp_alpha);
	}

	double eta_r(double f){
		if (f > 0.52){
			return 1.0/(1 + f);
		}
		else{
			return this->eta_o(f)/(this->eta_o(f) + 
				4*pow(this->C_o/this->C_r,2.0)*(this->eta_o(f)-1.0));
		}
	}

};

class gc_M2 : public KVOR{
public:
	gc_M2(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = 1/(1+f^3)";
	}

	double eta_o(double f){
		return 1.0/(1 + f*f*f);
	}
};


class gc_M7: public gc_M6{
public:
	gc_M7(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_M6(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR, eta_r = 1/(1+f)^b";
	}

	double eta_r(double f){
		return 1.0/pow(1 + f, gp_beta);
	}
};

class gc_M8: public KVOR{
public:
	gc_M8(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR, eta_r = KVOR^alpha";
		this->gp_alpha = 1.0;
	}
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}
	double eta_r(double f){
		return pow(abs(KVOR::eta_r(f)), this->gp_alpha);
	}

};

class gc_R1: public KVOR{
public:
	gc_R1(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR, eta_r = (eta_o/(2 - eta_o))^alpha";
		this->gp_alpha = 1.0;
	}
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}
	double eta_r(double f){
		return pow(this->eta_o(f)/(2 - this->eta_o(f)), this->gp_alpha);
	}

};

class gc_R2: public KVOR{
public:
	gc_R2(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR, eta_r = (eta_o/(2 + eta_o))^alpha";
		this->gp_alpha = 1.0;
	}
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}
	double eta_r(double f){
		return pow(this->eta_o(f)/(2 + this->eta_o(f)), this->gp_alpha);
	}

};

class gc_R3 : public gc_M5{
public:
	gc_R3(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_M5(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = (KVOR)^alpha, (eta_r = eta_o/(2-eta_o))^beta";
		this->gp_alpha = 1.0;
		this->gp_beta = 1.0;
	}
	
	double eta_o(double f){
		return pow( (1.0 + z*f0)/(1.0 + z*f), gp_alpha);
	}
	double eta_r(double f){
		return pow(this->eta_o(f)/(2 - this->eta_o(f)), gp_beta);
	}
};

class gc_OR1 : public KVOR{
public:
	gc_OR1(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	KVOR(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = (KVOR)^alpha, eta_r KVOR^beta";
		this->gp_alpha = 1.0;
		this->gp_beta = 1.0;
	}
	double gp_alpha, gp_beta;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}
	void set_gp_beta(double x){
		this->gp_beta = x;
	}

	double eta_o(double f){
		return pow(KVOR::eta_o(f), this->gp_alpha);
	}

	double eta_r(double f){
		return pow(abs(KVOR::eta_r(f)), this->gp_beta);
	}
};

class gc_OR2: public gc_OR1{
public:
	gc_OR2(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR1(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR^alpha, eta_r = (1+beta*f0)/(1+beta*f)";
	}

	double eta_r(double f){
		return (1+gp_beta*f0)/(1+gp_beta*f);
	}


};

class gc_OR3: public gc_OR2{
public:
	gc_OR3(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR2(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR^alpha, eta_r = ((1+beta*f0)/(1+beta*f))^2";
	}
	double gp_gamma;
	void set_gp_gamma(double x){
		this->gp_gamma = x;
	}
	double eta_r(double f){
		return pow(gc_OR2::eta_r(f), gp_gamma);
	}
};

class gc_OR4 : public gc_OR3{
public:
	gc_OR4(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR3(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR^alpha, eta_r = ((1+beta*f)/(1+beta*f0))^gamma";
	}

	double eta_r(double f){
		return pow((1 + gp_beta*f)/(1 + gp_beta*f0), gp_gamma);
	}
};

class gc_OR5 : public gc_OR4{
public:
	gc_OR5(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR4(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR(z)^alpha, eta_r = ((1+beta*f)/(1+beta*f0))^gamma";
		this->gp_z = 0.65;
	}

	double gp_z;
	void set_gp_z(double x){
		this->gp_z = x;
	}

	double eta_o(double f){
		return pow( (1.0 + gp_z*f0)/(1.0 + gp_z*f), gp_alpha);
	}
};

class gc_OR5_S : public gc_OR5{
public:
	gc_OR5_S(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR5(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR(z)^alpha + J, eta_r = ((1+beta*f)/(1+beta*f0))^gamma";
		this->gp_d = 0.0;
		this->gp_e = 0.0;
	}

	double gp_d, gp_e;

	void set_gp_d(double x){
		this->gp_d = x;
	}

	void set_gp_e(double x){
		this->gp_e = x;
	}

	double eta_s(double f){
		double c1 = this->c - 8*pow(this->C_s*this->b, 2.0) / 9;
		return pow(1.0 - 2*pow(C_s, 2.0)*b*f/3 - pow(C_s *f,2)*c1/2 + gp_d*pow(f,3)/3 + gp_e*pow(f,4)/4, -1);
	}

	double U(double f){
		return 0.0;
	}
};

class gc_OR5_S_J : public gc_OR5_S{
public:
	gc_OR5_S_J(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR5_S(Cs,Co,Cr,b,c,z){
		this->name = "O = KVOR(z)^alpha->a*th(w*(b-f)), R = ((1+beta*f)/(1+beta*f0))^gamma, S = ...";
	}

	double gp_a;
	double gp_b;
	double gp_fcut;
	double gp_w;

	void set_gp_a(double x){
		this->gp_a = x;
	}
	void set_gp_b(double x){
		this->gp_b = x;
	}
	void set_gp_fcut(double x){
		this->gp_fcut = x;
	}
	void set_gp_w(double x){
		this->gp_w = x;
	}
	double eta_o(double f){
		if (f < gp_fcut){
			return gc_OR5::eta_o(f);
		}
		else{
			//printf("val = %f \n", gp_a*(1 + tanh((f-gp_b)*(-10.0))));
			return gp_a*(1 + tanh((f-gp_b)*(-gp_w)));
		}
	}

};


class gc_OR5_S_Jcut : public gc_OR5_S_J{
public:
	gc_OR5_S_Jcut(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR5_S_J(Cs,Co,Cr,b,c,z){
		this->name = "O = KVOR(z)^alpha->a*th(w*(b-f)).. constant cut";
	}
	double gp_fup;
	void set_gp_fup(double x){
		this->gp_fup = x;
	}
	 
	double eta_o(double f){
		if (f > gp_fup){
			return gc_OR5_S_J::eta_o(gp_fup);
		}
		else{
			return gc_OR5_S_J::eta_o(f);
		}
	}
};

class gc_OR5_S_J2 : public gc_OR5_S_J{
public:
	gc_OR5_S_J2(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR5_S_J(Cs,Co,Cr,b,c,z){
		this->name = "O = KVOR(z)^alpha->a*th(w*(b-f)).. constant cut";
	}
	double gp_u;
	void set_gp_u(double x){
		this->gp_u = x;
	}

	double eta_o(double f){
		if (f < gp_fcut){
			return gc_OR5::eta_o(f);
		}
		else{
			//printf("val = %f \n", gp_a*(1 + tanh((f-gp_b)*(-10.0))));
			return gp_a*(gp_u + tanh((f-gp_b)*(-gp_w)));
		}
	}
};

class gc_OR5_J : public gc_OR5{
public:
	gc_OR5_J(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR5(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR(z)^alpha + J, eta_r = ((1+beta*f)/(1+beta*f0))^gamma";
		this->done = false;
	}

	double gp_fcut;
	void set_gp_fcut(double x){
		this->gp_fcut = x;
	}
	bool done;
	double eta_o(double f){
		double f1 = gp_fcut;
		double b = - gp_z*gp_alpha/(-2.0 - 2*f1*gp_z + f1*gp_z*gp_alpha);
		double a = gc_OR5::eta_o(f1)*pow(1.0 + b*f1,2.0);
		if (!done){
			printf("a = %f, b = %f \n", a, b);
			done = true;
		}
		if (f < f1){
			return gc_OR5::eta_o(f);
		}
		else{
			return a*pow(1.0/(1.0 + b*f),2.0);
		}
	}
};

class gc_OR5_J2 : public gc_OR5_J{
public:
	gc_OR5_J2(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR5_J(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR(z)^alpha + Jexp, eta_r = ((1+beta*f)/(1+beta*f0))^gamma";
		this->done = false;
	}

	double gp_a;
	double gp_b;
	void set_gp_a(double x){
		this->gp_a = x;
	}
	void set_gp_b(double x){
		this->gp_b = x;
	}

	double eta_o(double f){
		if (f < gp_fcut){
			return gc_OR5::eta_o(f);
		}
		else{
			//printf("val = %f \n", gp_a*(1 + tanh((f-gp_b)*(-10.0))));
			return gp_a*(1 + tanh((f-gp_b)*(-10.0)));
		}
	}
};

class gc_OR5_J3 : public gc_OR5_J2{
public:
	gc_OR5_J3(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_OR5_J2(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR(z)^alpha->a*th(w*(b-f)), eta_r = ((1+beta*f)/(1+beta*f0))^gamma";
		this->done = false;this->gp_w = 10.0;
	}
	double gp_w;
	void set_gp_w(double x){
		this->gp_w = x;
	}

	double eta_o(double f){
		if (f < gp_fcut){
			return gc_OR5::eta_o(f);
		}
		else{
			//printf("val = %f \n", gp_a*(1 + tanh((f-gp_b)*(-10.0))));
			return gp_a*(1 + tanh((f-gp_b)*(-gp_w)));
		}
	}

};


class gc_ORS1 : public gc_OR5_J2{
public:
	gc_ORS1(double Cs, double Co, double Cr,
		double b, double c, double z) :
	gc_OR5_J2(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR(z)^alpha + Jexp, eta_r=KVR, eta_s = (1+f0)/(1+f)";
		this->gp_a = 0.439652;
		this->gp_b = 0.684155;
	}

	double eta_s(double f){
		return (1+f0)/(1+f);
	}
};

class gc_ORS2 : public gc_ORS1{
public:
	gc_ORS2(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_ORS1(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR(z)^alpha + Jexp, eta_r=KVR, eta_s = 1 - f";
	}

	double eta_s(double f){
		return 1 - f;
	}

};


class gc_R4: public gc_R3{
public:
	gc_R4(double Cs, double Co, double Cr,
		double b, double c, double z) : 
	gc_R3(Cs,Co,Cr,b,c,z){
		this->name = "eta_o = KVOR, eta_r = (eta_o^2/(2 - eta_o))^alpha";
	}
	double eta_r(double f){
		return pow(eta_o(f), -gp_alpha);
	}
};

class gc_S1: public KVOR{
public:
	gc_S1(double cs, double co, double cr, 
		double b, double c, double z): KVOR(cs, co, cr ,b ,c, z){
			this->name = "KVOR eta_s = 1/(1+f)";
	}
	double eta_s(double f){
		return 1.0/(1+f);
	}
};

class gc_S2: public KVOR{
public:
	gc_S2(double cs, double co, double cr, 
		double b, double c, double z): KVOR(cs, co, cr ,b ,c, z){
			this->name = "KVOR eta_s = 1/(1+f)^a";
	}
	double gp_alpha;
	void set_gp_alpha(double x){
		this->gp_alpha = x;
	}


	double eta_s(double f){
		if (f > 0.52){
			return pow(1.52/(1+f),gp_alpha);
		}
		else{
			return 1.0;
		}
	}
};

class gc_J1: public KVOR{
public:
	gc_J1(double cs, double co, double cr, 
		double b, double c, double z): KVOR(cs, co, cr ,b ,c, z){
			this->name = "KVOR with omega jump at f*";
	}

	double eta_o(double f){
		if (f < 0.52){
			return KVOR::eta_o(f);
		}
		else{
			return 0.906478/(1 + pow(f - 0.260778, 2.0));
		}
	}
};

class gc_J2: public KVOR{
public:
	gc_J2(double cs, double co, double cr, 
		double b, double c, double z): KVOR(cs, co, cr ,b ,c, z){
			this->name = "KVOR with omega jump sqrt at f*";
	}

	double eta_o(double f){
		if (f < 0.52){
			return KVOR::eta_o(f);
		}
		else{
			return 1.3854/(1 + pow(f - 0.1218, 0.5));
		}
	}
};

class gc_J3: public gc_M6{
public:
	gc_J3(double cs, double co, double cr, 
		double b, double c, double z): gc_M6(cs, co, cr ,b ,c, z){
			this->name = "KVOR with omega jump alpha-dependent";
	}

	double eta_o(double f){

		if (f < 0.52){
			return KVOR::eta_o(f);
		}
		else{
			return (1.3854/(1 + pow(f - 0.1218, 0.5)))/(1 + pow(f - 0.52, gp_beta));
		}
	}
};

class gc_J4: public KVOR{
public:
	gc_J4(double cs, double co, double cr, 
		double b, double c, double z): KVOR(cs, co, cr ,b ,c, z){
			this->name = "KVOR with omega jump d-dependent";
			this->gp_d = 1.0; 
	}

	double gp_d;
	void set_gp_d(double x){
		this->gp_d = x;
	}

	double eta_o(double f){
		double d = this->gp_d;
		double a = 0.424701 +1.00521*sqrt(0.178506 +0.734897*d)-4.44089e-16*d;
		double b = (1.18343 *(-0.4225+sqrt(0.178506 +0.734897*d)-0.4303*d))/d;
		if (f < 0.52){
			return KVOR::eta_o(f);
		}
		else{
			return (a/(1 + pow(d*(f - b), 0.5)));
		}
	}
};

class Mod1 : public set_const{
public: Mod1(double Cs, double Co, double Cr,
		double b, double c) :
		set_const(Cs, Co, Cr, b, c, 0.0){};

		double eta_o(double f){
			return 1.0;
		}
		double eta_r(double f){
			return 1.0;
		}
		double eta_s(double f){
			return 1.0;
		}

		double phi_n(double f){
			//return 1.0 - f - f*f;
			return 1.0 / (1 + f);
		}
		double U(double f){
			return pow(m_n, 4.0)*(this->b * pow(f, 3.0) / 3.0 + this->c*pow(f, 4.0) / 4.0);
		}
};

class KVOR_enhanced : public set_const{
public: KVOR_enhanced(double Cs, double Co, double Cr,
		double b, double c, double z) :
		set_const(Cs, Co, Cr, b, c, z){
			this->f0 = 0.195;
			this->alpha = 1;
			this->beta = 1;
			this->gamma = 1;
};
		void set_z(double z){
			this->z = z;
		}
		void set_f0(double f){
			this->f0 = f;
		}

		double eta_o(double f){
			return pow((1.0 + this->z*this->f0) / (1.0 + this->z*f), 0.5);
		}
		double eta_r(double f){
			return this->eta_o(f) / (this->eta_o(f) +
				4 * pow(this->C_o / this->C_r, 2.0)*(this->eta_o(f) - 1.0));
		}
		double eta_s(double f){
			return 1.0;
		}

		double phi_n(double f){
			return (1.0 - f);///(1 + gamma*pow(f - f0, 3)*exp(alpha*(f - f0))*exp(beta*(f - 0.516775090581)));
		}
		double U(double f){
			return pow(m_n, 4.0)*(this->b * pow(f, 3.0) / 3.0 + this->c*pow(f, 4.0) / 4.0);
		}
		double dU(double f){
			return pow(m_n, 4.0)*(this->b * pow(f, 2.0) + this->c*pow(f, 3.0));
		}
		double f0;
		double alpha;
		double beta;
		double gamma;
};

class KVOR_enhanced2 : public set_const{
public: KVOR_enhanced2(double Cs, double Co, double Cr,
		double b, double c, double z) :
		set_const(Cs, Co, Cr, b, c, z){
			this->f0 = 0.195;
			this->alpha = 1;
			this->beta = 1;
			this->gamma = 1;
			this->gp_f_tr = 0.52;
			this->gp_width = 100.0;
			this->gp_power = 1;
}
		KVOR_enhanced2(std::string name, double Cs, double Co, double Cr,
			double b, double c, double z) :
		set_const(name, Cs, Co, Cr, b, c, z){
			this->f0 = 0.195;
			this->alpha = 1;
			this->beta = 1;
			this->gamma = 1;
			this->gp_f_tr = 0.52;
			this->gp_width = 100.0;
			this->gp_power = 1;
		}
		void set_thresh(double f){
			this->gp_f_tr = f;
		}

		void set_z(double z){
			this->z = z;
		}
		void set_f0(double f){
			this->f0 = f;
		}

		void set_power(double power){
			this->gp_power = power;
		}

		void set_width(double w){
			this->gp_width = w;
		}

		double eta_o(double f){
			return pow((1.0 + this->z*this->f0) / (1.0 + this->z*f), 1);
		}
		double eta_r(double f){
			return pow(this->eta_o(f) / (this->eta_o(f) +
				4 * pow(this->C_o / this->C_r, 2.0)*(this->eta_o(f) - 1.0)), 1.0) -
				(1 + tanh(gp_width*(f - gp_f_tr)))/pow(1.0 + f, gp_power);
		}
		double eta_s(double f){
			return 1.0;
		}

		double phi_n(double f){
			return (1.0 - f);///(1 + gamma*pow(f - f0, 3)*exp(alpha*(f - f0))*exp(beta*(f - 0.516775090581)));
		}
		double U(double f){
			return pow(m_n, 4.0)*(this->b * pow(f, 3.0) / 3.0 + this->c*pow(f, 4.0) / 4.0);
		}
		double dU(double f){
			return pow(m_n, 4.0)*(this->b * pow(f, 2.0) + this->c*pow(f, 3.0));
		}
		double f0;
		double alpha;
		double beta;
		double gamma;
		double gp_f_tr;
		double gp_width;
		double gp_power;
};

double (*fff)(double, double, set_const*) = &EoS::t_E;

double(*ff)(double, double, double, set_const*) = &EoS::t_E;

//double(*wr_feq)(double, double, set_const*) = &EoS::f_eq;
double wr_feq(double nn, double np, set_const * C){
	vec n;
	n.push_back(nn);
	n.push_back(np);
	return EoS::f_eq(n, C);
}


boost::python::list star1(double n, EoS_tab CT){
	double y[3];
	star(n, y, CT);
	boost::python::list l;
	for (int i = 0; i < 3; i++){
		l.append(y[i]);
	}
	return l;
}


boost::python::list star2(double n, set_const * C){
	double y[3];
	star(n, y, C);
	boost::python::list l;
	for (int i = 0; i < 3; i++){
		l.append(y[i]);
	}
	return l;
}

double f2(double f, double nn, double np, set_const * C){
	EoS::func_f_eq_params params = { nn, np, C };
	return EoS::func_f_eq(f, &params);
}

double np_test(double ntot, set_const * C){
	vec n;
	return EoS::np_eq(ntot, n, C);
}

boost::python::list wr_l0(double ntot, set_const * C){
	vec _in;
	_in.push_back(0.0);
	_in.push_back(0.0);
	_in.push_back(0.0);
	_in.push_back(0.0*ntot);
	vec res = Hyper_EoS::n_L0(ntot, _in, C);
	boost::python::list l;
	for (int i = 0; i < res.size(); i++){
		l.append(res[i]);
	}
	return l;
}

boost::python::list wr_xm(double ntot, set_const * C){
	vec _in;
	vec res = Hyper_EoS::n_Xm(ntot, _in, C);
	boost::python::list l;
	for (int i = 0; i < res.size(); i++){
		l.append(res[i]);
	}
	return l;
}

boost::python::list wr_sm(double ntot, set_const * C){
	vec _in;
	vec res = Hyper_EoS::n_Sm(ntot, _in, C);
	boost::python::list l;
	for (int i = 0; i < res.size(); i++){
		l.append(res[i]);
	}
	return l;
}

double wrap_feq(double x, double nn, double np, set_const * C){
	vec n;
	n.push_back(nn);
	n.push_back(np);
	EoS::func_f_eq_multi_params p = {n, C};
	return EoS::func_f_eq_multi(x, &p);
}

class gc_OR5_F : public gc_OR5{
public: gc_OR5_F(double Cs, double Co, double Cr,
			double b, double c, double z) :
		gc_OR5(Cs, Co, Cr, b, c, z){
			this->name = "gc_KVO(RR)_Phi";
			this->phi_c = 0.0;
			this->phi_cpow=5.0;
			this->phi_cfcut=0.45;
		};

		double gp_fcut;
		double gp_a;
		double gp_b;
		double gp_w;

		void set_gp_fcut(double x){
			this->gp_fcut = x;
		}
		
		void set_gp_a(double x){
			this->gp_a = x;
		}
		void set_gp_b(double x){
			this->gp_b = x;
		}

		void set_gp_w(double x){
			this->gp_w = x;
		}

		double phi_c;
		void set_phi_c(double x){
			this->phi_c = x;
		}

		double phi_cpow;
		void set_phi_cpow(double x){
					this->phi_cpow = x;
				}

		double phi_cfcut;
		void set_phi_cfcut(double x){
			this->phi_cfcut = x;
		}
		double phi_n(double f){
			double res = gc_OR5::phi_n(f);
			if (f > this->gp_fcut){
				res -= gp_a*pow(f-gp_fcut,3)/pow(cosh(gp_w*(f-gp_fcut)),2);
			}
			if (f > this->phi_cfcut){
				res += phi_c*pow(f - phi_cfcut, phi_cpow);
			}
			return res;
		}
};

class gc_OR5_S_F : public gc_OR5_F{
public: gc_OR5_S_F(double Cs, double Co, double Cr,
			double b, double c, double z) :
		gc_OR5_F(Cs, Co, Cr, b, c, z){
			this->name = "gc_KVO(RR)_Phi";
		};
	double gp_d;
	void set_gp_d(double x){
		this->gp_d = x;
	}

	double gp_e;
	void set_gp_e(double x){
		this->gp_e = x;
	}

	double eta_s(double f){
		double c1 = this->c - 8*pow(this->C_s*this->b, 2.0) / 9;
		return pow(1.0 - 2*pow(C_s, 2.0)*b*f/3 - pow(C_s *f,2)*c1/2 + gp_d*pow(f,3)/3 + gp_e*pow(f,4)/4, -1);
	}

	double U(double f){
		return 0.0;
	}

};

class gc_OR5_S_F_J : public gc_OR5_S_F{
public: gc_OR5_S_F_J(double Cs, double Co, double Cr,
			double b, double c, double z) :
		gc_OR5_S_F(Cs, Co, Cr, b, c, z){
			this->name = "gc_KVO(RR)_Phi";
			this->omega_c = 0.0;

		};

	double fcut_omega;
	void set_fcut_omega(double x){
		this->fcut_omega = x;
	}

	double omega_a;
	void set_omega_a(double x){
		this->omega_a = x;
	}

	double omega_b;

	double rho_a, rho_b;
	double fcut_rho;
	void set_fcut_rho(double x){
		this->fcut_rho = x;
	}
	void set_rho_a(double x){
		this->rho_a = x;
	}

	void set_rho_b(double x){
			this->rho_b = x;
	}

	double f_tilde;
	void set_f_tilde(double x){
		this->f_tilde = x;
	}


	void set_omega_b(double x){
			this->omega_b = x;
	}

	double omega_c;
	void set_omega_c(double x){
		this->omega_c = x;
	}

	double omega_a_l, fcut_omega_l;
	double sigma_a_l, fcut_sigma_l;
	double phi_a_l, fcut_phi_l;

	void set_phi_a_l(double x){
		this->phi_a_l = x;
	}

	void set_fcut_phi_l(double x){
		this->fcut_phi_l = x;
	}
	void set_omega_a_l(double x){
		this->omega_a_l = x;
	}

	void set_fcut_omega_l(double x){
			this->fcut_omega_l = x;
	}

	void set_sigma_a_l(double x){
		this->sigma_a_l = x;
	}

	void set_fcut_sigma_l(double x){
		this->fcut_sigma_l = x;
	}

	double eta_o(double f){
		if (f < fcut_omega_l){
			return gc_OR5::eta_o(f) + omega_a_l*pow(fcut_omega_l - f, 3.0);
		}
		if ((f < fcut_omega) && (f > fcut_omega_l)){
			return gc_OR5::eta_o(f);
		}
		else{
			return gc_OR5::eta_o(f) - omega_a*pow(f-fcut_omega,3)*cosh(omega_b*(f-fcut_omega)) +
					omega_c * pow(f - fcut_omega, 5.0);
		}
	}


	double eta_r(double f){
		if (f < fcut_rho){
			return gc_OR5::eta_r(f);
		}
		else{
			return gc_OR5::eta_r(f) - rho_a*pow(f-fcut_rho,3);
		}
	}

	double eta_s(double f){
		double c1 = this->c - 8*pow(this->C_s*this->b, 2.0) / 9;
		double res = pow(1.0 - 2*pow(C_s, 2.0)*b*f/3 - pow(C_s *f,2)*c1/2 + gp_d*pow(f,3)/3 + gp_e*pow(f/f_tilde,6)/6, -1);

		if (f < fcut_sigma_l){
			res += sigma_a_l*pow(f - fcut_sigma_l, 3.0);
		}
		return res;
	}
	double phi_n(double f){
		if (f > this->fcut_phi_l ){
			return gc_OR5_F::phi_n(f);
		}
		else{
			return gc_OR5_F::phi_n(f) - phi_a_l*f*pow((f - fcut_phi_l), 3);
		}
	}
};


class gc_OR5_Fexp: public gc_OR5_F{
public: gc_OR5_Fexp(double Cs, double Co, double Cr,
			double b, double c, double z) :
		gc_OR5_F(Cs, Co, Cr, b, c, z){
			this->name = "gc_KVO(RR)_Phi";
		};
};


BOOST_PYTHON_MODULE(eosCore){
	class_<set_const>("set_const", init<double, double, double, double, double, double>())
		.def(init<std::string, double, double, double, double, double, double>())
		.add_property("C_s", &set_const::C_s)
		.add_property("C_o", &set_const::C_o)
		.add_property("C_r", &set_const::C_r)
		.add_property("b", &set_const::b)
		.add_property("c", &set_const::c)
		.add_property("a", &set_const::a)
		.add_property("z", &set_const::z)
		.add_property("name", &set_const::name)
		.def("phi_n", &set_const::phi_n)
		.def("eta_s", &set_const::eta_s)
		.def("eta_o", &set_const::eta_o)
		.def("eta_r", &set_const::eta_r)
		.def("__repr__", &set_const::repr)
		.def("init", &set_const::init)
		.def("set_name", &set_const::set_name)
	;

	class_<KVOR, bases<set_const> >("KVOR", init<double, double, double, double, double, double>())
		.def(init<std::string, double, double, double, double, double, double>())
		.def("set_f0", &KVOR::set_f0)
		.def("dU", &KVOR::dU)
		.add_property("f0", &KVOR::f0)
		;
	class_<KVOR_enhanced, bases<set_const> >("KVOR_enhanced", init<double, double, double, double, double, double>())
		.def("set_f0", &KVOR_enhanced::set_f0)
		.def("dU", &KVOR_enhanced::dU)
		.def("set_z", &KVOR_enhanced::set_z)
		.add_property("f0", &KVOR_enhanced::f0)
		.add_property("alpha", &KVOR_enhanced::alpha)
		.add_property("beta", &KVOR_enhanced::beta)
		.add_property("gamma", &KVOR_enhanced::gamma)
		;
	class_<KVOR_enhanced2, bases<set_const> >("KVOR_enhanced2", init<double, double, double, double, double, double>())
		.def(init<std::string, double, double, double, double, double, double>())
		.def("set_f0", &KVOR_enhanced2::set_f0)
		.def("dU", &KVOR_enhanced2::dU)
		.def("set_z", &KVOR_enhanced2::set_z)
		.add_property("f0", &KVOR_enhanced2::f0)
		.add_property("alpha", &KVOR_enhanced2::alpha)
		.add_property("beta", &KVOR_enhanced2::beta)
		.add_property("gamma", &KVOR_enhanced2::gamma)
		.add_property("gp_f", &KVOR_enhanced2::gp_f_tr)
		.def("set_gp_f", &KVOR_enhanced2::set_thresh)
		.add_property("gp_width", &KVOR_enhanced2::gp_width)
		.def("set_gp_width", &KVOR_enhanced2::set_width)
		.add_property("gp_power", &KVOR_enhanced2::gp_power)
		.def("set_gp_power", &KVOR_enhanced2::set_power)
		;
	class_<Walecka, bases<set_const> >("Walecka", init<double, double, double, double, double, double>())
		.def(init<std::string, double, double, double, double, double, double>())
		.def("set_f0", &Walecka::set_f0)
		.add_property("f0", &Walecka::f0)
		;

	class_<gc_test, bases<set_const> >("gc_test", init<double, double, double, double, double, double>())
		.def(init<std::string, double, double, double, double, double, double>())
		.def(init<>())
		.def("set_f0", &gc_test::set_f0)
		.add_property("f0", &gc_test::f0)
		.add_property("gp_alpha", &gc_test::gp_alpha)
		.def("set_gp_alpha", &gc_test::set_gp_alpha)
		.add_property("gp_beta", &gc_test::gp_beta)
		.def("set_gp_beta", &gc_test::set_gp_beta)
		;
	class_<gc_M,bases<KVOR> >("gc_M", init<double, double, double, double, double, double>())
	//.def(init<>())
	.def("set_gp_alpha", &gc_M::set_gp_alpha)
	.add_property("gp_alpha", &gc_M::gp_alpha)
	.def("set_gp_beta", &gc_M::set_gp_beta)
	.add_property("gp_beta", &gc_M::gp_beta)
	;
	class_<gc_M2,bases<KVOR> >("gc_M2", init<double, double, double, double, double, double>())
		//.def(init<>())
		;
	class_<Mod1, bases<set_const> >("Mod1", init<double, double, double, double, double>());

	def("E", EoS::E);
	def("P", EoS::P);
	def("t_E", *fff);
	def("t_E2", *ff);
	def("np_eq", np_test);
	def("f_eq", *wr_feq);
	def("APR_fit", APR_fit::fit);
	def("APR_fit_n", APR_fit::fit_n);
	def("APR_fit_mix", APR_fit::fit_mix);
	def("p_f", EoS::p_f);
	def("mu_e", EoS::mu_e);
//	def("mu_e", EoS::mu_e);
// 	class_<EoS::func_np_params>("np_params")
// 		.def("set", &EoS::func_np_params::set)
// 		.add_property("n", &EoS::func_np_params::n)
// 		.add_property("C", &EoS::func_np_params::C)
// 	;

	class_<EoS_tab>("EoS_tab", init<set_const *,char*,double, double, double,bool>());
	def("star", star1);
	def("star_notab", star2);
	def("SaturationDensity", Auxillary::SaturationDensity);
	def("solve", solve);
	def("f2", f2);
	def("Asymm", Auxillary::Asymm);
	def("Compr", Auxillary::Compr);
	def("L", Auxillary::L);
	def("test_l0", wr_l0);
	def("test_xm", wr_xm);
	def("test_sm", wr_sm);
	def("func_feq", EoS::func_f_eq);
	def("wrap_feq", &wrap_feq);
	class_<gc_M4,bases<KVOR> >("gc_M4", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_M4::set_gp_alpha)
		.add_property("gp_alpha", &gc_M4::gp_alpha)
		.def("set_gp_beta", &gc_M4::set_gp_beta)
		.add_property("gp_beta", &gc_M4::gp_beta)
		;
	class_<gc_M5,bases<KVOR> >("gc_M5", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_M5::set_gp_alpha)
		.add_property("gp_alpha", &gc_M5::gp_alpha)
		.def("set_gp_beta", &gc_M5::set_gp_beta)
		.add_property("gp_beta", &gc_M5::gp_beta)
		;

	class_<gc_M6,bases<KVOR> >("gc_M6", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_M6::set_gp_alpha)
		.add_property("gp_alpha", &gc_M6::gp_alpha)
		.def("set_gp_beta", &gc_M6::set_gp_beta)
		.add_property("gp_beta", &gc_M6::gp_beta)
		;
	class_<gc_M7,bases<gc_M6> >("gc_M7", init<double, double, double, double, double, double>());
	class_<gc_S1,bases<KVOR> >("gc_S1", init<double, double, double, double, double, double>());
	class_<gc_S2,bases<KVOR> >("gc_S2", init<double, double, double, double, double, double>())
		.add_property("gp_alpha", &gc_S2::gp_alpha)
		.def("set_gp_alpha", &gc_S2::set_gp_alpha)
		;
	class_<gc_M8,bases<KVOR> >("gc_M8", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_M8::set_gp_alpha)
		.add_property("gp_alpha", &gc_M8::gp_alpha)
		;
	class_<gc_R1,bases<KVOR> >("gc_R1", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_R1::set_gp_alpha)
		.add_property("gp_alpha", &gc_R1::gp_alpha)
		;
	class_<gc_R2,bases<KVOR> >("gc_R2", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_R2::set_gp_alpha)
		.add_property("gp_alpha", &gc_R2::gp_alpha)
		;
	class_<gc_R3,bases<gc_M5> >("gc_R3", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_R3::set_gp_alpha)
		.add_property("gp_alpha", &gc_R3::gp_alpha)
		.add_property("bp_beta", &gc_R3::gp_beta)
		.def("set_gp_beta", &gc_R3::set_gp_beta)
		;
	class_<gc_R4,bases<gc_R3> >("gc_R4", init<double, double, double, double, double, double>())
	;
	class_<gc_J1, bases<KVOR> >("gc_J1", init<double, double, double, double, double, double>())
		;
	class_<gc_J2, bases<KVOR> >("gc_J2", init<double, double, double, double, double, double>())
		;
	class_<gc_J3, bases<gc_M6> >("gc_J3", init<double, double, double, double, double, double>())
		;
	class_<gc_J4, bases<KVOR> >("gc_J4", init<double, double, double, double, double, double>())
		.def("set_gp_d", &gc_J4::set_gp_d)
		.add_property("gp_d", &gc_J4::gp_d)
		;
	class_<gc_OR1,bases<KVOR> >("gc_OR1", init<double, double, double, double, double, double>())
		//.def(init<>())
		.def("set_gp_alpha", &gc_OR1::set_gp_alpha)
		.add_property("gp_alpha", &gc_OR1::gp_alpha)
		.add_property("gp_beta", &gc_OR1::gp_beta)
		.def("set_gp_beta", &gc_OR1::set_gp_beta)
		;
	class_<gc_OR2, bases<gc_OR1> >("gc_OR2", init<double, double, double, double, double, double>())
		;
	class_<gc_OR3, bases<gc_OR2> >("gc_OR3", init<double, double, double, double, double, double>())
		.add_property("gp_gamma", &gc_OR3::gp_gamma)
		.def("set_gp_gamma", &gc_OR3::set_gp_gamma)
		;
	class_<gc_OR4, bases<gc_OR3> >("gc_OR4", init<double, double, double, double, double, double>())
		;
	class_<gc_OR5, bases<gc_OR4> >("gc_OR5", init<double, double, double, double, double, double>())
		.add_property("gp_z", &gc_OR5::gp_z)
		.def("set_gp_z", &gc_OR5::set_gp_z)
		;
	class_<gc_OR5_J, bases<gc_OR5> >("gc_OR5_J", init<double, double, double, double, double, double>())
		.add_property("gp_fcut", &gc_OR5_J::gp_fcut)
		.def("set_gp_fcut", &gc_OR5_J::set_gp_fcut)
		;
	class_<gc_OR5_J2, bases<gc_OR5_J> >("gc_OR5_J2", init<double, double, double, double, double, double>())
		.add_property("gp_a", &gc_OR5_J2::gp_a)
		.def("set_gp_a", &gc_OR5_J2::set_gp_a)
		.add_property("gp_b", &gc_OR5_J2::gp_b)
		.def("set_gp_b", &gc_OR5_J2::set_gp_b)
		;
	class_<gc_ORS1, bases<gc_OR5_J2> >("gc_ORS1", init<double, double, double, double, double, double>())
		;
	class_<gc_ORS2, bases<gc_ORS1> >("gc_ORS2", init<double, double, double, double, double, double>())
		;
	class_<gc_OR5_J3, bases<gc_OR5_J2> >("gc_OR5_J3", init<double, double, double, double, double, double>())
		.add_property("gp_w", &gc_OR5_J3::gp_w)
		.def("set_gp_w", &gc_OR5_J3::set_gp_w)
		;
	class_<gc_OR5_S, bases<gc_OR5> >("gc_OR5_S", init<double, double, double, double, double, double>())
		.add_property("gp_d", &gc_OR5_S::gp_d)
		.def("set_gp_d", &gc_OR5_S::set_gp_d)
		.add_property("gp_e", &gc_OR5_S::gp_e)
		.def("set_gp_e", &gc_OR5_S::set_gp_e)
		;
	class_<gc_OR5_S_J, bases<gc_OR5_S> >("gc_OR5_S_J", init<double, double, double, double, double, double>())
		.add_property("gp_a", &gc_OR5_S_J::gp_a)
		.def("set_gp_a", &gc_OR5_S_J::set_gp_a)
		.add_property("gp_b", &gc_OR5_S_J::gp_b)
		.def("set_gp_b", &gc_OR5_S_J::set_gp_b)
		.add_property("gp_fcut", &gc_OR5_S_J::gp_fcut)
		.def("set_gp_fcut", &gc_OR5_S_J::set_gp_fcut)
		.add_property("gp_w", &gc_OR5_S_J::gp_w)
		.def("set_gp_w", &gc_OR5_S_J::set_gp_w)
		;
	class_<gc_OR5_S_Jcut, bases<gc_OR5_S_J> >("gc_OR5_S_Jcut", init<double, double, double, double, double, double>())
		.add_property("gp_fup", &gc_OR5_S_Jcut::gp_fup)
		.def("set_gp_fup", &gc_OR5_S_Jcut::set_gp_fup)
		;
	class_<gc_OR5_S_J2, bases<gc_OR5_S_J> >("gc_OR5_S_J2", init<double, double, double, double, double, double>())
		.add_property("gp_u", &gc_OR5_S_J2::gp_u)
		.def("set_gp_u", &gc_OR5_S_J2::set_gp_u)
		;
	class_<gc_OR5_F, bases<gc_OR5> >("gc_OR5_F", init<double, double, double, double, double, double>())
		.add_property("gp_a", &gc_OR5_F::gp_a)
		.def("set_gp_a", &gc_OR5_F::set_gp_a)
		.add_property("gp_b", &gc_OR5_F::gp_b)
		.def("set_gp_b", &gc_OR5_F::set_gp_b)
		.add_property("gp_fcut", &gc_OR5_F::gp_fcut)
		.def("set_gp_fcut", &gc_OR5_F::set_gp_fcut)
		.add_property("gp_w", &gc_OR5_F::gp_w)
		.def("set_gp_w", &gc_OR5_F::set_gp_w)
		.add_property("phi_c", &gc_OR5_F::phi_c)
		.def("set_phi_c", &gc_OR5_F::set_phi_c)
		.add_property("phi_cpow", &gc_OR5_F::phi_cpow)
		.def("set_phi_cpow", &gc_OR5_F::set_phi_cpow)
		.add_property("phi_cfcut", &gc_OR5_F::phi_cfcut)
		.def("set_phi_cfcut", &gc_OR5_F::set_phi_cfcut)
		;
	class_<gc_OR5_Fexp, bases<gc_OR5_F> >("gc_OR5_Fexp", init<double, double, double, double, double, double>());
	class_<gc_OR5_S_F, bases<gc_OR5_F > >("gc_OR5_S_F", init<double, double, double, double, double, double>())
			.add_property("gp_d", &gc_OR5_S_F::gp_d)
			.def("set_gp_d", &gc_OR5_S_F::set_gp_d)
			.add_property("gp_e", &gc_OR5_S_F::gp_e)
			.def("set_gp_e", &gc_OR5_S_F::set_gp_e)
			;
	class_<gc_OR5_S_F_J, bases<gc_OR5_S_F > >("gc_OR5_S_F_J", init<double, double, double, double, double, double>())
			.add_property("fcut_omega", &gc_OR5_S_F_J::fcut_omega)
			.def("set_fcut_omega", &gc_OR5_S_F_J::set_fcut_omega)
			.add_property("omega_a", &gc_OR5_S_F_J::omega_a)
			.def("set_omega_a", &gc_OR5_S_F_J::set_omega_a)
			.add_property("omega_b", &gc_OR5_S_F_J::omega_b)
			.def("set_omega_b", &gc_OR5_S_F_J::set_omega_b)
			.add_property("rho_a", &gc_OR5_S_F_J::rho_a)
			.def("set_rho_a", &gc_OR5_S_F_J::set_rho_a)
			.add_property("rho_b", &gc_OR5_S_F_J::rho_b)
			.def("set_rho_b", &gc_OR5_S_F_J::set_rho_b)
			.add_property("f_tilde", &gc_OR5_S_F_J::f_tilde)
			.def("set_f_tilde", &gc_OR5_S_F_J::set_f_tilde)

			.add_property("omega_c", &gc_OR5_S_F_J::omega_c)
			.def("set_omega_c", &gc_OR5_S_F_J::set_omega_c)
			.add_property("fcut_rho", &gc_OR5_S_F_J::fcut_rho)
			.def("set_fcut_rho", &gc_OR5_S_F_J::set_fcut_rho)
			.add_property("fcut_omega_l", &gc_OR5_S_F_J::fcut_omega_l)
			.def("set_fcut_omega_l", &gc_OR5_S_F_J::set_fcut_omega_l)
			.add_property("fcut_sigma_l", &gc_OR5_S_F_J::fcut_sigma_l)
			.def("set_fcut_sigma_l", &gc_OR5_S_F_J::set_fcut_sigma_l)
			.add_property("fcut_phi_l", &gc_OR5_S_F_J::fcut_phi_l)
			.def("set_fcut_phi_l", &gc_OR5_S_F_J::set_fcut_phi_l)

			.add_property("omega_a_l", &gc_OR5_S_F_J::omega_a_l)
			.def("set_omega_a_l", &gc_OR5_S_F_J::set_omega_a_l)
			.add_property("sigma_a_l", &gc_OR5_S_F_J::sigma_a_l)
			.def("set_sigma_a_l", &gc_OR5_S_F_J::set_sigma_a_l)
			.add_property("phi_a_l", &gc_OR5_S_F_J::phi_a_l)
			.def("set_phi_a_l", &gc_OR5_S_F_J::set_phi_a_l)
			;
}


