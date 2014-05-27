#include "EoS_tab.h"
#include "set_const.h"
#include <fstream>
#include <iostream>
#include "EoS.h"
#include <math.h>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include "constants.h"

using namespace std;
EoS_tab::EoS_tab(set_const * C, char* filename,double start, double range, double step, bool forsed):table(boost::extents[ceil(range/step)][4])
{
	this->filename = filename;
	this->C = C;
	this->range = range;
	this->start = start;
	this->step = step;
	this->num = ceil(range/step);
	if (!isInitialized()){
		Init();
	}
	else{
		if (forsed){
			Init();
		}
		else{
			cout << this->C->name << "Initialized with filename " << this->filename << endl;
		}
	}

//	Init();
}


EoS_tab::~EoS_tab(void)
{
}


double EoS_tab::E(double n)
{
	int min = 0;
	double ncurrent = 0.0;
	double nprev = 10.0;
	for (int i = 0; i < num-1;i++){
		//cout << abs(this->table[i][0] - n) <<"  "<< nprev<<endl;
		if (abs(this->table[i][0] - n) < nprev)
		{
			min = i;
			nprev = abs(this->table[i][0] - n);
		}
	}
	return this->table[min][1];
}


double EoS_tab::P(double n)
{
	int min = 0;
	double ncurrent = 0.0;
	double nprev = 10.0;
	for (int i = 0; i < num-1;i++){
		//cout << abs(this->table[i][0] - n) <<"  "<< nprev<<endl;
		if (abs(this->table[i][0] - n) < nprev)
		{
			min = i;
			nprev = abs(this->table[i][0] - n);
		}
	}
	return this->table[min][2];
}


double EoS_tab::PofE(double E)
{
	int min = 0;
	double Ecurrent = 0.0;
	double Eprev = 1e20;
	for (int i = 0; i < num-1;i++){
		//cout << abs(this->table[i][0] - n) <<"  "<< nprev<<endl;
		if (abs(this->table[i][1] - E) < Eprev)
		{
			min = i;
			Eprev = abs(this->table[i][1] - E);
		}
	}
	return this->table[min][2];
}
double EoS_tab::NofE(double E){
	int min = 0;
	double Ecurrent = 0.0;
	double Eprev = 1e20;
	for (int i = 0; i < num-1;i++){
		//cout << abs(this->table[i][0] - n) <<"  "<< nprev<<endl;
		if (abs(this->table[i][1] - E) < Eprev)
		{
			min = i;
			Eprev = abs(this->table[i][1] - E);
		}
	}
	return this->table[min][0];
}

double EoS_tab::EofP(double P)
{
	int min = 0;
	double Pcurrent = 0.0;
	double Pprev = 1e20;
	for (int i = 0; i < num-1;i++){
		//cout << abs(this->table[i][0] - n) <<"  "<< nprev<<endl;
		if (abs(this->table[i][2] - P) < Pprev)
		{
			min = i;
			Pprev = abs(this->table[i][2] - P);
		}
	}
	return this->table[min][1];
}

bool EoS_tab::isInitialized(){
	vector<string> lines;
	string line;
	if (access(this->filename,0) != -1){
		ifstream ifs(filename);
		getline(ifs, line);
		cout << "Reading EoS " << this->C->name  << " from file " 
			<<this->filename << "..." << endl;
		cout << "NAME : " << line << endl;
		getline(ifs, line);
		cout << "CONSTANTS : " << line << endl;
		vector<string> SplitVector;
		boost::split(SplitVector, line, boost::is_any_of(","));
		double Cs = boost::lexical_cast<double>(SplitVector[0]);
		double Co = boost::lexical_cast<double> (SplitVector[1]);
		double Cr = boost::lexical_cast<double>(SplitVector[2]);
		double b = boost::lexical_cast<double>(SplitVector[3]);
		double c = boost::lexical_cast<double>(SplitVector[4]);
		if ((abs(Cs - C->C_s) / abs(Cs) > 1e-4) ||
			(abs(Co - C->C_o) / abs(Co) > 1e-4) ||
			(abs(Cr - C->C_r) / abs(Cr) > 1e-4) ||
			(abs(b - C->b) / abs(b) > 1e-4) ||
			(abs(c - C->c) / abs(c) > 1e-4))
		{
			cout << "Constants mismatch; reinitializing" << endl;
			return false;
		}
		else{
			cout << "Constants OK" << endl;
		}
		
		for (int i = 0; i < this->num - 2; i++){
			if (!ifs.eof()){
				getline(ifs, line);
				lines.push_back(line);
			}
			else{
				cout << "Initialization not finished; reinitializing" << endl;
				return false;
			}
		}

		
		vector<string>::iterator i_lines = lines.begin();
		//while(i_lines != lines.end()){
		//	cout << *i_lines++ << endl;
		//}

		i_lines = lines.begin();
		vector<string>::iterator i_split;

		typedef Table::index_range range;
		for (int i = 0; i < this->num - 2; i++){
			//cout << i << endl;
			line = *i_lines;
			boost::split(SplitVector, line, boost::is_any_of(","));	
			i_split = SplitVector.begin();
			while (i_split != SplitVector.end()){
				//cout << "- "<<*i_split<<endl;
				++i_split;
			}
			i_split = SplitVector.begin();
			Table::array_view<1>::type view = 
				this->table[boost::indices[i][range()]];
			Table::array_view<1>::type::iterator i_view = view.begin();
			while ((i_view != view.end())&&(i_split != SplitVector.end())){
				//cout << "--" << *i_split << endl;
				*i_view = boost::lexical_cast<double>(*i_split);
				++i_view;
				++i_split;
			}
			++i_lines;
		}
		cout << "DONE" << endl;
		return true;
	}
	return false;
}

int EoS_tab::Init(){
	cout << "Initializing" << endl;
	int num = ceil(this->range/this->step);

	ofstream ofs(this->filename);
	ofs << this->C->name << endl;
	ofs << this->C->C_s << "," << this->C->C_o<< "," << this->C->C_r
		<< "," << this->C->b<< "," << this->C->c << endl;
	vec n_ext;
	double np_eq;
	double n;
//#pragma omp parallel for private(n)
	for (int i = ceil(start/range); i< num - 1; i++){
		n = step*(i + 1);
		this->table[i][0] = n;
		np_eq = 0;//EoS::np_eq(n, n_ext,C);
		//this->table[i][1] = EoS::t_E(n-np_eq,np_eq, C);
		this->table[i][1] = EoS::E(n, C);
		this->table[i][2] = EoS::P(n, C);
		this->table[i][3] = np_eq;
		ofs << n << "," << this->table[i][1] << "," << this->table[i][2] << "," << this->table[i][3] << ","
			<< np_eq / n << endl;
//		cout << n << "," << this->table[i][1] << "," << this->table[i][2] << "," << this->table[i][3] << ","
//			<< np_eq/n << endl;
	}
	/*for (int i = ceil(start / range); i < num - 1; i++){
		n = step*(i + 1);
		ofs << n << "," << this->table[i][1] << "," << this->table[i][2] << "," << this->table[i][3] << endl;
		cout << n << "," << this->table[i][1] << "," << this->table[i][2] << "," << this->table[i][3] << endl;
	}*/
	cout << "Done" << endl;
	ofs.close();
	return 0;
}


