#include "set_const.h"


#ifndef EOS_TABULATED_H_
#define EOS_TABULATED_H_

class EoS_tab{
public: double E(double);
		double P(double);
		double EofP(double);
		double PofE(double);
		char * filename;
		double range;
		double step;
		char * name;
		EoS_tab(set_const, char*, double, double);
private: bool isInitialized(set_const);
		 bool init(set_const);

};

#endif /* EOS_TABULATED_H_ */