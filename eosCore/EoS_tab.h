#include "set_const.h"
#include <boost/multi_array.hpp>
#ifndef EOS_TAB_H_
#define EOS_TAB_H_


class EoS_tab
{
public:
	EoS_tab(set_const * C, char* filename, double start, double range, double step);
	~EoS_tab(void);
	double E(double n);
	double P(double n);
	double PofE(double E);
	double EofP(double P);
	double NofE(double E);
	char *filename;
	set_const *C;
	double step;
	double range;
	double start;
	int num;
	typedef boost::multi_array<double, 2> Table;
	Table table;
	int Init();



private: bool isInitialized();
		
};
#endif /* EOS_TAB_H_ */
