#ifndef _MIX_FIT_H

#define _MIX_FIT_H
#include "APR_fit.h"

namespace APR_fit{

	struct data_mix{
		size_t n_snm;
		size_t n_pnm;
		double * y;
		double * t;
		set_const* Init;
	};
	double fit_mix(set_const* Init, int nn);

}
#endif