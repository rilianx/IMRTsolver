/*
 * IntensityILS.h
 *
 *  Created on: 2 ago. 2018
 *      Author: iaraya
 */

#include "ILS.h"

#ifndef ILS_INTENSITYILS_H_
#define ILS_INTENSITYILS_H_

namespace imrt {

class IntensityILS  : public ILS {
public:
	IntensityILS(int bsize, int vsize, bool maxdelta, bool maxratio, bool alpha, bool beta) :
		ILS(bsize, vsize), maxdelta(maxdelta), maxratio(maxratio), alpha(alpha), beta(beta) { };
	virtual ~IntensityILS() { }

	virtual double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P);

	virtual bool acceptanceCriterion(double new_eval, double prev_eval){
		return false;
	}

	private:

	bool maxdelta;
	bool maxratio;
	bool alpha;
	bool beta;

};

}
#endif /* SRC_ILS_INTENSITYILS_H_ */
