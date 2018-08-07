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
	IntensityILS(int step_intensity, int bsize, int vsize, int maxdelta, int maxratio, double alpha, double beta) :
		ILS(bsize, vsize), step_intensity(step_intensity), maxdelta(maxdelta), maxratio(maxratio), alpha(alpha), beta(beta) { };
	virtual ~IntensityILS() { }

	virtual double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P);

	virtual bool acceptanceCriterion(double new_eval, double prev_eval){
		return false;
	}

	virtual void undoLast(Plan& p){
		p.undoLast2();
	}

	private:

	int step_intensity;

	double maxdelta;
	double maxratio;
	double alpha;
	double beta;



};

}
#endif /* SRC_ILS_INTENSITYILS_H_ */
