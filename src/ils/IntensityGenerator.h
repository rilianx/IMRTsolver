/*
 * IntensityGenerator.h
 *
 *  Created on: Apr 23, 2019
 *      Author: iaraya
 */

#ifndef ILS_INTENSITYGENERATOR_H_
#define ILS_INTENSITYGENERATOR_H_

#include "Plan.h"

namespace imrt {

class IntensityGenerator {
public:
	typedef struct potencial_Gap{
		int GAP;
		int Suplier1; //biggest to the rigth
		int Suplier2; //biggest to the left
	}Gap;

	IntensityGenerator();

	void generate(Plan& P,double alpha, int max_intensity);
	void IntensityRepair(Plan& P);
	void changeworst(double* intensities, list<Gap> gaps);
	virtual ~IntensityGenerator();
};

} /* namespace imrt */

#endif /* ILS_INTENSITYGENERATOR_H_ */
