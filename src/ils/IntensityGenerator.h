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
	IntensityGenerator();

	void generate(Plan& P);

	virtual ~IntensityGenerator();
};

} /* namespace imrt */

#endif /* ILS_INTENSITYGENERATOR_H_ */
