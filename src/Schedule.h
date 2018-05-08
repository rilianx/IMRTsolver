/*
 * Schedule.h
 *
 *  Created on: 7 may. 2018
 *      Author: iaraya
 */

#ifndef SCHEDULE_H_
#define SCHEDULE_H_

#include "IMRTProblem.h"

namespace imrt {


/**
 *
 */

class Schedule {
public:
	Schedule();
	virtual ~Schedule();

	double eval(IMRT_Problem& p){
		
	}

private:
	Station stations;
};

} /* namespace imrt */

#endif /* SCHEDULE_H_ */
