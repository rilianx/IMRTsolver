/*
 * Schedule.h
 *
 *  Created on: 7 may. 2018
 *      Author: iaraya
 */

#ifndef SCHEDULE_H_
#define SCHEDULE_H_

namespace imrt {


class Station {
	double angle;
	int ap_mid[N];
	int ap_ext[N];
};

/**
 *
 */

class Schedule {
public:
	Schedule();
	virtual ~Schedule();

	double eval();

private:
	Station stations;
};

} /* namespace imrt */

#endif /* SCHEDULE_H_ */
