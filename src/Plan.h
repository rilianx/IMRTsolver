/*
 * Schedule.h
 *
 *  Created on: 7 may. 2018
 *      Author: iaraya
 */

#ifndef PLAN_H_
#define PLAN_H_

#include "Station.h"

namespace imrt {


/* An IMRT plan
 *
 * It consists in a list of stations, each station corresponds to an
 * angle, an aperture and an intensity
 */

class Plan {
public:
	Plan() {};
	virtual ~Plan() {};

	// Adds a new station to the plan
	void add_station(Station& s){
		stations.push_back(&s);
	}

	//The list of stations
	list<Station*> stations;
};

} /* namespace imrt */

#endif /* PLAN_H_ */
