/*
 * Schedule.h
 *
 *  Created on: 7 may. 2018
 *      Author: iaraya
 */

#ifndef PLAN_H_
#define PLAN_H_

#include "EvaluationFunction.h"
#include "Station.h"

namespace imrt {


/* An IMRT plan
 *
 * It consists in a list of stations, each station corresponds to an
 * angle, and a set of apertures (or a matrix of intensities)
 */

class Plan {
public:
	Plan(EvaluationFunction &ev) : ev(ev) {};
	virtual ~Plan() {};

	// Adds a new station to the plan
	void add_station(Station& s){
		stations.push_back(&s);
	}

	double eval(vector<double>& w, vector<double>& Zmin, vector<double>& Zmax){
		return ev.eval(*this,w,Zmin,Zmax);
	}

	const list<Station*>& get_stations() const{
		return stations;
	}

private:
	//The list of stations
	list<Station*> stations;

	EvaluationFunction& ev;
};

} /* namespace imrt */

#endif /* PLAN_H_ */
