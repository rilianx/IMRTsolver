/*
 * SearchStrategy.h
 *
 *  Created on: 01-06-2017
 *      Author: ignacio
 */

#include <time.h>
#include <stdio.h>
#include <iostream>
#include <list>
#include "Schedule.h"

#ifndef SEARCHSTRATEGY_H_
#define SEARCHSTRATEGY_H_

using namespace std;

namespace imrt {


class SearchStrategy {
public:
	SearchStrategy() : timelimit(0.0), begin_time(clock()), best_shedule(NULL) {} ;

	virtual ~SearchStrategy() {

	}

	double get_time(){
		return (double(clock()-begin_time)/double(CLOCKS_PER_SEC));
    }

	/*
	 * Initialize the variables of the specific strategy
	 */
	virtual void initialize()=0;

	/**
	 * Run the strategy
	 */
	virtual double run(double timelimit=99999.9, clock_t bt=clock())=0;

	virtual double get_best_value() const {
		if(best_shedule)
			return best_shedule->get_value();
		else return 0;
	}

	virtual const Schedule* get_best_shedule() {
		return best_shedule;
	}


private:
	double timelimit;
	Schedule* best_shedule;
	clock_t begin_time;


};

} /* namespace clp */

#endif /* SEARCHSTRATEGY_H_ */
