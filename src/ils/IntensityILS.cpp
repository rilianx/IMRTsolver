/*
 * IntensityILS.cpp
 *
 *  Created on: 2 ago. 2018
 *      Author: iaraya
 */

#include "IntensityILS.h"

namespace imrt {

	double IntensityILS::localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P){
		Station*s = target_beam.second.first; int beamlet=target_beam.second.second;
		bool sign=target_beam.first; //impact in F (+ or -)

		//double delta_intensity= rand()%3+1;

		double delta_intensity= (maxdelta>0.5)? rand()%int(maxdelta + 0.5) : 1;
		maxdelta = maxdelta*alpha;

		if(sign) delta_intensity*=-1;

		double ratio= (maxratio>0.5)? rand()%int(maxratio + 0.5) : 0;
		maxratio = maxratio*beta;

		//cout << maxdelta << "," << maxratio << endl;

		auto diff=s->increaseIntensity_repair(beamlet,delta_intensity,ratio);
		double eval=P.incremental_eval(*s,diff);
		//F.incremental_eval(*s,w,Zmin,Zmax, diff);

		return eval;
	}
}
