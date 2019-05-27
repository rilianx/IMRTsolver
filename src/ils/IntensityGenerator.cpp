/*
 * IntensityGenerator.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: iaraya
 */

#include "IntensityGenerator.h"
#include "Plan.h"
#include "Station.h"


namespace imrt {

void IntensityGenerator::generate(Plan& P, double alpha){
	Plan Pprime(P);
	cout <<"OldPlan:"<< Pprime.eval() << endl  ;
	Start:	
		Pprime.Zsavepoint();
		multimap < double, pair<Station*, int>, MagnitudeCompare > sorted_beamlets;
		Pprime.get_vc_sorted_beamlets(sorted_beamlets);
		for (auto elem:sorted_beamlets){
			Station* s = elem.second.first;
			int b = elem.second.second;
			pair<double,double> vc=Pprime.get_value_cost(s,b);
			if (Pprime.Zupdate(s,b,int(alpha/vc.second+0.5), true)){
				alpha = alpha * 0.9;//alpha decrease
				Pprime.Zrollback();
				goto Start;
			}
		} 
		cout <<"NewPlan:"<< Pprime.eval() << endl  ;
		cout << endl;
  		for(int i=0;i<5;i++)
  		P.printIntensity(i);
  		cout << endl;

}
// preguntar return_if_unfeasible
// delta_intensity dentro de plan.h es el equivalente al alpha?
// Por que se definieron las funciones el .h y no en el .cpp?

IntensityGenerator::IntensityGenerator() {
	// TODO Auto-generated constructor stub

}

IntensityGenerator::~IntensityGenerator() {
	// TODO Auto-generated destructor stub
}

} /* namespace imrt */
