/*
 * IntensityGenerator.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: iaraya
 */
#include <algorithm>

#include "IntensityGenerator.h"
#include "Plan.h"
#include "Station.h"


namespace imrt {

void IntensityGenerator::generate(Plan& P, double alpha, int max_intensity){
	Plan Pprime(P);
	cout <<"OldPlan:"<< Pprime.eval() << endl  ;
	Start:
		Pprime.Zsavepoint();
		multimap < double, pair<Station*, int>, MagnitudeCompare > sorted_beamlets;
		Pprime.get_vc_sorted_beamlets(Pprime,sorted_beamlets); //ordeandos v/cx
		for (auto elem:sorted_beamlets){
			Station* s = elem.second.first;
			int b = elem.second.second;
			pair<double,double> vc=Pprime.get_value_cost(s,b);

			int i=s->beam2pos[b].first; int j=s->beam2pos[b].second;
			s->change_intensity(i, j, s->I(i,j) + std::min(int(alpha/vc.second+0.5),max_intensity));

			//Pprime.Zupdate(s,b,int(alpha/vc.second+0.5), false);
			/*
			if ( !Pprime.Zupdate(s,b,int(alpha/vc.second+0.5), false)){
				alpha = alpha * 0.9;//alpha decrease
				Pprime.Zrollback();
				goto Start;
				//cout << alpha << endl;
			}*/
		}

    //se modifican las intensidades en la solucion (arreglar)

		cout << alpha << endl;
		cout <<"NewPlan:"<< Pprime.eval() << endl  ;
		cout << endl;
  		for(int i=0;i<5;i++)
  		Pprime.printIntensity(i);
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
