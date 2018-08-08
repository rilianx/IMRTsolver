/*
 * IntensityILS.h
 *
 *  Created on: 2 ago. 2018
 *      Author: iaraya
 */

#include "ILS.h"

#ifndef ILS_INTENSITYILS_H_
#define ILS_INTENSITYILS_H_

namespace imrt {

class IntensityILS  : public ILS {
public:
	IntensityILS(int step_intensity, int bsize, int vsize, int maxdelta, int maxratio, double alpha, double beta, int perturbation_size=0) :
		ILS(bsize, vsize), step_intensity(step_intensity), maxdelta(maxdelta), maxdelta0(maxdelta), maxratio0(maxratio), maxratio(maxratio), alpha(alpha), beta(beta),
		perturbation_size(perturbation_size){ };
	virtual ~IntensityILS() { }

	virtual double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P);

	virtual bool acceptanceCriterion(double new_eval, double prev_eval){
		return false;
	}

	virtual void undoLast(Plan& p){
		p.undoLast2();
	}

   // Consiste en una reduccion de las aperturas en matrices de intensidad
   // seleccionadas aleatoriamente

	virtual double perturbation(Plan& P) {
		double eval=P.getEvaluation();
		for(int i=0; i<perturbation_size; i++){
			int r=rand()%P.get_stations().size();
			Station* s=P.get_station(r);
			list< pair< int, double > > diff;
			s->reduce_apertures(diff);
			eval=P.incremental_eval(*s,diff);
		}
		maxdelta=maxdelta0;
		maxratio=maxratio0;

		return(eval);
	};

	virtual bool perturbate(int no_improvement, int iteration) {

		return (perturbation_size>0 && no_improvement>100);
	};

	private:

	int step_intensity;
	int perturbation_size;

	double maxdelta0;
	double maxratio0;

	double maxdelta;
	double maxratio;
	double alpha;
	double beta;



};

}
#endif /* SRC_ILS_INTENSITYILS_H_ */
