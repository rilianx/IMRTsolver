/*
 * IntensityILS.cpp
 *
 *  Created on: 2 ago. 2018
 *      Author: iaraya
 */
 #include <algorithm>    // std::random_shuffle

#include "IntensityILS.h"

namespace imrt {

	// random generator function:
	int myrandom (int i) { return std::rand()%i;}

   vector < pair<int, int> > IntensityILS::getShuffledCells(Station* s){
		 vector < pair<int, int> > cells;
		 for (int i=0; i<s->I.nb_rows();i++)
			 for (int j=0; j<s->I.nb_cols(); j++)
				 if(s->I(i,j)!=-1) cells.push_back(make_pair(i,j));
		 std::random_shuffle ( cells.begin(), cells.end(), myrandom);
		 return cells;
	 }

   double IntensityILS::iLocalSearch(Plan& P, bool verbose){
		 const list<Station*> stations=P.get_stations();
		 list< pair< int, double > > diff;
		 double eval =P.getEvaluation();

		 for(auto s:stations){
			  vector < pair<int, int> > cells = getShuffledCells(s);
				for(auto cell:cells){
					int a=rand()%2;
					for(int b=0; b<2;b++){
					  int i=cell.first;
					  int j=cell.second;
						if(a==b) diff=s->intensityUp(i,j);
						else diff=s->intensityDown(i,j);
						if(diff.empty()) continue;

						double new_eval=P.incremental_eval(*s,diff);
						if(new_eval<eval)
							return eval;
						undoLast(P);
					}
				}
		 }
	 }

	double IntensityILS::localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P){
		Station*s = target_beam.second.first; int beamlet=target_beam.second.second;
		bool sign=target_beam.first; //impact in F (+ or -)

		//double delta_intensity= rand()%3+1;

		double delta_intensity= rand()%int(maxdelta-step_intensity+1)+step_intensity; //random entre step_intensity y maxdelta
		delta_intensity = (int)  (delta_intensity/step_intensity) * step_intensity;

		maxdelta = maxdelta*alpha;
		if(maxdelta < step_intensity) maxdelta=step_intensity;

		if(sign) delta_intensity*=-1;

		double ratio= (maxratio>0.5)? rand()%int(maxratio + 0.5) : 0;
		maxratio = maxratio*beta;

		auto diff=s->increaseIntensity_repair(beamlet,delta_intensity,ratio);
		double eval=P.incremental_eval(*s,diff);
		//F.incremental_eval(*s,w,Zmin,Zmax, diff);

		return eval;
	}
}
