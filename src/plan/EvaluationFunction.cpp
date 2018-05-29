/*
 * EvaluationFunction.cpp
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include "EvaluationFunction.h"


namespace imrt {

EvaluationFunction::EvaluationFunction(vector<Volume>& volumes) :
	nb_organs(volumes.size()), nb_voxels(volumes.size()) {

	for(int i=0; i<nb_organs; i++){
		nb_voxels[i]=volumes[i].getNbVoxels();
		Z.insert(Z.end(), vector<double>(nb_voxels[i]));
	}
}

EvaluationFunction::~EvaluationFunction() { }


void EvaluationFunction::generate_Z(const Plan& p){
	const list<Station*>& stations=p.get_stations();


	for(auto station : stations){
		for(int o=0; o<nb_organs; o++)
		 	 std::fill(Z[o].begin(), Z[o].end(), 0.0);

		//considering 2*Xmid, Xext
		//we update the dose distribution matrices Z with the dose delivered by the station
		for(int o=0; o<nb_organs; o++){
				Matrix&  D = station->getDepositionMatrix(o);
				for(int k=0; k<nb_voxels[o]; k++){
			  double dose=0.0;
			  for(int b=0; b<station->getNbBeamlets(); b++){
					dose +=  D(k,b)*station->getIntensity( b );
				}
			  Z[o][k] += dose;
			}
		}
	}
}


double EvaluationFunction::eval(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax){
	generate_Z(p);

	double F=0.0;
  double max_pen=0.0;
	pair<int,int> worst_voxel=make_pair(-1,-1);
	for(int o=0; o<nb_organs; o++)
		for(int k=0; k<nb_voxels[o]; k++){
			double pen=0.0;
			if(Z[o][k] < Zmin[o] )
				 pen = w[o] * ( pow(Zmin[o]-Z[o][k], 2) );

			if(Z[o][k] > Zmax[o] )
				 double pen = w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
			F+= pen;

			if(pen>max_pen){
				max_pen=pen;
				worst_voxel=make_pair(o,k);
			}
		}

	return F;
}




} /* namespace imrt */
