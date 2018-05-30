/*
 * EvaluationFunction.cpp
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include "EvaluationFunction.h"
#include <string>

namespace imrt {

EvaluationFunction::EvaluationFunction(vector<Volume>& volumes) :
	nb_organs(volumes.size()), nb_voxels(volumes.size()), voxel_dose(volumes.size(), vector<double>(100)) {

	for(int i=0; i<nb_organs; i++){
		nb_voxels[i]=volumes[i].getNbVoxels();
		Z.insert(Z.end(), vector<double>(nb_voxels[i]));
	}
}

EvaluationFunction::~EvaluationFunction() { }


void EvaluationFunction::generate_Z(const Plan& p){
	const list<Station*>& stations=p.get_stations();

	for(int o=0; o<nb_organs; o++)
	 	 std::fill(Z[o].begin(), Z[o].end(), 0.0);

	for(auto station : stations){


		//considering 2*Xmid, Xext
		//we update the dose distribution matrices Z with the dose delivered by the station
		for(int o=0; o<nb_organs; o++){
			 const Matrix&  D = station->getDepositionMatrix(o);
			 for(int k=0; k<nb_voxels[o]; k++){
			   double dose=0.0;
			   for(int b=0; b<station->getNbBeamlets(); b++){
				   if(station->getIntensity( b )>0.0)
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

	worst_voxel=make_pair(-1,-1);
	for(int o=0; o<nb_organs; o++)
		for(int k=0; k<nb_voxels[o]; k++){
			double pen=0.0;
			if(Z[o][k] < Zmin[o] )
				 pen = w[o] * ( pow(Zmin[o]-Z[o][k], 2) );

			if(Z[o][k] > Zmax[o] )
				 pen = w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
			F+= pen;

			if(pen>max_pen && black_list.find(make_pair(o,k)) == black_list.end()){
				max_pen=pen;
				worst_voxel=make_pair(o,k);
			}
		}

	cout << max_pen << endl;

	return F;
}

pair<int,int> EvaluationFunction::get_worst_voxel(){
	return worst_voxel;
}


int EvaluationFunction::best_beamlet(const Station& s, pair<int,int>& voxel){
	const Matrix&  D = s.getDepositionMatrix(voxel.first);
	int k=voxel.second;
	int bestb=-1; double max_dose=0.0;

	for(int b=0; b<s.getNbBeamlets(); b++){
		double dose=D(k,b);
		if(dose>max_dose){
			max_dose=dose;
			bestb=b;
		}
	}


	if(bestb==-1)
		black_list.insert(voxel);


	return bestb;
}

void EvaluationFunction::generate_voxel_dose_functions (){
	for(int o=0; o<nb_organs; o++){
		std::fill(voxel_dose[o].begin(), voxel_dose[o].end(), 0.0);
		for(int k=0; k<nb_voxels[o]; k++){
			if(Z[o][k]<100)
				voxel_dose[o][(int) Z[o][k]]+=1;
		}
	}


	for(int o=0; o<nb_organs; o++){
		ofstream myfile;
		myfile.open ("plotter/organ"+std::to_string(o)+".txt");

		double cum=0.0;
		for(int k=99; k>=0; k--){
			cum+= voxel_dose[o][k];
			myfile << k+1 << "," << cum << endl;
		}

		myfile.close();
	}
}




} /* namespace imrt */
