/*
 * EvaluationFunction.cpp
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include "EvaluationFunction.h"
#include <string>

namespace imrt {

EvaluationFunction::EvaluationFunction(vector<Volume>& volumes) : last_F(0.0),
	nb_organs(volumes.size()), nb_voxels(volumes.size()), voxel_dose(volumes.size(), vector<double>(100)) {

	for(int i=0; i<nb_organs; i++){
		nb_voxels[i]=volumes[i].getNbVoxels();
		Z.insert(Z.end(), vector<double>(nb_voxels[i]));
		P.insert(P.end(), vector<double>(nb_voxels[i]));
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
			   for(int b=0; b<station->getNbBeamlets(); b++)
					    dose +=  D(k,b)*station->getIntensity( b );
			   Z[o][k] += dose;
			 }
		}
		station->changed_lets.clear();
	}
}


double EvaluationFunction::eval(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax){
	sorted_voxels = priority_queue< pair <double, pair<int,int> > >();

	for(int o=0; o<nb_organs; o++)
	 	 std::fill(P[o].begin(), P[o].end(), 0.0);

	generate_Z(p);

	double F=0.0;

	for(int o=0; o<nb_organs; o++)
		for(int k=0; k<nb_voxels[o]; k++){
			double pen=0.0;
			if(Z[o][k] < Zmin[o] )
				 pen = w[o] * ( pow(Zmin[o]-Z[o][k], 2) );

			if(Z[o][k] > Zmax[o] )
				 pen = w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
			F+= pen;

			sorted_voxels.push(make_pair(pen,make_pair(o,k)));
			P[o][k]=pen;
		}

	last_F=F;
	return F;
}

double EvaluationFunction::incremental_eval(Station& station, vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax){
	/* Another option is recomputing F with the modified Z instead of computing delta_F */

  double delta_F=0.0;

  //for each voxel we compute the change produced by the modified beamlets
  //while at the same time we compute the variation in the function F produced by all these changes

  for(int o=0; o<nb_organs; o++){
			const Matrix&  D = station.getDepositionMatrix(o);
	for(int k=0; k<nb_voxels[o]; k++){

		//we compute the change in the delivered dose in voxel k of the organ o
		double delta=0.0;

		//cout << station.changed_lets.size() << endl;
		for (auto let:station.changed_lets){
		    int b=station.pos2beam[let.first];
			  if(D(k,b)==0.0) continue;
				delta+= D(k,b)*let.second;
		}


		if(delta==0.0) continue; //no change in the voxel

		double pen=0.0;
		//with the change in the dose of a voxel we can incrementally modify the value of F
		if(Z[o][k] < Zmin[o] && Z[o][k] + delta < Zmin[o]) //update the penalty
			pen += w[o]*delta*(delta+2*(Z[o][k]-Zmin[o]));
		else if(Z[o][k] < Zmin[o]) //the penalty disappears
			pen -=  w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
		else if(Z[o][k] + delta < Zmin[o]) //the penalty appears
			pen +=  w[o] * ( pow(Zmin[o]-(Z[o][k]+delta), 2) );

		if(Z[o][k] > Zmax[o] && Z[o][k] + delta > Zmax[o]) //update the penalty
			pen += w[o]*delta*(delta+2*(-Zmax[o] + Z[o][k]));
		else if(Z[o][k] > Zmax[o]) //the penalty disappears
			pen -=  w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
		else if(Z[o][k] + delta > Zmax[o]) //the penalty appears
			pen +=  w[o] * ( pow(Z[o][k]+delta - Zmax[o], 2) );

		delta_F += pen;
		P[o][k] +=pen;
		sorted_voxels.push(make_pair(P[o][k],make_pair(o,k)));

		Z[o][k]+=delta;
	}
  }


  last_F+=delta_F;
  return last_F;

  //return eval(p, false);

}

pair<int,int> EvaluationFunction::get_worst_voxel(){
	while(sorted_voxels.size()>0){
		 auto voxel=sorted_voxels.top();
		 if(voxel.first != P[voxel.second.first][voxel.second.second])
			 	sorted_voxels.pop();
			else return voxel.second;
	}
}

void EvaluationFunction::pop_worst_voxel(){
	sorted_voxels.pop();
}


int EvaluationFunction::max_beamlet_dose(const Station& s, pair<int,int>& voxel){
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
