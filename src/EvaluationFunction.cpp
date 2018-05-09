/*
 * EvaluationFunction.cpp
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include "EvaluationFunction.h"

namespace imrt {

EvaluationFunction::EvaluationFunction(int nb_organs, int nb_beamlets, vector<int>& nb_voxels, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax) :
	w(w), Zmin(Zmin), Zmax(Zmax), nb_organs(nb_organs), nb_voxels(nb_voxels), nb_beamlets(nb_beamlets), last_eval(0.0) {

	for(int i=0; i<nb_organs; i++)
		Z.insert(Z.end(), vector<double>(nb_voxels[i]));

}

EvaluationFunction::~EvaluationFunction() { }

//TODO: this should be automatically computed in the future
bool EvaluationFunction::set_deposition_matrix(double angle, int id_organ, string file){

	if(D.find(angle)==D.end())
		D[angle]=vector<Matrix>(nb_organs);

	string line;
	ifstream myfile (file);

	if (myfile.is_open()){
		//myfile.clear();
		//myfile.seekg(0, ios::beg);

		D[angle][id_organ] = Matrix(nb_voxels[id_organ],nb_beamlets);

		int i=0, j=0;
		getline (myfile,line);
		while (getline (myfile,line, '\t') )	{
			if(j==0){	j++; continue; }
			D[angle][id_organ](i,j-1)=atof(line.c_str());
			j++;
			if(j==nb_beamlets+1){j=0; i++; getline (myfile,line);}
		}

		myfile.close();

		return true;
	}

	return false;
}

void EvaluationFunction::generate_Z(const Plan& p){
	const list<Station*>& stations=p.stations;

	for(auto station : stations){

		if(station->beamlets.size() != nb_beamlets){
			cout << "the number of beamlets in the station ("<< station->beamlets.size() ;
			cout << ") does not match with the beamlets in the evaluation function (" << nb_beamlets << ")" << endl;
			exit(1);
		}


		if(D.find(station->get_angle()) == D.end()){
			//generate_deposition_matrices(angle, D[angle]);
		}else{
			for(int o=0; o<nb_organs; o++)
				std::fill(Z[o].begin(), Z[o].end(), 0.0);
		}



		//considering 2*Xmid, Xext
		//we update the dose distribution matrices Z with the dose delivered by the station


			for(int o=0; o<nb_organs; o++)
			for(int k=0; k<nb_voxels[o]; k++){
			  double dose=0.0;
			  for(int b=0; b<nb_beamlets; b++){
				if(station->beamlets[b])
					dose +=  D[station->get_angle()][o](k,b);

			  }
			  Z[o][k] += dose* station->get_intensity();
			}



		station->switched_lets.clear();
	}
}


double EvaluationFunction::eval(const Plan& p, bool genZ){
	if(genZ) generate_Z(p);

	double F=0.0;

	for(int o=0; o<nb_organs; o++)
		for(int k=0; k<nb_voxels[o]; k++){
			if(Z[o][k] < Zmin[o] ) F+= w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
			if(Z[o][k] > Zmax[o] ) F+= w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
		}

	last_eval=F;
	return F;
}


double EvaluationFunction::incremental_eval(Station& station){
	/* Another option is recomputing F with the modified Z instead of computing delta_F */

  double delta_F=0.0;

  //for each voxel we compute the change produced by the modified beamlets
  //while at the same time we compute the variation in the function F produced by all these changes

  for(int o=0; o<nb_organs; o++)
	for(int k=0; k<nb_voxels[o]; k++){

		//we compute the change in the delivered dose in voxel k of the organ o
		double delta=0.0;

	    double delta_st=0.0;
	    for(auto b : station.switched_lets){
				//we compute the change (delta) of the dose in the voxel k
		  if(station.beamlets[b]) //the beamlet turned on
			delta_st+=D[station.get_angle()][o](k, b);
		  else //the beamlet turned off
				delta_st-=D[station.get_angle()][o](k, b);
		}
		delta += delta_st*station.get_intensity();

		if(delta==0.0) continue; //no change in the voxel

		//with the change in the dose of a voxel we can incrementally modify the value of F
		if(Z[o][k] < Zmin[o] && Z[o][k] + delta < Zmin[o]) //update the penalty
			delta_F += w[o]*delta*(delta+2*(Z[o][k]-Zmin[o]));
		else if(Z[o][k] < Zmin[o]) //the penalty disappears
			delta_F -=  w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
		else if(Z[o][k] + delta < Zmin[o]) //the penalty appears
			delta_F +=  w[o] * ( pow(Zmin[o]-(Z[o][k]+delta), 2) );

		if(Z[o][k] > Zmax[o] && Z[o][k] + delta > Zmax[o]) //update the penalty
			delta_F += w[o]*delta*(delta+2*(Zmax[o] - Z[o][k]));
		else if(Z[o][k] > Zmax[o]) //the penalty disappears
			delta_F -=  w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
		else if(Z[o][k] + delta > Zmax[o]) //the penalty appears
			delta_F +=  w[o] * ( pow(Z[o][k]+delta - Zmax[o], 2) );

		Z[o][k]+=delta;
	}

  station.switched_lets.clear();

  last_eval+=delta_F;
  return last_eval;

  //return eval(p, false);

}


} /* namespace imrt */
