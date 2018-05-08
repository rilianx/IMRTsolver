/*
 * EvaluationFunction.h
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include <map>
#include <vector>
#include <list>
#include <iterator>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "Matrix.h"
#include "Station.h"

#ifndef EVALUATIONFUNCTION_H_
#define EVALUATIONFUNCTION_H_

using namespace std;
using namespace maths;


namespace imrt {

class EvaluationFunction {
public:

	/**
	 * Constructor of the evaluator.
	 * w_min[tumor] is the penalization if we deliver a dose lower than Zmin[tumor] in the tumor
	 * w_max[organ] is the penalization if we deliver a dose greater than Zmax[organ] in the organ
	 */

	EvaluationFunction(int nb_organs, int nb_cols, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax) :
		w(w), Zmin(Zmin), Zmax(Zmax), nb_organs(nb_organs), nb_rows(nb_organs), nb_cols(nb_cols), last_eval(0.0) {

		for(int i=0; i<nb_organs; i++)
			Z.insert(Z.begin(), vector<double>(nb_rows[i]));

	}

	virtual ~EvaluationFunction() { }

	bool set_deposition_matrix(double angle, int id_organ, string file){
		string line;
		ifstream myfile (file);

		if (myfile.is_open()){
			nb_rows[id_organ]=std::count(std::istreambuf_iterator<char>(myfile),
			             std::istreambuf_iterator<char>(), '\n')-1;
			myfile.clear();
			myfile.seekg(0, ios::beg);
			D[angle][id_organ]=Matrix(nb_rows[id_organ],nb_cols);

			cout << nb_rows[id_organ] << " rows" << endl;

			int i=0, j=0;
			getline (myfile,line);
			while (getline (myfile,line, '\t') )	{
				if(j==0){	j++; continue; }
				D[angle][id_organ](i,j-1)=atof(line.c_str());
				j++;
				if(j==nb_cols+1){j=0; i++;}
			}

			myfile.close();

			return true;
		}

		return false;
	}

	/*
	 * Generate the dose distribution matrices Z for each organ
	 */
	void generate_Z(const list<Station>& stations){

		for(auto station : stations){
			if(D.find(station.get_angle()) == D.end()){
				//generate_deposition_matrices(angle, D[angle]);
			}else{
				for(int o=0; o<nb_organs; o++)
					std::fill(Z[o].begin(), Z[o].end(), 0.0);
			}

			//considering 2*Xmid, Xext
			//we update the dose distribution matrices Z with the dose delivered by the station
			for(int o=0; o<nb_organs; o++)
			for(int k=0; k<nb_rows[o]; k++)
			for(int b=0; b<nb_cols; b++)
				Z[o][k] +=  get_dose_voxel(k, D[station.get_angle()][o], station);


		}
	}

	/**
	 * Eval the cost F based on the dose deposition matrix Z
	 */
	double eval(const list<Station>& stations){
		generate_Z(stations);

		double F=0.0;

		for(int o=0; o<nb_organs; o++)
			for(int k=0; k<nb_rows[o]; k++){
				if(Z[o][k] < Zmin[o] ) F+= w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
				if(Z[o][k] > Zmax[o] ) F+= w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
			}

		last_eval=F;
		return F;
	}

	double incremental_eval(const list<Station>& stations){

		last_eval += get_delta_F(stations);
		return last_eval;

	}

	double get_delta_F(const list<Station>& stations, bool updateZ=true){
	  double delta_F=0.0;

	  //for each voxel we compute the change produced by the modified beamlets
	  //while at the same time we compute the variation in the function F produced by all these changes

	  for(int o=0; o<nb_organs; o++)
		for(int k=0; k<nb_rows[o]; k++){

			//we compute the change in the delivered dose in voxel k of the organ o
			double delta=0.0;
			for(auto station : stations)
			  for(auto change : station.last_changes){
				int b=change.first;
				if(change.second != station.beamlets[b]) continue; //no change in the beamlet

					//we compute the change (delta) of the dose in the voxel k
				if(change.second) //the beamlet turned on
					delta+=D[station.get_angle()][o](k, b);
				else if(!change.second) //the beamlet turned off
					delta-=D[station.get_angle()][o](k, b);


			  }
			if(delta==0.0) continue; //no change in the voxel

			//with the change in the dose of a voxel we can incrementally modify the value of F

			if(Z[o][k] < Zmin[o] && Z[o][k] + delta < Zmin[o]) //update the penalty
				delta_F += w[o]*delta*(delta+2*(Z[o][k]-Zmin[o]));
			else if(Z[o][k] < Zmin[o]) //the penalty disappears
				delta_F -=  w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
			else if(Z[o][k] + delta < Zmin[o]) //the penalty appears
				delta_F +=  w[o] * ( pow(Zmin[o]-Z[o][k]+delta, 2) );

			if(Z[o][k] > Zmax[o] && Z[o][k] + delta > Zmax[o]) //update the penalty
				delta_F += w[o]*delta*(delta+2*(Zmax[o] - Z[o][k]));
			else if(Z[o][k] > Zmax[o]) //the penalty disappears
				delta_F -=  w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
			else if(Z[o][k] + delta > Zmax[o]) //the penalty appears
				delta_F +=  w[o] * ( pow(Z[o][k]+delta - Zmax[o], 2) );

			if(updateZ) Z[o][k]+=delta;
			}

	  return delta_F;
	}

	/**
	 * dose delivered to the voxel k of the organ by the given station
	 */
	double get_dose_voxel(int k, Matrix& d, const Station& station){
		double dose=0.0;
		for(int b=0; b<nb_cols; b++)
			if(station.beamlets[b])	dose+=d(k, b);
		return dose;
	}

private:
	/*
	 * Dose deposition matrices for each angle and organ
	 * D[a][organ_id](k,b): Dose delivered to voxel k of the organ by the beamlet b and angle a
	 */
	map<double, vector<Matrix>> D;

	/*
	 * Dose distribution vectors for each organ
	 */
	vector< vector<double> > Z;

	/* penalizations for delivering less (w_min) or greater (w_max) than the dose required for each organ*/
	vector<double> w;

	/* minimum and maximum dose required by each organ. Zmax=0 for the tumor and Zmin=0 for the organs */
	vector<double> Zmin;
	vector<double> Zmax;

	double last_eval;

	int nb_organs;
	int nb_cols;
	vector<int> nb_rows;
};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
