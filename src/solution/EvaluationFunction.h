/*
 * EvaluationFunction.h
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <iterator>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <queue>
#include <functional>

#include "Plan.h"
#include "Volume.h"
#include "Matrix.h"
#include "Station.h"

#ifndef EVALUATIONFUNCTION_H_
#define EVALUATIONFUNCTION_H_

using namespace std;
using namespace maths;


namespace imrt {

class Plan;

/* The Evaluation Function
 *
 * The evaluation function based on [this paper](https://drive.google.com/file/d/1YfMNk4GhBK97gSQ0nvpJAnyM6A3EPv61/view).
 * Basically it penalizes voxels with doses larger than Zmax (healthy organs), and with doses smaller than Zmin (tumor)
 */

struct pair_hash
{
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};


 struct MagnitudeCompare
 {
     bool operator()(const double& lhs, const double& rhs) const
     {
         return abs(lhs) > abs(rhs);
     }

     bool operator()(const pair< double, pair<int,int> >& lhs, const pair< double, pair<int,int> >& rhs) const
     {
    	 if(abs(lhs.first) > abs(rhs.first)) return true;
    	 if(abs(lhs.first) == abs(rhs.first)){
    		 //if(pair_hash()(lhs.second)< pair_hash()(rhs.second)) return true;
    		 if(lhs.second.first < rhs.second.first) return true;
    		 if(lhs.second.first == rhs.second.first &&  lhs.second.second < rhs.second.second) return true;
    	 }
         return false;
     }
};


class EvaluationFunction {

private: //singleton implementation
  static EvaluationFunction* instance;
	EvaluationFunction(EvaluationFunction& e) : volumes(e.volumes) {};

	//Constructor of the evaluator.
	EvaluationFunction(vector<Volume>& volumes, const Collimator& collimator);

public:

  //EvaluationFunction(const EvaluationFunction& F);

	/* Static access method. */
	static EvaluationFunction& getInstance(vector<Volume>& volumes, const Collimator& collimator){
		if(!instance) instance = new EvaluationFunction(volumes, collimator);
		return *instance;
	}

	static EvaluationFunction& getInstance(){
		return *instance;
	}

  void create_voxel2beamlet_list(vector<Volume>& volumes, const Collimator& collimator);

	virtual ~EvaluationFunction();

	// Generate the dose distribution matrices Z for each organ
	void generate_Z(const Plan& p);

	// Eval the cost F based on the dose deposition matrix Z
	double eval(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	double incremental_eval(Station& station, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax,
			list< pair< int, double > >& diff);

  double get_delta_eval(int angle, int b, double delta_intensity,
      	const vector<double>& w, const vector<double>& Zmin, const vector<double>& Zmax, int n_voxels=999999) const;

  double get_delta_eval(list< pair< int, double > >& diff, double angle,
                        const vector<double>& w, const vector<double>& Zmin, const vector<double>& Zmax, int n_voxels=999999) const;

	//Update D (partial derivative of F w.r.t. the voxels)
	//And resorts the set of voxels
	//Must be called every time that Z[o][k] changes
	void update_sorted_voxels(vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax, int o, int k);


  //Additional functions

	// Generate files in the plotter directory with the voxel_dose functions for each organ
	void generate_voxel_dose_functions ();

	//Return the beamlets sorted by impact on F taking into account the nv worst voxels.
	//Each returned beamlet is a pair eval,(station, beamlet)
  	multimap < double, pair<int, int>, MagnitudeCompare > sorted_beamlets(Plan& p, double vsize=0.01);


	void generate_linear_system(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	static int n_evaluations;

	int n_volumes;

	vector< vector<double> >& get_Z(){
		return Z;
	};

private:
	//dose distribution vectors for each organ
	vector< vector<double> > Z;

	//voxel_dose[o][d] indicates the number of voxels in organ o with a dose between d and d+1
	vector<vector<double> > voxel_dose;

  vector<Volume>& volumes;

	//number of organs, including the tumor
	int nb_organs;

	//number of voxels for each organ
	vector<int> nb_voxels;

	double F;

  //Extra data
	// all the tumor voxels sorted by derivative (may be useful for algorithms)
	// TODO:DEBERIA SER MULTISET??
   set< pair< double, pair<int,int> >, MagnitudeCompare >  voxels;

  //voxels[angle][o,k] --> beamlet list (b)
  unordered_map < int, unordered_map < pair<int, int>, list <int> , pair_hash > > voxel2beamlet_list;

  //beamlet [angle][b] --> lista de voxels (o,k)
  unordered_map < int, unordered_map <int, multimap<double, pair<int,int> > > > beamlet2voxel_list;

	//Matrix of derivatives for each organ and voxel (may be useful for algorithms)
	//How much increase/decrease F increasing the voxel in one unity.
	vector< vector<double> > D;




};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
