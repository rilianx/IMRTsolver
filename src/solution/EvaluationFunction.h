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
     bool operator()(const double& lhs, const double& rhs)
     {
         return abs(lhs) > abs(rhs);
     }

     bool operator()(const pair< double, pair<int,int> >& lhs, const pair< double, pair<int,int> >& rhs)
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
      	vector<double>& w, vector<double>& Zmin, vector<double>& Zmax, int n_voxels=999999) const;

  double get_delta_eval(list< pair< int, double > >& diff, double angle,
                        vector<double>& w, vector<double>& Zmin, vector<double>& Zmax, int n_voxels=999999) const;

  double get_impact_beamlet(int angle, int b);

  //return a kind of evaluation of the beamlets: (impact o the tumor)/(impact to the organs)
  //unlike get_impact_beamlet, it is not based on the evaluation function (only Z and constraints)
  double get_ratio_beamlet(vector<double>& w,
  	vector<double>& Zmin, vector<double>& Zmax, int angle, int b);



	void undo_last_eval(vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	//Update D (partial derivative of F w.r.t. the voxels)
	//And resorts the set of voxels
	//Must be called every time that Z[o][k] changes
	void update_sorted_voxels(vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax, int o, int k);

  //Update the impacts og each beamlet affected by the voxel [o][k]
	void update_beamlets_impact(int o, int k, double prev_Dok=0.0);


  //Additional functions

	// Generate files in the plotter directory with the voxel_dose functions for each organ
	void generate_voxel_dose_functions ();

	//Return the beamlets sorted by impact on F taking into account the nv worst voxels.
	//Each returned beamlet is a pair eval,(station, beamlet)
  multimap < double, pair<int, int>, MagnitudeCompare >
  best_beamlets(Plan& p, int nv);

  //returns a map of beamlets sorted by their impact in F (derivative, (station, beamlet))
	multimap < double, pair<Station*, int>, MagnitudeCompare > get_sorted_beamlets(Plan& p);

	void generate_linear_system(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	//update Z by increasing the intensity of beamlet (angle,b) in delta_intensity
	//if return_if_unfeasible=true, then it returns when some organ voxel surpasses Zmax
	//return false if some voxel surpasses Zmax
	bool Zupdate(int angle, int b, double delta_intensity, bool return_if_unfeasible,
				 vector<double>& Zmax);

	//regresa al savepoint para Z
	void Zrollback();

	void Zsavepoint(){
		Z_diff.clear();
	}

	pair<double,double> get_value_cost(int angle, int b, vector<double>& Zmin, vector<double>& Zmax);

	void get_vc_sorted_beamlets(Plan& p, vector<double>& Zmin, vector<double>& Zmax,
			multimap < double, pair<Station*, int>, MagnitudeCompare >& sorted_set);

	static int n_evaluations;

	int n_volumes;

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

    //Evaluation of F before the last incremental evaluation
	double prev_F;

  //Extra data
	// all the tumor voxels sorted by derivative (may be useful for algorithms)
   set< pair< double, pair<int,int> >, MagnitudeCompare >  voxels;

  //Beamlets (angle,b) sorted by impact (partial derivative magnitude)
  //multimap < double, pair<int, int>, MagnitudeCompare > sorted_beamlets;
	unordered_map <pair<int, int>, double, pair_hash> beamlet_impact;

  //voxels (o,k) --> beamlet list (angle,b)
  unordered_map < pair<int, int>, list< pair<int,int> >, pair_hash > voxel2beamlet_list;

  //beamlet (angle,b) --> lista de voxels (o,k)
  unordered_map < pair<int, int>, multimap<double, pair<int,int> >, pair_hash > beamlet2voxel_list;

	//Matrix of derivatives for each organ and voxel (may be useful for algorithms)
	//How much increase/decrease F increasing the voxel in one unity.
	vector< vector<double> > D;

  //Store the changes in Z (dose distribution) after an incremental evaluation
	//It is used to undo the changes by using the undo_last_eval function
	list < pair < pair<int,int>, double > > Z_diff;


};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
