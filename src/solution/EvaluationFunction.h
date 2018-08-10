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

class EvaluationFunction {
public:

	//Constructor of the evaluator.
	EvaluationFunction(vector<Volume>& volumes);
  
  EvaluationFunction(const EvaluationFunction& F);

	virtual ~EvaluationFunction();
	
	EvaluationFunction& operator=(const EvaluationFunction & ef);

	// Generate the dose distribution matrices Z for each organ
	void generate_Z(const Plan& p);

	// Eval the cost F based on the dose deposition matrix Z
	double eval(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	double incremental_eval(Station& station, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax,
			list< pair< int, double > >& diff);

 // double delta_eval(Station& station, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax, 
 //     list< pair< int, double > >& diff);

	void undo_last_eval(vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	void update_sorted_voxels(vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax, int o, int k, bool erase=true);
  //Additional functions

	// Generate files in the plotter directory with the voxel_dose functions for each organ
	void generate_voxel_dose_functions ();

	//Return the n beamlets with most impact in F taking into account the nv worst voxels.
	//mode=1: only considers beamlets with positive impact in F
	//mode=-1: only considers beamlets with negative impact in F
	//Each returned beamlet is a pair (eval, sign),(station, beamlet)
	set < pair< pair<double,bool>, pair<Station*, int> >,
	std::greater< pair< pair<double,bool>, pair<Station*, int> > > >
	best_beamlets(Plan& p, int n, int nv, int mode=0);

	void generate_linear_system(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	static int n_evaluations;
	
	int n_volumes;

private:
	//dose distribution vectors for each organ
	vector< vector<double> > Z;

	//voxel_dose[o][d] indicates the number of voxels in organ o with a dose between d and d+1
	vector<vector<double> > voxel_dose;

	//number of organs, including the tumor
	int nb_organs;

	//number of voxels for each organ
	vector<int> nb_voxels;

	double F;

    //Evaluation of F before the last incremental evaluation
	double prev_F;

  //Extra data

	// all the tumor voxels sorted by derivative (may be useful for algorithms)
  set< pair <double, pair<int,int> > > tumor_voxels;

  set< pair <double, pair<int,int> >, std::greater< pair <double, pair<int,int> > > > voxels;

  set< pair <double, pair<int,int> >, std::greater< pair <double, pair<int,int> > > > o_voxels;

	//Matrix of derivatives for each organ and voxel (may be useful for algorithms)
	//How much increase/decrease F increasing the voxe in one unity.
	vector< vector<double> > D;

	list < pair < pair<int,int>, double > > Z_diff;


};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
