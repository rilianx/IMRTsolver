/*
 * EvaluatorF.h
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
#include "EvaluationFunction.h"

#ifndef Evaluator_H_
#define Evaluator_H_

using namespace std;
using namespace maths;


namespace imrt {


class FluenceMap{

	FluenceMap(vector<Volume>& volumes, const Collimator& collimator);
	virtual ~FluenceMap();
	void create_voxel2beamlet_list(vector<Volume>& volumes, const Collimator& collimator);

	// Generate the dose distribution matrices Z for each organ
	void computeFM(const Plan& p);

	//Get deltaZ related to the list of changes
	std::vector<list < pair<int, double>>>& compute_deltaFM(list< pair< int, double > >& changes, double angle) const;
	std::vector<list < pair<int, double>>>& get_deltaFM() const {return deltaFM;}
	void updateFM(list< pair< int, double > >& changes, double angle, std::vector<list < pair<int, double>>>& deltaFM);

private:
	//dose distribution vectors for each organ
	vector< vector<double> > FM;
	mutable std::vector<list < pair<int, double>>> deltaFM;

	//voxel_dose[o][d] indicates the number of voxels in organ o with a dose between d and d+1
	vector<vector<double> > voxel_dose;

  	vector<Volume>& volumes;

	//number of organs, including the tumor
	int nb_organs;

	//number of voxels for each organ
	vector<int> nb_voxels;

  	//voxels[angle][o,k] --> beamlet list (b)
  	unordered_map < int, unordered_map < pair<int, int>, list <int> , pair_hash > > voxel2beamlet_list;

  	//beamlet [angle][b] --> lista de voxels (o,k)
  	unordered_map < int, unordered_map <int, multimap<double, pair<int,int> > > > beamlet2voxel_list;

};


/* Abstract Evaluator
 */


class Evaluator{


public:

	//Constructor of the evaluator.
	Evaluator(EvaluationFunction& z_structure, vector<double>& w, vector<double>& Zmin, 
    vector<double>& Zmax) : z_structure(z_structure), Zmin(Zmin), Zmax(Zmax), w(w), Z(z_structure.get_Z()) { };


    EvaluationFunction& z_structure;
    vector< vector<double> >& Z;
    vector<double>& w;
    vector<double>& Zmin;
    vector<double>& Zmax;


	virtual ~Evaluator() { }

	// Eval the cost F based on the dose deposition matrix Z
	virtual double eval(const Plan& p)=0;

	virtual double incremental_eval(list< pair< int, double > >& changes, double angle)=0;

    double get_delta_eval(int angle, int b, double delta_intensity) const{
			  list< pair< int, double > > changes;
			  changes.push_back(make_pair(b,delta_intensity));
			  return  get_delta_eval(changes, angle);
	}

    virtual double get_delta_eval(list< pair< int, double > >& changes, double angle) const=0;

	vector< vector<double> >& get_Z(){
		return Z;
	};


};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
