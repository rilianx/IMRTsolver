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
#include "Evaluator.h"

#ifndef EvaluatorF_H_
#define EvaluatorF_H_

using namespace std;
using namespace maths;


namespace imrt {

class Plan;

 struct MagnitudeCompare2
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

/* The Evaluator F
 *
 * The evaluation function based on [this paper](https://drive.google.com/file/d/1YfMNk4GhBK97gSQ0nvpJAnyM6A3EPv61/view).
 * Basically it penalizes voxels with doses larger than Zmax (healthy organs), and with doses smaller than Zmin (tumor)
 */


class EvaluatorF : public Evaluator{



public:


private: //singleton implementation
    static EvaluatorF* instance;
	EvaluatorF(EvaluatorF& e) : Evaluator(e.z_structure,e.w,e.Zmin,e.Zmax) {}; //no copy

	//Constructor of the evaluator.
	EvaluatorF(EvaluationFunction& z_structure, vector<double>& w, 
    vector<double>& Zmin, vector<double>& Zmax) : Evaluator(z_structure,w,Zmin,Zmax), F(0.0){
		for(int i=0; i<Z.size(); i++)
			D.insert(D.end(), vector<double>(Z[i].size()));
		
	 };

public:

	/* Static access method. */
	static EvaluatorF& getInstance(EvaluationFunction& z_structure, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax){
		if(!instance) instance = new EvaluatorF(z_structure,w,Zmin,Zmax);
		return *instance;
	}

	static EvaluatorF& getInstance(){
		return *instance;
	}


	virtual ~EvaluatorF() { }


	// Eval the cost F based on the dose deposition matrix Z
	virtual double eval(const Plan& p);

	virtual double incremental_eval(list< pair< int, double > >& changes, double angle);

    virtual double get_delta_eval(list< pair< int, double > >& changes, double angle) const;

	virtual double get_evaluation(){ return F; }



  //Additional functions

    //Return the beamlets sorted by impact on F taking into account the nv worst voxels.
    //Each returned beamlet is a pair eval,(station, beamlet)
    multimap < double, pair<int, int>, MagnitudeCompare2 > sorted_beamlets(const Plan& p, double vsize);

    void update_sorted_voxels(int o, int k);

	static int n_evaluations;
    double F;
    multiset< pair< double, pair<int,int> >, MagnitudeCompare2 >  voxels;

    //Matrix of derivatives for each organ and voxel (may be useful for algorithms)
	//How much increase/decrease F increasing the voxel in one unity.
	vector< vector<double> > D;

	vector< vector<double> >& get_Z(){
		return Z;
	};


};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
