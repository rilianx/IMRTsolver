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

class EvalGscore : public EvaluationFunction{

    virtual double eval(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	virtual double incremental_eval(Station& station, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax,
			list< pair< int, double > >& diff);

    virtual double get_delta_eval(list< pair< int, double > >& diff, double angle,
                        const vector<double>& w, const vector<double>& Zmin, const vector<double>& Zmax, 
						int n_voxels=999999, map< pair< int, int >, double >* Z_diff=NULL) const;

}

}
