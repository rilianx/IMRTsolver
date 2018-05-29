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

#include "Plan.h"
#include "Volume.h"
#include "Matrix.h"

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

	virtual ~EvaluationFunction();

	// Generate the dose distribution matrices Z for each organ
	void generate_Z(const Plan& p);

	// Eval the cost F based on the dose deposition matrix Z
	double eval(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);




private:
	//dose distribution vectors for each organ
	vector< vector<double> > Z;

	//number of organs, including the tumor
	int nb_organs;

	//number of voxels for each organ
	vector<int> nb_voxels;


};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
