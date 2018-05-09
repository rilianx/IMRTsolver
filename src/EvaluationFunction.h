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
#include "Plan.h"

#ifndef EVALUATIONFUNCTION_H_
#define EVALUATIONFUNCTION_H_

using namespace std;
using namespace maths;


namespace imrt {

/* The Evaluation Function
 *
 * The evaluation function based on [this paper](https://drive.google.com/file/d/1YfMNk4GhBK97gSQ0nvpJAnyM6A3EPv61/view).
 * Basically it penalizes voxels with doses larger than Zmax (healthy organs), and with doses smaller than Zmin (tumor)
 */

class EvaluationFunction {
public:

	/*Constructor of the evaluator.
	 * @nb_organs 	number of organs (includes the tumor)
	 * @nb_beamlets number of beamleats in each station
	 * @nb_voxels 	number of voxels of organ i
	 * @w			penalization related to the organ i
	 * @Zmin		minimum acceptable dose for the organ i (0.0 if it is a healthy organ)
	 * @Zmax		maximum acceptable dose for the organ i
	 */

	EvaluationFunction(int nb_organs, int nb_beamlets, vector<int>& nb_voxels, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	virtual ~EvaluationFunction();

	/* Initialize the deposition matrix Z[id] from a file
	 *
	 * This matrices should be generated automatically from the original organs
	 * and the given angle
	 */
	bool set_deposition_matrix(double angle, int id_organ, string file);

	// Generate the dose distribution matrices Z for each organ
	void generate_Z(const Plan& p);

	// Eval the cost F based on the dose deposition matrix Z
	double eval(const Plan& p, bool genZ=true);

	//Incremental evaluation of F based on the modified beamlets of the station S
	double incremental_eval(Station& S);


private:

	/* Dose deposition matrices for each angle and organ
	 *
	 * D[a][o](k,b): Dose delivered to voxel k of the organ o by the beamlet b and angle a
	 */
	map<double, vector<Matrix>> D;


	//dose distribution vectors for each organ
	vector< vector<double> > Z;

	//penalizations for delivering less (w_min) or greater (w_max) than the dose required for each organ
	vector<double> w;

	// minimum and maximum dose required by each organ. Zmin=0 for the organs
	vector<double> Zmin;
	vector<double> Zmax;

	//last performed evaluation, used for incremental evaluations
	double last_eval;

	//number of organs, including the tumor
	int nb_organs;

	//number of beamlets of each station
	int nb_beamlets;

	//number of voxels for each organ
	vector<int> nb_voxels;
};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
