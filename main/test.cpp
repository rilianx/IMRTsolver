/*
 * main.cpp
 *
 *  Created on: 3 may. 2018
 *      Author: iaraya
 */

#include <iostream>
#include "EvaluationFunction.h"

using namespace std;
using namespace imrt;

int main(){

	vector<double> w={-1,1,1};
	vector<double> Zmin={70,0,0};
	vector<double> Zmax={90,20,20};

	EvaluationFunction F(3, 331, w, Zmin, Zmax);

	F.set_deposition_matrix(0,0,"data/65429427DDM_LECHOPROST.dat");
	F.set_deposition_matrix(0,1,"data/65429427DDM_RECTO.dat");
	F.set_deposition_matrix(0,2,"data/65429427DDM_VEJIGA.dat");
	return 0;
}
