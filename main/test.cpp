/*
 * main.cpp
 *
 *  Created on: 3 may. 2018
 *      Author: iaraya
 */

#include <iostream>
#include "EvaluationFunction.h"
#include "Plan.h"
#include "args.hxx"


using namespace std;
using namespace imrt;

int main(int argc, char** argv){

	args::ArgumentParser parser("********* IMRT-Solver *********", "An IMRT Solver.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	//args::ValueFlag<string> _format(parser, "string", "Format: (BR, BRw, 1C)", {'f'});
	//args::ValueFlag<double> _min_fr(parser, "double", "Minimum volume occupied by a block (proportion)", {"min_fr"});
	//args::Flag trace(parser, "trace", "Trace", {"trace"});
	args::Positional<std::string> _file(parser, "instance", "Instance");

	try
	{
		parser.ParseCLI(argc, argv);

	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	string file=_file.Get();


	vector<double> w={1,1,1};
	vector<double> Zmin={70,0,0};
	vector<double> Zmax={90,20,20};
	vector<int> nb_voxels={13081,8500,19762};

	EvaluationFunction F(3, 330, nb_voxels, w, Zmin, Zmax);
	F.set_deposition_matrix(0,0,"data/65429427DDM_LECHOPROST.dat");
	F.set_deposition_matrix(0,1,"data/65429427DDM_RECTO.dat");
	F.set_deposition_matrix(0,2,"data/65429427DDM_VEJIGA.dat");

	Station S(33,10);
	Plan P;
	P.add_station(S);
	S.set_intensity(10);

	cout << "ev:" << F.eval(P) << endl;

	S.set_aperture(10,4,6);  //xxxx--xxxx
	cout << "incr_ev:" << F.incremental_eval(S) << endl;
	cout << "ev:" << F.eval(P) << endl;

	S.set_aperture(11,4,9);  //xxxx-----x
	cout << "incr_ev:" << F.incremental_eval(S) << endl;
	cout << "ev:" << F.eval(P) << endl;

	S.close_aperture();
	cout << "incr_ev:" << F.incremental_eval(S) << endl;
	cout << "ev:" << F.eval(P) << endl;

	return 0;
}
