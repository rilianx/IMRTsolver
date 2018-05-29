/*
 * main.cpp
 *
 *  Created on: 3 may. 2018
 *      Author: iaraya
 */

#include <iostream>

#include "../src/plan/EvaluationFunction.h"
#include "../src/plan/Plan.h"
#include "Collimator.h"
#include "Volume.h"
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




	vector< pair<int, string> > coord_files(5);
  coord_files[0]=(make_pair(0,"data/CERR_Prostate/CoordinatesBeam_0.txt"));
  coord_files[1]=(make_pair(70,"data/CERR_Prostate/CoordinatesBeam_70.txt"));
  coord_files[2]=(make_pair(140,"data/CERR_Prostate/CoordinatesBeam_140.txt"));
  coord_files[3]=(make_pair(210,"data/CERR_Prostate/CoordinatesBeam_210.txt"));
  coord_files[4]=(make_pair(280,"data/CERR_Prostate/CoordinatesBeam_280.txt"));

  vector<string> organ_files;
  organ_files.push_back("data/CERR_Prostate/DAO_DDM_BLADDER.dat");
  organ_files.push_back("data/CERR_Prostate/DAO_DDM_PTVHD.dat");
  organ_files.push_back("data/CERR_Prostate/DAO_DDM_RECTUM.dat");

  	Collimator collimator(coord_files);
  	collimator.printAxisValues();
  	collimator.printActiveBeam();

	vector<Volume> volumes;

  for (int i=0; i<organ_files.size(); i++) {
	  volumes.push_back(Volume(collimator, organ_files[i]));
	}

   //volumes[0].print_deposition();

	 cout << volumes[0].getDepositionMatrix(0)(0,0) << endl;


   Station station(collimator,volumes, 0, 10);
   station.generateIntensity();
   station.printIntensity();

   Station station2(collimator,volumes, 70, 10);
   station2.generateIntensity();
   station2.printIntensity();

   EvaluationFunction F(volumes);

   Plan P(F);

   P.add_station(station);
   P.add_station(station2);

	vector<double> w={1,1,1};
	vector<double> Zmin={70,0,0};
	vector<double> Zmax={90,20,20};

	cout << "ev:" << P.eval(w,Zmin,Zmax) << endl;


	return 0;

}
