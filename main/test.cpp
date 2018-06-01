/*
 * main.cpp
 *
 *  Created on: 3 may. 2018
 *      Author: iaraya
 */

#include <iostream>
#include <iterator>

#include "EvaluationFunction.h"
#include "Plan.h"
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
	organ_files.push_back("data/CERR_Prostate/DAO_DDM_RECTUM.dat");
	organ_files.push_back("data/CERR_Prostate/DAO_DDM_PTVHD.dat");

  	Collimator collimator(coord_files);
  //	collimator.printAxisValues();
  //	collimator.printActiveBeam();

	vector<Volume> volumes;

  for (int i=0; i<organ_files.size(); i++) {
	  volumes.push_back(Volume(collimator, organ_files[i]));
	}

   //volumes[0].print_deposition();

   vector<Station*> stations(5);
   Station* station;
   for(int i=0;i<5;i++){
	   station = new Station(collimator,volumes, i*70, 4);
	   station->generateIntensity();
	   station->printIntensity();
	   stations[i]=station;
   }

   EvaluationFunction F(volumes);
   Plan P(F);
   for(int i=0;i<5;i++)
	   P.add_station(*stations[i]);


	vector<double> w={1,1,1};
	vector<double> Zmin={0,0,70};
	vector<double> Zmax={50,50,1000};

	double best_eval=F.eval(P,w,Zmin,Zmax);
	cout << "ev:" << best_eval << endl;

	for(int i=0;i<100;i++){
		auto sb=F.best_beamlets(P, 10, 10);
		auto it=sb.begin();
		std::advance(it,rand()%sb.size());

		Station*s = it->second.first; int beamlet=it->second.second;
		bool sign=it->first.second;

		//cout << eval_beamlet << endl;

		double intensity=rand()%3;
		if(sign) intensity*=-1;

		double ratio= rand()%5 ;
		//if(eval_beamlet<0) intensity=-intensity;


		auto diff=s->increaseIntensity(beamlet,intensity,ratio);
		double eval=F.incremental_eval(*s,w,Zmin,Zmax, diff);

		if(eval > best_eval){ //reject the move
			s->revert(diff);
			F.undo_last_eval(w,Zmin,Zmax);
		}else{ //accept the move
			cout << eval << endl;
			best_eval=eval;
		}
	}

	for(int i=0;i<5;i++)
		stations[i]->printIntensity();

	cout << "best_eval:" << best_eval << endl;

	best_eval=F.eval(P,w,Zmin,Zmax);
		cout << "ev:" << best_eval << endl;

  F.generate_voxel_dose_functions ();
  system("python plotter/plot.py");

	return 0;

}
