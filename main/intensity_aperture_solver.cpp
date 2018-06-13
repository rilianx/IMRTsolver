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

  int vsize=20;
  int bsize=5;
  double int0=4.0;
  int maxiter=100;
  int max_apertures=4;

	args::ArgumentParser parser("********* IMRT-Solver (Intensity-aperture solver) *********", "An IMRT Solver.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	//args::ValueFlag<string> _format(parser, "string", "Format: (BR, BRw, 1C)", {'f'});
	args::ValueFlag<int> _bsize(parser, "int", "Number of considered beamlets for selection ("+to_string(bsize)+")", {"bsize"});
	args::ValueFlag<int> _vsize(parser, "int", "Number of considered worst voxels ("+to_string(vsize)+")", {"vsize"});
  args::ValueFlag<int> _int0(parser, "int", "Initial intensity for beams  ("+to_string(int0)+")", {"int0"});
  args::ValueFlag<int> _max_apertures(parser, "int", "Initial intensity for the station  ("+to_string(max_apertures)+")", {"max_ap"});
  args::ValueFlag<int> _maxiter(parser, "int", "Number of iterations ("+to_string(maxiter)+")", {"max_iter"});
	//args::Flag trace(parser, "trace", "Trace", {"trace"});
	//args::Positional<std::string> _file(parser, "instance", "Instance");

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

  if(_bsize) bsize=_bsize.Get();
  if(_vsize) vsize=_vsize.Get();
  if(_int0) int0=_int0.Get();
  if(_maxiter) maxiter=_maxiter.Get();
  if(_max_apertures) max_apertures=_max_apertures.Get();

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
	   station = new Station(collimator,volumes, i*70, max_apertures);
     station->setUniformIntensity(int0);
	   station->printIntensity();
	   stations[i]=station;
   }

   EvaluationFunction F(volumes);
   Plan P(F);
   for(int i=0;i<5;i++)
	   P.add_station(*stations[i]);


	vector<double> w={1,1,1};
	vector<double> Zmin={0,0,76};
	vector<double> Zmax={65,60,1000};

	double best_eval=F.eval(P,w,Zmin,Zmax);
	cout << "ev:" << best_eval << endl;

	for(int i=0;i<maxiter;i++){
		auto sb=F.best_beamlets(P, bsize, vsize);

		auto it=sb.begin();
		std::advance(it,rand()%sb.size());

		Station*s = it->second.first; int beamlet=it->second.second;
		bool sign=it->first.second; //impact in F (+ or -)

		double delta_intensity= rand()%3+1;
		if(sign) delta_intensity*=-1;

		double ratio= rand()% int(3) ;

		auto diff=s->increaseIntensity_repair(beamlet,delta_intensity,ratio);
		double eval=F.incremental_eval(*s,w,Zmin,Zmax, diff);

		if(eval > best_eval){ //reject the move
			s->revert(diff);
			F.undo_last_eval(w,Zmin,Zmax);
		}else{ //accept the move
			cout << eval << endl;
			best_eval=eval;
		}
	}

	cout << endl;
	for(int i=0;i<5;i++){
		stations[i]->printIntensity(true);
        //cout << "nb_apertures:" << stations[i]->int2nb.size() << endl;
    }
	cout << endl;

	cout << "best_eval:" << best_eval << endl;

	best_eval=F.eval(P,w,Zmin,Zmax);
		cout << "ev:" << best_eval << endl;

  F.generate_voxel_dose_functions ();
  system("python plotter/plot.py");

	return 0;

}
