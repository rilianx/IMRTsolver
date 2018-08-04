/*
 * main.cpp
 *
 *  Created on:
 *      Author: leslie
 */

#include <iostream>
#include <iterator>

#include "EvaluationFunction.h"
#include "Plan.h"
#include "Collimator.h"
#include "Volume.h"
#include "ApertureILS.h"
#include "IntensityILS.h"
#include "args.hxx"


using namespace std;
using namespace imrt;


vector<Volume> createVolumes (string organ_filename, Collimator& collimator){
  ifstream organ_file(organ_filename.c_str(), ios::in);
  vector<string> organ_files;
  vector<Volume> volumes;
  string line;
  
  if (! organ_file) 
    cerr << "ERROR: unable to open instance file: " << organ_filename << ", stopping! \n";
  
  cout << "##Reading volume files." << endl;
  while (organ_file) {
    getline(organ_file, line);
    if (line.empty()) continue;
    cout << "##  " << line << endl;
    //Assuming one data point
    organ_files.push_back(line);
  }
  organ_file.close();
  cout << "##  Read " << organ_files.size() << " files"<< endl;
  
  for (int i=0; i<organ_files.size(); i++) 
    volumes.push_back(Volume(collimator, organ_files[i]));
  
  return(volumes);
}


int main(int argc, char** argv){

  int vsize=50;
  int bsize=20;
  double int0=4.0;
  int maxiter=5000;
  int maxtime=0;
  int max_apertures=5;
  double alpha=1.0;
  double beta=1.0;
  double maxdelta=5.0;
  double maxratio=3.0;
  bool search_aperture=false;
  bool search_intensity=false;
  string strategy="dao_ls";

  int initial_intensity=2;
  int max_intensity=28;
  int step_intensity=2;
  bool acceptance=ILS::ACCEPT_NONE;

  // ls params
  double prob_intensity=0.2;
  double temperature, initial_temperature=10;
  double min_temperature=0;
  double alphaT=0.95;
  int perturbation=2;

  int seed=time(NULL);

  // Tabu list <<beam,station*>, open> for intensity
  vector<pair<pair<int,Station*>,  bool > > tabu_list_inten;
  // Tabu list <<beam,station*>, <aperture,open>> for aperture
  vector<pair<pair<int,Station*>, pair<int, bool> > > tabu_list_aper;
  int tabusize=10;


	args::ArgumentParser parser("********* IMRT-Solver (Aperture solver) *********", "Example.\n../AS -s ibo_ls --maxiter=400 --maxdelta=8 --maxratio=6 --alpha=0.999 --beta=0.999 --bsize=5 --vsize=20 --max_ap=4 --seed=0 --int0=1 --open-setup");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	//args::ValueFlag<string> _format(parser, "string", "Format: (BR, BRw, 1C)", {'f'});
  args::ValueFlag<string> _strategy(parser, "string", "Strategy  (dao_ls|ibo_ls)", {'s', "strategy"});
	args::ValueFlag<int> _bsize(parser, "int", "Number of considered beamlets for selection ("+to_string(bsize)+")", {"bsize"});
	args::ValueFlag<int> _vsize(parser, "int", "Number of considered worst voxels ("+to_string(vsize)+")", {"vsize"});
  args::ValueFlag<int> _int0(parser, "int", "Initial intensity for beams  ("+to_string(int0)+")", {"int0"});
  args::ValueFlag<int> _max_apertures(parser, "int", "Initial intensity for the station  ("+to_string(max_apertures)+")", {"max_ap"});
  args::ValueFlag<int> _maxdelta(parser, "int", "Max delta  ("+to_string(maxdelta)+")", {"maxdelta"});
  args::ValueFlag<int> _maxratio(parser, "int", "Max ratio  ("+to_string(maxratio)+")", {"maxratio"});
  args::ValueFlag<double> _alpha(parser, "double", "Initial temperature for intensities  ("+to_string(alpha)+")", {"alpha"});
  args::ValueFlag<double> _beta(parser, "double", "Initial temperature for ratio  ("+to_string(beta)+")", {"beta"});
  args::ValueFlag<int> _maxiter(parser, "int", "Number of iterations ("+to_string(maxiter)+")", {"maxiter"});
  args::ValueFlag<int> _maxtime(parser, "int", "Maximum time in seconds ("+to_string(maxtime)+")", {"maxtime"});
  args::ValueFlag<int> _seed(parser, "int", "Seed  ("+to_string(seed)+")", {"seed"});

  args::Group dao_ls (parser, "Direct aperture local search:", args::Group::Validators::DontCare);
  args::Flag open_setup(parser, "open_setup", "Initialize apertures as open", {"open-setup"});
  
  args::ValueFlag<int> _initial_intensity(parser, "int", "Initial value aperture intensity  ("+to_string(initial_intensity)+")", {"initial-intensity"});
  args::ValueFlag<int> _max_intensity(parser, "int", "Max value aperture intensity  ("+to_string(max_intensity)+")", {"max-intensity"});
  args::ValueFlag<int> _step_intensity(parser, "int", "Step size for aperture intensity  ("+to_string(step_intensity)+")", {"step-intensity"});
  
  args::Flag ls_aperture(dao_ls, "ls_apertures", "Apply aperture local search", {"ls-aperture"});
  args::Flag ls_intensity(dao_ls, "ls_intensity", "Apply intensity local search", {"ls-intensity"});
  args::ValueFlag<double> _prob_intensity(dao_ls, "double", "Probability to search over intensity  ("+to_string(prob_intensity)+")", {"prob-intensity"});

  args::Group accept (parser, "Acceptance criterion:", args::Group::Validators::AtMostOne);
  args::Flag accept_best(accept, "accept-best", "Accept only improvement", {"accept-best"});
  args::Flag accept_sa(accept, "accept-sa", "Accept as simulated annealing", {"accept-sa"});

  args::ValueFlag<double> _temperature(parser, "double", "Temperature for acceptance criterion  ("+to_string(temperature)+")", {"temperature"});
  args::ValueFlag<double> _alphaT(parser, "double", "Reduction rate of the temperature  ("+to_string(alphaT)+")", {"alphaT"});
  args::ValueFlag<int> _perturbation(parser, "int", "Perturbation size  ("+to_string(perturbation)+")", {"perturbation-size"});
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

  if(_strategy) strategy=_strategy.Get();
  if(_bsize) bsize=_bsize.Get();
  if(_vsize) vsize=_vsize.Get();
  if(_maxdelta) maxdelta=_maxdelta.Get();
  if(_maxratio) maxratio=_maxratio.Get();
  if(_alpha) alpha=_alpha.Get();
  if(_beta) beta=_beta.Get();
  if(_int0) int0=_int0.Get();
  if(_maxiter) maxiter=_maxiter.Get();
  if(_maxtime) maxtime=_maxtime.Get();
  if(_max_apertures) max_apertures=_max_apertures.Get();
  if(_seed) seed=_seed.Get();
  if(_initial_intensity) initial_intensity=_initial_intensity.Get();
  if(_max_intensity) max_intensity=_max_intensity.Get();
  if(_step_intensity) step_intensity=_step_intensity.Get();
  if(_prob_intensity) prob_intensity=_prob_intensity.Get();
  if(_temperature) temperature=initial_temperature=_temperature.Get();
  if(_alphaT) alphaT=_alphaT.Get();
  if (ls_aperture) search_aperture=true;
  if (ls_intensity) search_intensity=true;
  if(!ls_aperture && !ls_intensity){
    search_aperture=true;
    search_intensity=true;
  }
  if (accept_sa) acceptance=ILS::ACCEPT_SA;
  if (accept_best) acceptance=ILS::ACCEPT_NONE;
  if (_perturbation) perturbation=_perturbation.Get();

  cout << "##**************************************************************************"<< endl;
  cout << "##**************************************************************************"<< endl;
  if(strategy=="dao_ls")
    cout << "##******** IMRT-Solver (Direct Aperture Optimization Local Search) *********"<< endl;
  else if(strategy=="ibo_ls")
    cout << "##******** IMRT-Solver (Intensity-based Optimization Local Search) *********"<< endl;
  cout << "##**************************************************************************"<< endl;
  cout << "##**************************************************************************"<< endl;
  
  vector<double> w={1,1,1};
  vector<double> Zmin={0,0,76};
  vector<double> Zmax={65,60,1000};
  
  Collimator collimator("data/test_instance_coordinates.txt");
  vector<Volume> volumes= createVolumes ("data/test_instance_organs.txt", collimator);
  
  Plan P(w, Zmin, Zmax, collimator, volumes, max_apertures, max_intensity, initial_intensity, open_setup);
  double best_eval=P.getEvaluation();
  
  
  /*vector<Station*> stations(5);
  Station* station;
  for(int i=0;i<5;i++){
	  station = new Station(collimator,volumes, i*70, max_apertures, max_intensity, initial_intensity, open_setup);
    //station = new Station(collimator,volumes, i*70, max_apertures);
	  station->generateIntensity();
	  stations[i]=station;
  }*/
  
  
  cout << "##" << endl << "##**************************************************************************"<< endl;
  cout << "##*********************************** INFO *********************************"<< endl;
  cout << "##**************************************************************************"<< endl;
  cout << "##" << endl << "## Solver: "<< endl;
  cout << "##   Iterations: " << maxiter << endl;
  cout << "##   Time: " << maxtime << endl;
  cout << "##   Seed: " << seed << endl;
  cout << "##   Temperature: " << temperature << endl;
  cout << "##   alpha: " << alpha << endl;
  if (search_aperture)
    cout << "##   Searching: aperture pattern" << endl;
  if (search_intensity)
    cout << "##   Searching: intensity" << endl;
  cout << "##   Probability intensity ls: " << prob_intensity << endl;
  cout << "##   Open initial setup: " << open_setup << endl;
  cout << "##   Initial intensity: " << initial_intensity << endl;
  cout << "##   Max intensity: " << max_intensity << endl;
  cout << "##   Step intensity: " << step_intensity << endl;
  cout << "##   Perturbation size: " << perturbation << endl;
  cout << "##" << endl << "## Colimator configuration: "<< endl;
  cout << "##   Stations: " << collimator.getNbAngles() << endl;
  cout << "##   Angles: ";
  for (int i=0; i<collimator.getNbAngles();i++) cout << collimator.getAngle(i) << " ";
  cout << endl;
  cout << "##   Apertures: " << max_apertures << endl;
  
 


  cout << "##" << endl << "## Instance information: "<< endl;
  cout << "##   Volumes: " << volumes.size() << endl;

  cout << "##" << endl << "##**************************************************************************"<< endl;
  cout << "##********************************** SEARCH ********************************"<< endl;
  cout << "##**************************************************************************"<< endl;


 
  cout << "## Initial solution: " << best_eval << endl;
  cout  << "##" << endl;
  ILS* ils;
  if(strategy=="dao_ls")    
    ils = new ApertureILS(bsize, vsize, search_intensity, search_aperture, prob_intensity, step_intensity, initial_temperature, alphaT, perturbation, acceptance);
  else if(strategy=="ibo_ls")
    ils = new IntensityILS(bsize, vsize, maxdelta, maxratio, alpha, beta);

  ils->search(P, maxtime, maxiter);

  cout << "##**************************************************************************"<< endl;
  cout << "##******************************* RESULTS **********************************"<< endl;
  cout << "##**************************************************************************"<< endl;
  cout << "##"<<endl;
  cout << "## Best solution found: " <<  P.getEvaluation() << endl; 
/*
	//cout << endl;
	//for(int i=0;i<5;i++){
	//	//stations[i]->printIntensity();
	//	stations[i]->printIntensity(false);
  //      //cout << "nb_apertures:" << stations[i]->int2nb.size() << endl;
  //  }
	//cout << endl;


	cout << "********   Summary of the results    *********"<< endl;
	best_eval = F.eval(P,w,Zmin,Zmax);
	cout << "Final solution: " << best_eval << endl << endl;

  F.generate_voxel_dose_functions ();
  system("python plotter/plot.py");*/

	return 0;

}
