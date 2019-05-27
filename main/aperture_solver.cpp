/*
 * main.cpp
 *
 *  Created on:
 *      Author: leslie
 */

#include <iostream>
#include <iterator>
#include <set>
#include <stack>
#include <unistd.h>
#include <string.h>

#include "IntensityGenerator.h"
#include "EvaluationFunction.h"
#include "Plan.h"
#include "Collimator.h"
#include "Volume.h"
#include "ApertureILS.h"
#include "IntensityILS.h"
#include "args.hxx"


using namespace std;
using namespace imrt;

set<int> get_angles(string file, int n){
  ifstream _file(file.c_str(), ios::in);
  string line;

  if (! _file)
    cerr << "ERROR: unable to open instance file: " << file << ", stopping! \n";

  //cout << "##Reading volume files." << endl;
  getline(_file, line);
  _file.close();

  set<int> angles;
  stack<string> q;
  char delim=' ';
  std::size_t current, previous = 0;
  current = line.find(delim);
  //cout << line << endl;
  while (current != std::string::npos) {
    angles.insert(atoi(line.substr(previous, current - previous).c_str()));
    previous = current + 1;
    current = line.find(delim, previous);
  }
  angles.insert(atoi(line.substr(previous, current - previous).c_str()));

  return angles;
}


vector<Volume> createVolumes (string organ_filename, Collimator& collimator){
  ifstream organ_file(organ_filename.c_str(), ios::in);
  vector<string> organ_files;
  vector<Volume> volumes;
  string line;

  if (! organ_file)
    cerr << "ERROR: unable to open instance file: " << organ_filename << ", stopping! \n";

  cout << "##Reading volume files." << endl;
  getline(organ_file, line);
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

  // Target heuristic
  int vsize = 50;
  int bsize = 20;
  double maxdelta = 5.0;
  double maxratio = 3.0;

  // Budget and execution variables
  int maxiter = 5000;
  int maxtime = 0;
  int seed = time(NULL);

  // Intensity local search parameters
  double alpha = 1.0;
  double beta = 1.0;

  // Aperture and intensity configuration
  int max_apertures = 5;
  int open_apertures = -1;
  int initial_setup;
  int initial_intensity = 2;
  int max_intensity = 28;
  int step_intensity = 2;

  // Type of local search
  string strategy = "dao_ls";
  int ls_type = ApertureILS::FIRST_IMPROVEMENT;
  bool acceptance = ILS::ACCEPT_NONE;
  bool targeted_search = false;

  // Aperture representation local search
  bool search_aperture = false;
  bool search_intensity = false;
  double prob_intensity = 0.2;
  int perturbation = 2;
  double temperature;
  double initial_temperature = 10;
  double min_temperature = 0;
  double alphaT = 0.95;
  
  //Intensities Generator parameter
  double alpha2 = 0.5;

  string file="data/testinstance_0_70_140_210_280.txt";
  string file2="data/test_instance_coordinates.txt";
  char* file3=NULL;

  string path=".";

  // Tabu list <<beam,station*>, open> for intensity
  vector<pair<pair<int,Station*>,  bool > > tabu_list_inten;
  // Tabu list <<beam,station*>, <aperture,open>> for aperture
  vector<pair<pair<int,Station*>, pair<int, bool> > > tabu_list_aper;
  int tabusize=10;


	args::ArgumentParser parser("********* IMRT-Solver (Aperture solver) *********",
                             "Example.\n./AS -s ibo_ls --maxiter=400 --maxdelta=8 --maxratio=6 --alpha=0.999 --beta=0.999 --bsize=5 --vsize=20 --max-apertures=4 --seed=0 --open-apertures=1 --initial-intensity=4 --step-intensity=1 --file-dep=data/Equidistantes/equidist00.txt --file-coord=data/Equidistantes/equidist-coord.txt");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	//args::ValueFlag<string> _format(parser, "string", "Format: (BR, BRw, 1C)", {'f'});
  args::ValueFlag<string> _strategy (parser, "string", "Strategy  (dao_ls|ibo_ls)", {'s', "strategy"});

  args::ValueFlag<int> _bsize(parser, "int", "Number of considered beamlets for selection ("+to_string(bsize)+")", {"bsize"});
	args::ValueFlag<int> _vsize(parser, "int", "Number of considered worst voxels ("+to_string(vsize)+")", {"vsize"});
 // args::ValueFlag<int> _int0(parser, "int", "Initial intensity for beams  ("+to_string(int0)+")", {"int0"});
  args::ValueFlag<int> _maxdelta(parser, "int", "Max delta  ("+to_string(maxdelta)+")", {"maxdelta"});
  args::ValueFlag<int> _maxratio(parser, "int", "Max ratio  ("+to_string(maxratio)+")", {"maxratio"});

  args::ValueFlag<int> _maxiter(parser, "int", "Number of iterations ("+to_string(maxiter)+")", {"maxiter"});
  args::ValueFlag<int> _maxtime(parser, "int", "Maximum time in seconds ("+to_string(maxtime)+")", {"maxtime"});
  args::ValueFlag<int> _seed(parser, "int", "Seed  ("+to_string(seed)+")", {"seed"});

  args::ValueFlag<double> _alpha(parser, "double", "Initial temperature for intensities  ("+to_string(alpha)+")", {"alpha"});
  args::ValueFlag<double> _beta(parser, "double", "Initial temperature for ratio  ("+to_string(beta)+")", {"beta"});

  args::Group setup (parser, "Initial solution setup (these override all provided configurations):", args::Group::Validators::DontCare);
  args::Flag open_max(setup, "open_max", "Open aperture setup with max intensity", {"open-max-setup"});
  args::Flag open_min(setup, "open_min", "Open aperture setup with min intensity", {"open-min-setup"});
  args::Flag closed_min(setup, "closed_min", "Closed aperture setup with min intensity", {"closed-min-setup"});
  args::Flag closed_max(setup, "closed_max", "Closed aperture setup with max intensity", {"closed-max-setup"});
  args::Flag all_rand(setup, "all_rand", "Random aperture setup with random intensity", {"rand-setup"});

  args::ValueFlag<int> _max_apertures(parser, "int", "Initial intensity for the station  ("+to_string(max_apertures)+")", {"max-apertures"});
  args::ValueFlag<int> _open_apertures(parser, "int", "Number of initialized open apertures (-1: all, default:"+to_string(open_apertures)+")", {"open-apertures"});
  args::ValueFlag<int> _initial_intensity(parser, "int", "Initial value aperture intensity  ("+to_string(initial_intensity)+")", {"initial-intensity"});
  args::ValueFlag<int> _max_intensity(parser, "int", "Max value aperture intensity  ("+to_string(max_intensity)+")", {"max-intensity"});
  args::ValueFlag<int> _step_intensity(parser, "int", "Step size for aperture intensity  ("+to_string(step_intensity)+")", {"step-intensity"});

  args::Group dao_ls (parser, "Direct aperture local search:", args::Group::Validators::DontCare);
  args::Flag _targeted_search(dao_ls, "targeted_search", "Apply targeted local search", {"targeted"});
  args::Flag ls_aperture(dao_ls, "ls_apertures", "Apply aperture local search", {"ls-aperture"});
  args::Flag ls_intensity(dao_ls, "ls_intensity", "Apply intensity local search", {"ls-intensity"});
  args::ValueFlag<double> _prob_intensity(dao_ls, "double", "Probability to search over intensity  ("+to_string(prob_intensity)+")", {"prob-intensity"});
  args::Flag best_improv(dao_ls, "best_improv", "Apply best improvement when searching intensity or apertures", {"best-improvement"});
  args::Flag first_improv(dao_ls, "first_improv", "Apply first improvement when searching intensity or apertures", {"first-improvement"});
  args::Flag do_perturbate(dao_ls, "do_perturbate", "Perturbate after a selected criterion is triggered", {"perturbate"});
  args::ValueFlag<int> _perturbation(parser, "int", "Perturbation size  ("+to_string(perturbation)+")", {"perturbation-size"});

  args::Group accept (parser, "Acceptance criterion:", args::Group::Validators::AtMostOne);
  args::Flag accept_best(accept, "accept-best", "Accept only improvement", {"accept-best"});
  args::Flag accept_sa(accept, "accept-sa", "Accept as simulated annealing", {"accept-sa"});

  args::ValueFlag<double> _temperature(parser, "double", "Temperature for acceptance criterion  ("+to_string(temperature)+")", {"temperature"});
  args::ValueFlag<double> _alphaT(parser, "double", "Reduction rate of the temperature  ("+to_string(alphaT)+")", {"alphaT"});

	//args::Flag trace(parser, "trace", "Trace", {"trace"});
  args::ValueFlag<string> _file(parser, "string", "File with the deposition matrix", {"file-dep"});
  args::ValueFlag<string> _file2(parser, "string", "File with the beam coordinates", {"file-coord"});
  args::ValueFlag<string> _file3(parser, "string", "File with initial intensities", {"file-sol"});
  args::ValueFlag<string> _path(parser, "string", "Absolute path of the executable (if it is executed from other directory)", {"path"});
  args::Flag _plot(parser, "bool", "Generate plot and save in file", {"plot"});

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

	if(_maxiter) maxiter=_maxiter.Get();
	if(_maxtime) maxtime=_maxtime.Get();
	if(_seed) seed=_seed.Get();

	if(_alpha) alpha=_alpha.Get();
	if(_beta) beta=_beta.Get();

  if(_max_apertures) max_apertures=_max_apertures.Get();
  if(_open_apertures) open_apertures=_open_apertures.Get();

  if(_initial_intensity) initial_intensity=_initial_intensity.Get();
  if(_max_intensity) max_intensity=_max_intensity.Get();
  if(_step_intensity) step_intensity=_step_intensity.Get();
  if(_prob_intensity) prob_intensity=_prob_intensity.Get();

  if (_targeted_search) targeted_search = true;
  if (ls_aperture) search_aperture = true;
  if (ls_intensity) search_intensity = true;
  if(!ls_aperture && !ls_intensity){
    search_aperture=true;
    search_intensity=true;
  }
  if (_perturbation) perturbation=_perturbation.Get();

  if (first_improv) ls_type=ApertureILS::FIRST_IMPROVEMENT;
  if (best_improv) ls_type=ApertureILS::BEST_IMPROVEMENT;
  if (accept_sa) acceptance=ILS::ACCEPT_SA;
  if (accept_best) acceptance=ILS::ACCEPT_NONE;

  if(_temperature) temperature=initial_temperature=_temperature.Get();
  if(_alphaT) alphaT=_alphaT.Get();

  //This should be at the end given that will modify values
  if(open_max) {
    initial_setup=Station::OPEN_MAX_SETUP;
    initial_intensity=max_intensity;
  } else if (open_min) {
    initial_setup=Station::OPEN_MIN_SETUP;
    initial_intensity=0;
  } else if (closed_min) {
    initial_setup=Station::CLOSED_MIN_SETUP;
    initial_intensity=0;
  } else if (closed_max) {
    initial_setup=Station::CLOSED_MAX_SETUP;
    initial_intensity=max_intensity;
  } else if (all_rand) {
    if(strategy=="dao_ls")
      initial_setup=Station::RAND_RAND_SETUP;
    else if (strategy=="ibo_ls")
      initial_setup=Station::RAND_INTENSITIES;
  } else {
    initial_setup=Station::MANUAL_SETUP;
  }

  srand(seed);
  if (_file) file=_file.Get();
  if (_file2) file2=_file2.Get();
  if (_file3) file3=strdup(_file3.Get().c_str()); //intensidades de partida
  if (_path) path=_path.Get();

  chdir(path.c_str());

  cout << "##**************************************************************************"<< endl;
  cout << "##**************************************************************************"<< endl;
  if(strategy=="dao_ls") {
    cout << "##******** IMRT-Solver (Direct Aperture Optimization Local Search) *********"<< endl;
  } else if(strategy=="ibo_ls") {
    cout << "##******** IMRT-Solver (Intensity-based Optimization Local Search) *********"<< endl;
  }
  cout << "##**************************************************************************"<< endl;
  cout << "##**************************************************************************"<< endl;

  vector<double> w={1,1,1};
  vector<double> Zmin={0,0,76};
  vector<double> Zmax={65,60,1000};

  Collimator collimator(file2, get_angles(file, 5));
  vector<Volume> volumes= createVolumes (file, collimator);

  Plan P(w, Zmin, Zmax, collimator, volumes, max_apertures, max_intensity, initial_intensity, step_intensity, open_apertures, initial_setup,file3);

  double best_eval=P.getEvaluation();

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
  if (strategy=="dao_ls" && ls_type==ApertureILS::FIRST_IMPROVEMENT)
    cout << "##   Local search: first improvement"  << endl;
  if (strategy=="dao_ls" && ls_type==ApertureILS::BEST_IMPROVEMENT)
    cout << "##   Local search: best improvement"  << endl;
  if (strategy=="dao_ls" && targeted_search)
    cout << "##   Targeted search: yes"  << endl;
  if (strategy=="dao_ls" && !targeted_search)
    cout << "##   Targeted search: no"  << endl;

  cout << "##   Open initial setup: " ;
  if (initial_setup==Station::OPEN_MAX_SETUP) cout << "open max intensity" << endl;
  else if (initial_setup==Station::OPEN_MIN_SETUP) cout << "open min intensity" << endl;
  else if (initial_setup==Station::CLOSED_MAX_SETUP) cout << "closed max intensity" << endl;
  else if (initial_setup==Station::CLOSED_MIN_SETUP) cout << "closed min intensity" << endl;
  else if (initial_setup==Station::RAND_RAND_SETUP) cout << "random" << endl;
  else if (initial_setup==Station::RAND_INTENSITIES) cout << "random" << endl;
  else cout << "manual, " << open_apertures << " open apertures" << endl;
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
  
  cout << endl;
  for(int i=0;i<5;i++)
    P.printIntensity(i);
  cout << endl;

  ILS* ils;
  if(strategy=="dao_ls"){
    ils = new ApertureILS(bsize, vsize, search_intensity, search_aperture,
                          prob_intensity, step_intensity, initial_temperature,
                          alphaT, do_perturbate, perturbation, acceptance, ls_type);
   /* if (!targeted_search)*/
      ils-> notTargetedSearch(P, maxtime, maxiter);
    /*else
      ils-> beamTargetedSearch(P, maxtime, maxiter);*/
  }else if(strategy=="ibo_ls"){
/*
	    for(auto s:P.get_stations()){
	      for (int i=0; i<collimator.getXdim();i++) {
				  for (int j=0; j<collimator.getYdim(); j++) {
	          if(s->I(i,j)!=-1){
	        	  pair<double, double> vc=EvaluationFunction::getInstance().get_value_cost(
	        	              		s->getAngle(), s->pos2beam[make_pair(i,j)], Zmin, Zmax);
	            printf("%5.1f ",
	            vc.first);
	          }else cout << "   -1 " ;
	        }
	        cout << endl;
	      }
	      cout << endl;
	    }

	    for(auto s:P.get_stations()){
	      for (int i=0; i<collimator.getXdim();i++) {
				  for (int j=0; j<collimator.getYdim(); j++) {
	          if(s->I(i,j)!=-1){
	        	  pair<double, double> vc=EvaluationFunction::getInstance().get_value_cost(
	        	              		s->getAngle(), s->pos2beam[make_pair(i,j)], Zmin, Zmax);
	            printf("%6.1f ",
	            1/vc.second);
	          }else cout << "    -1 " ;
	        }
	        cout << endl;
	      }
	      cout << endl;
	    }
*/

    ils = new IntensityILS(step_intensity, bsize, vsize, maxdelta, maxratio, alpha, beta, perturbation);
    ils->beamTargetedSearch(P, maxtime, maxiter);


    
    cout << endl;
  	for(int i=0;i<5;i++)
  		P.printIntensity(i);
  	cout << endl;



    for(int i=0;i<50;i++) cout << ils->iLocalSearch(P, false) << endl;
    cout << P.eval() << endl  ;
  }else if(strategy=="intgen"){
    cout << "entra al if" << endl  ;
	  IntensityGenerator intgen;
    cout <<"Crea intgen" << endl  ;
	  intgen.generate(P,alpha2);
    exit(0);
  }




  cout << "##**************************************************************************"<< endl;
  cout << "##******************************* RESULTS **********************************"<< endl;
  cout << "##**************************************************************************"<< endl;



  cout << "##"<<endl;
  cout << "## Best solution found: " <<  P.getEvaluation() << endl;

  cout <<  P.getEvaluation() << " ";

  const list<Station*> stations=P.get_stations();

  int tot_alpha=0;
  for(auto s:stations){
    int alpha=s->get_sum_alpha(strategy);
    cout << alpha << " " ;
    tot_alpha+=alpha;
  }
  cout << tot_alpha << " ";

  int nb_apertures=0;
  for(auto s:stations){
    int ap=s->get_nb_apertures(strategy);
    cout << ap << " " ;
    nb_apertures+=ap;
  }
  cout << nb_apertures << endl;


	cout << endl;
	for(int i=0;i<5;i++)
		P.printIntensity(i);

	cout << endl;

	/*cout << "********   Summary of the results    *********"<< endl;
	best_eval = F.eval(P,w,Zmin,Zmax);
	cout << "Final solution: " << best_eval << endl << endl;

  F.generate_voxel_dose_functions ();*/

  set<int> l = get_angles(file, 5);
  if(_plot){
	  std::stringstream ss;
	  ss << "python plotter/plot.py " << *l.begin()/5 << "_" << strategy << "_" << initial_intensity << "_" << initial_setup << "_" << maxtime << "_" << seed <<".pdf";
	  std::string s = ss.str();
	  system(s.c_str());
  }

	return 0;

}
