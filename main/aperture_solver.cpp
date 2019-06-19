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

  int seed = time(NULL);

  // Budget and execution variables
  int maxiter = 5000;
  int maxtime = 0;
  int maxeval = 0;

  // Intensity local search parameters
  double maxdelta = 5.0;
  double maxratio = 3.0;
  double alpha = 1.0;
  double beta = 1.0;

  // Aperture and intensity configuration
  int initial_setup;
  int max_apertures = 5;
  int initial_intensity = 2;
  int max_intensity = 28;
  int step_intensity = 1;

  // Type of local search
  string strategy = "dao_ls";
  LSType ls_type = LSType::first;
  NeighborhoodType neighborhood = NeighborhoodType::mixed;
  double prob_intensity = 0.2;

  // Target beamlet heuristic
  bool targeted_search = false;
  int vsize = 50;
  int bsize = 20;

  // Perturbation
  int perturbation = 2;

  // Files
  string path = ".";
  string file = "data/testinstance_0_70_140_210_280.txt";
  string file2 = "data/test_instance_coordinates.txt";
  char* file3 = NULL;

  

  args::ArgumentParser parser("********* IMRT-Solver (Aperture solver) *********",
                             "Example.\n./AS -s ibo_ls --maxiter=400 --maxdelta=8 --maxratio=6 --alpha=0.999 --beta=0.999 --bsize=5 --vsize=20 --max-apertures=4 --seed=0 --open-apertures=1 --initial-intensity=4 --step-intensity=1 --file-dep=data/Equidistantes/equidist00.txt --file-coord=data/Equidistantes/equidist-coord.txt");

  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<int>    _seed     (parser, "int", "Seed  (" + 
                                     to_string(seed)+")", {"seed"});

  // Search strategies
  args::Group strat (parser, "Strategy options:");
  args::ValueFlag<string> _strategy (strat, "string", 
                                     "Strategy  (dao_ls|ibo_ls)", 
                                     {'s', "strategy"});
  args::ValueFlag<string> _ls (strat, "string", 
                               "Local search strategy (best|first)", 
                               {'l', "ls"});

  // Execution parameters
  args::Group budget (parser, "Budget options:");
  args::ValueFlag<int>    _maxiter  (budget, "int", 
                                    "Number of iterations (" + 
                                     to_string(maxiter)+ ")", {"maxiter"});
  args::ValueFlag<int>    _maxtime  (budget, "int", 
                                    "Maximum time in seconds (" + 
                                     to_string(maxtime)+")", {"maxtime"});
  args::ValueFlag<int>    _maxeval  (budget, "int", 
                                    "Number of evaluations (" + 
                                     to_string(maxiter)+ ")", {"maxiter"});

  
  // Initial collimator setup (initial solution: aperture, intensity)
  args::Group isetup (parser, "Initial collimator setup:");
  args::ValueFlag<string> _setup (isetup, "string", 
                                  "Initial setup  (open_max|open_min|closed_max|closed_min|random|manual)", 
                                   {'t', "setup"});
  args::ValueFlag<int> _max_apertures     (isetup, "int", 
                                          "Number of apertures per angle (station) (" + 
                                           to_string(max_apertures)+")", {"max-apertures"});
  args::ValueFlag<int> _initial_intensity (isetup, "int", 
                                          "Initial value aperture intensity  (" + 
                                           to_string(initial_intensity) + ")", 
                                          {"initial-intensity"});

  args::Group intensity (parser, "Intensity options:");
  args::ValueFlag<int> _max_intensity     (intensity, "int", 
                                          "Max value aperture intensity  (" + 
                                           to_string(max_intensity)+")", 
                                          {"max-intensity"});
  args::ValueFlag<int> _step_intensity    (intensity, "int", 
                                          "Step size for aperture intensity  (" + 
                                           to_string(step_intensity)+")", 
                                          {"step-intensity"});


  // Neighborhood parameters
  args::Group neighborhoodsel (parser, "Neighborhood selection:", 
                        args::Group::Validators::AtMostOne);
  args::ValueFlag<string> _nsimple (neighborhoodsel , "string", 
                               "Simple neighborhood in local search  (aperture|intensity|mixed)", 
                               {"ls_simple"});
  args::ValueFlag<string> _nseq (neighborhoodsel , "string", 
                               "Sequential neighborhood in local search starting by (intensity|aperture)", 
                               {"ls_sequential"});
  args::ValueFlag<double> _nprob (neighborhoodsel , "double", 
                               "Probabilistic sequential neighborhood in local search [0,1]", 
                               {"ls_sequentialp"});

  // Target beamlet heuristic parameters
  args::Group heur (parser, "Heuristic options:");
  args::Flag _targeted_search (heur, "targeted_search", 
                                "Apply targeted local search", {"targeted"});
  args::ValueFlag<int>    _bsize    (heur, "int", 
                                    "Number of considered beamlets for selection (" + 
                                     to_string(bsize)+")", {"bsize"});
  args::ValueFlag<int>    _vsize    (heur, "int", 
                                    "Number of considered worst voxels (" + 
                                     to_string(vsize)+")", {"vsize"});
  
  // Perturbation parameters
  args::Group perargs (parser, "Heuristic options:");
  args::Flag do_perturbate                (perargs, "do_perturbate", 
                                          string("Perturbation is triggered  after a ") + 
                                          "selected criterion ", {"perturbate"});
  args::ValueFlag<int> _perturbation      (perargs, "int", 
                                          "Perturbation size  (" + 
                                          to_string(perturbation)+")", 
                                          {"perturbation-size"});

  // Intensity matrix representation ls parameters 
  // Aperture representation local search parameters
  args::Group ibo_ls (parser, "Intensity matrix local search:");
  args::ValueFlag<int>    _maxdelta (ibo_ls, "int", 
                                    "Max delta  (" + 
                                     to_string(maxdelta)+")", {"maxdelta"});
  args::ValueFlag<int>    _maxratio (ibo_ls, "int", 
                                    "Max ratio  (" + 
                                     to_string(maxratio)+")", {"maxratio"});
  args::ValueFlag<double> _alpha    (ibo_ls, "double", 
                                    "Initial temperature for intensities  (" + 
                                     to_string(alpha)+")", {"alpha"});
  args::ValueFlag<double> _beta     (ibo_ls, "double", 
                                    "Initial temperature for ratio  (" + 
                                     to_string(beta)+")", {"beta"});

  // Problem file parameters
  args::ValueFlag<string> _file  (parser, "string", 
                                 "File with the deposition matrix", {"file-dep"});
  args::ValueFlag<string> _file2 (parser, "string", 
                                 "File with the beam coordinates", {"file-coord"});
  args::ValueFlag<string> _file3 (parser, "string", 
                                 "File with initial intensities", {"file-sol"});
  args::ValueFlag<string> _path  (parser, "string", 
                                 string("Absolute path of the executable ") + 
                                 "(if it is executed from other directory)", {"path"});
  args::Flag _plot               (parser, "bool",  
                                 "Generate plot and save in file", {"plot"});

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

	if(_maxiter) maxiter=_maxiter.Get();
	if(_maxtime) maxtime=_maxtime.Get();
  if(_maxeval) maxeval=_maxeval.Get();
	if(_seed) seed=_seed.Get();

	if(_maxdelta) maxdelta=_maxdelta.Get();
	if(_maxratio) maxratio=_maxratio.Get();
	if(_alpha) alpha=_alpha.Get();
	if(_beta) beta=_beta.Get();

	if(_bsize) bsize=_bsize.Get();
	if(_vsize) vsize=_vsize.Get();

  if(_max_apertures) max_apertures=_max_apertures.Get();
  if(_initial_intensity) initial_intensity=_initial_intensity.Get();
  if(_max_intensity) max_intensity=_max_intensity.Get();
  if(_step_intensity) step_intensity=_step_intensity.Get();

  // Initial collimator setup
  if (_setup) {
    string setup = _setup.Get();
    if (setup == "open_max"){
      initial_setup = Station::OPEN_MAX_SETUP;
      initial_intensity = max_intensity;
    } else if (setup == "open_min"){
      initial_setup = Station::OPEN_MIN_SETUP;
      initial_intensity = 0;
    } else if (setup == "closed_max") {
      initial_setup = Station::CLOSED_MAX_SETUP;
      initial_intensity = max_intensity;
    } else if (setup == "closed_min") {
      initial_setup = Station::CLOSED_MIN_SETUP;
      initial_intensity = 0;
    } else if (setup == "random") {
      //TODO: Why these ones are different??
      if(strategy=="dao_ls")
        initial_setup = Station::RAND_RAND_SETUP;
      else if (strategy=="ibo_ls")
        initial_setup = Station::RAND_INTENSITIES;
    } else if (setup == "manual") {
      initial_setup = Station::MANUAL_SETUP;
    } else {
      cout << "Error: setup "<< setup << " not recognized"<< endl;
    }
  }
 
  // Local search strategy
  if (_ls) {
    string nn = _ls.Get();
    if (nn == "first") ls_type = LSType::first;
    else if (nn == "best") ls_type = LSType::best;
    else {
      cout << "Error: local search strategy " << nn << 
              " not recognized."<< endl;
    }
  }

  // Neighborhood
  if (_nsimple) {
    string nn = _nsimple.Get();
    if (nn == "intensity") 
      neighborhood = NeighborhoodType::intensity;
    else if (nn == "aperture") 
      neighborhood = NeighborhoodType::aperture;
    else if (nn == "mixed")
      neighborhood = NeighborhoodType::mixed;
    else 
      cout << "Neighborhood operator " << 
              neighborhood << "not recognized"<<endl;
  } else if (_nseq) {
    //TODO: implement aperture first!
    neighborhood = NeighborhoodType::sequential;
  } else if (_nprob) {
    prob_intensity = _nprob.Get();
    neighborhood = NeighborhoodType::sequentialp;
  }

  if (_targeted_search) targeted_search = true;
  if (_perturbation) perturbation=_perturbation.Get();

  // Archivos del problem
  if (_file) file=_file.Get();
  if (_file2) file2=_file2.Get();
  if (_file3) file3=strdup(_file3.Get().c_str()); //intensidades de partida
  if (_path) path=_path.Get();
  chdir(path.c_str());

  cout << "##**************************************************************************" 
       << endl;
  cout << "##**************************************************************************"
       << endl;
  if(strategy=="dao_ls") {
    cout << "##******** IMRT-Solver (Direct Aperture Optimization Local Search) *********"
         << endl;
  } else if(strategy=="ibo_ls") {
    cout << "##******** IMRT-Solver (Intensity-based Optimization Local Search) *********"
         << endl;
  }
  cout << "##**************************************************************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;

  // Iniciar generador de numeros aleatorios
  srand(seed);

  vector<double> w={1,1,1};
  vector<double> Zmin={0,0,76};
  vector<double> Zmax={65,60,1000};

  // Create colimator object and volumes
  Collimator collimator (file2, get_angles(file, 5));
  vector<Volume> volumes = createVolumes (file, collimator);
  
  // Create an initial plan
  Plan P (w, Zmin, Zmax, collimator, volumes, max_apertures, 
          max_intensity, initial_intensity, step_intensity, 
          initial_setup, file3);

  double best_eval = P.getEvaluation();

  cout << "##" << endl 
       << "##**************************************************************************"
       << endl;
  cout << "##*********************************** INFO *********************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;
  cout << "##" << endl << "## Solver: "<< endl;
  cout << "##   Iterations: " << maxiter << endl;
  cout << "##   Time: " << maxtime << endl;
  cout << "##   Evaluations: " << maxeval << endl;
  cout << "##   Seed: " << seed << endl;

  cout << "##   Open initial setup: " ;
  if (initial_setup==Station::OPEN_MAX_SETUP) cout << "open max intensity" << endl;
  else if (initial_setup==Station::OPEN_MIN_SETUP) cout << "open min intensity" << endl;
  else if (initial_setup==Station::CLOSED_MAX_SETUP) cout << "closed max intensity" << endl;
  else if (initial_setup==Station::CLOSED_MIN_SETUP) cout << "closed min intensity" << endl;
  else if (initial_setup==Station::RAND_RAND_SETUP) cout << "random" << endl;
  else if (initial_setup==Station::RAND_INTENSITIES) cout << "random" << endl;
  cout << "##   Apertures: " << max_apertures << endl;
  cout << "##   Initial intensity: " << initial_intensity << endl;
  cout << "##   Max intensity: " << max_intensity << endl;
  cout << "##   Step intensity: " << step_intensity << endl;

  if (ls_type == LSType::first)
    cout << "##   Local search: first improvement"  << endl;
  if (ls_type==LSType::best)
    cout << "##   Local search: best improvement"  << endl;

  if (neighborhood == NeighborhoodType::intensity)
    cout << "##   Neighborhood: intensity" << endl;
  else if (neighborhood == NeighborhoodType::aperture)
    cout << "##   Neighborhood: aperture" << endl;
  else if (neighborhood == NeighborhoodType::mixed)
    cout << "##   Neighborhood: mixed" << endl;
  else if (neighborhood == NeighborhoodType::sequential)
    cout << "##   Neighborhood: sequential" << endl;
  else if (neighborhood == NeighborhoodType::sequentialp) {
    cout << "##   Neighborhood: sequential probabilisty intensity: " << 
            prob_intensity << endl;
  }

  if (targeted_search)
    cout << "##   Targeted search: yes"  << endl;
  else
    cout << "##   Targeted search: no"  << endl;

  cout << "##   Perturbation size: " << perturbation << endl;

  cout << "##" << endl << "## Colimator configuration: "<< endl;
  cout << "##   Stations: " << collimator.getNbAngles() << endl;
  cout << "##   Angles: ";
  for (int i=0; i<collimator.getNbAngles();i++) 
    cout << collimator.getAngle(i) << " ";
  cout << endl;
 

  cout << "##" << endl << "## Instance information: "<< endl;
  cout << "##   Volumes: " << volumes.size() << endl;

  cout << "##" << endl 
       << "##**************************************************************************"
       << endl;
  cout << "##********************************** SEARCH ********************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;



  cout << "## Initial solution: " << best_eval << endl;
  cout  << "##" << endl;
  cout << endl;
  for(int i=0;i<5;i++)
    P.printIntensity(i);
  cout << endl;
  getchar();

  ILS* ils;
  if (strategy=="dao_ls") {
    ils = new ApertureILS(bsize, vsize, prob_intensity, step_intensity);
    ils->iteratedLocalSearch(P, maxtime, maxeval, ls_type, neighborhood, LSTarget::none);
  }else if(strategy=="ibo_ls"){

//    ils = new IntensityILS(step_intensity, bsize, vsize, maxdelta, maxratio, alpha, beta, perturbation);
//    ils->iteratedLocalSearch(P, maxtime, maxeval,LSType::first,NeighborhoodType::mixed,
      //                           LSTarget::none);
    
    cout << endl;
  	for(int i=0;i<5;i++)
  		P.printIntensity(i);
  	cout << endl;



    for(int i=0;i<50;i++) cout << ils->iLocalSearch(P, false) << endl;
    cout << P.eval() << endl  ;
  }else if(strategy=="intgen"){
	  IntensityGenerator intgen;
	  intgen.generate(P);
  }




  cout << "##**************************************************************************"
       << endl;
  cout << "##******************************* RESULTS **********************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;
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
