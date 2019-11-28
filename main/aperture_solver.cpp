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
#include "IntensityILS2.h"
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
    cerr << "ERROR: unable to open instance file: " <<
             organ_filename << ", stopping! \n";

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
  int ibo_evals = 2500;
  int maxtime = 0;
  int maxeval = 0;

  // Intensity local search parameters
  double maxdelta = 5.0;
  double maxratio = 3.0;
  double alpha = 1.0;
  double beta = 1.0;

  // Aperture and intensity configuration
  StationSetup initial_setup = StationSetup::open_all_min;
  int max_apertures = 5;
  int initial_intensity = 2;
  int max_intensity = 28;
  int step_intensity = 1;

  // Type of local search
  string strategy = "dao_ls";
  LSType ls_type = LSType::first;
  NeighborhoodType neighborhood = NeighborhoodType::mixed;
  double prob_intensity = 0.2;
  int tabu_size = 0;

  // Target beamlet heuristic
  bool targeted_search = false;
  LSTargetType target_type = LSTargetType::target_none;
  int vsize = 50;
  int bsize = 20;

  // Perturbation
  PerturbationType perturbation_type = PerturbationType::p_mixed;
  int perturbation_size = 2;

  // Files
  string path = ".";
  string file = "data/testinstance_0_70_140_210_280.txt";
  string file2 = "data/test_instance_coordinates.txt";
  char* file3 = NULL;
  string convergence_file = "";
  string output_file = "";
  string json_file = "";



  args::ArgumentParser parser("********* IMRT-Solver (Aperture solver) *********",
                             "Example.\n../AS -s ibo_ls --maxeval=4000 --ls_sequential=intensity --setup=open_min --seed=2 --ls=first  --max-intensity=20");

  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<int>    _seed     (parser, "int", "Seed  (" +
                                     to_string(seed)+")", {"seed"});

  // Search strategies
  args::Group strat (parser, "Strategy options:");
  args::ValueFlag<string> _strategy (strat, "string",
                                     "Strategy  (dao_ls|ibo_ls|ibo+dao)",
                                     {'s', "strategy"});
  args::ValueFlag<string> _ls (strat, "string",
                               "Local search strategy (best|first)",
                               {'l', "ls"});
  args::ValueFlag<int>    _tabu_size (strat, "int",
                                    "Tabu list size(" +
                                     to_string(tabu_size)+")", {"tabu-size"});


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
                                     to_string(maxiter)+ ")", {"maxeval"});


  // Initial collimator setup (initial solution: aperture, intensity)
  args::Group isetup (parser, "Initial collimator setup:");
  args::ValueFlag<string> _setup (isetup, "string",
                                  "Initial setup  (open_max|open_min|closed_max|closed_min|random|manual|open_min_min|open_min_k). * open_min_k initializes intensity of open apertures  with the value of initial-intensity " + to_string(initial_intensity), {'t', "setup"});
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
                               "Simple neighborhood in local search  (aperture|intensity|mixed|imixed)",
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
  args::Group perargs (parser, "Perturbation:");
  args::ValueFlag<string> _perturbation_type (perargs , "string",
                                      "Type of perturbation to be applied (intensity|aperture|mixed)",
                                      {"perturbation"});
  args::ValueFlag<int> _perturbation_size (perargs, "int",
                                          "Perturbation size  (" +
                                          to_string(perturbation_size)+")",
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

  // ibo+dao parameter
  args::Group ibo_dao (parser, "IBO+DAO options:");
  args::ValueFlag<int>    _ibo_evals  (ibo_dao, "int",
                                    "Number of iterations ibo (" +
                                     to_string(maxiter)+ ")", {"ibo_evals"});

  // Problem file parameters
  args::Group io_opt (parser, "Input output options:");
  args::ValueFlag<string> _file  (io_opt, "string",
                                 "File with the deposition matrix", {"file-dep"});
  args::ValueFlag<string> _file2 (io_opt, "string",
                                 "File with the beam coordinates", {"file-coord"});
  args::ValueFlag<string> _file3 (io_opt, "string",
                                 "File with initial intensities", {"file-sol"});
  args::ValueFlag<string> _path  (io_opt, "string",
                                 string("Absolute path of the executable ") +
                                 "(if it is executed from other directory)", {"path"});
  args::Flag _plot               (io_opt, "bool",
                                 "Generate plot and save in file", {"plot"});

  // Output file parameters
  args::ValueFlag<string> _convergence_file (io_opt, "string",
                                 "File to output convergence", {"convergence"});

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

	if(_strategy) strategy = _strategy.Get();

	if(_maxiter) maxiter = _maxiter.Get();
	if(_maxtime) maxtime = _maxtime.Get();
  if(_maxeval) maxeval = _maxeval.Get();
	if(_seed) seed = _seed.Get();

	if(_maxdelta) maxdelta = _maxdelta.Get();
	if(_maxratio) maxratio = _maxratio.Get();
	if(_alpha) alpha = _alpha.Get();
	if(_beta) beta = _beta.Get();

	if(_bsize) bsize = _bsize.Get();
	if(_vsize) vsize = _vsize.Get();

  if(_max_apertures) max_apertures=_max_apertures.Get();
  if(_initial_intensity) initial_intensity=_initial_intensity.Get();
  if(_max_intensity) max_intensity=_max_intensity.Get();
  if(_step_intensity) step_intensity=_step_intensity.Get();

  if(_ibo_evals) ibo_evals = _ibo_evals.Get();

  // Initial collimator setup
  if (_setup) {
    string setup = _setup.Get();
    if (setup == "open_max"){
      initial_setup = StationSetup::open_all_max;
      initial_intensity = max_intensity;
    } else if (setup == "open_min"){
      initial_setup = StationSetup::open_all_min;
      initial_intensity = 0;
    } else if (setup == "closed_max") {
      initial_setup = StationSetup::closed_all_max;
      initial_intensity = max_intensity;
    } else if (setup == "closed_min") {
      initial_setup = StationSetup::closed_all_min;
      initial_intensity = 0;
    } else if (setup == "open_min_min") {
      initial_setup = StationSetup::open_min_min;
      initial_intensity = 0;
    } else if (setup == "open_min_k") {
      initial_setup = StationSetup::open_min_k;
    } else if (setup == "random") {
      //TODO: Why these ones are different??
      if(strategy=="dao_ls")
        initial_setup = StationSetup::rand_all_rand;
      else if (strategy=="ibo_ls" || strategy=="ibo+dao")
        initial_setup = StationSetup::rand_int;
    } else if (setup == "manual") {
      initial_setup = StationSetup::manual_all_manual;
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

  if(_tabu_size) tabu_size = _tabu_size.Get();

  // Neighborhood
  if (_nsimple) {
    string nn = _nsimple.Get();
    if (nn == "intensity")
      neighborhood = NeighborhoodType::intensity;
    else if (nn == "aperture")
      neighborhood = NeighborhoodType::aperture;
    else if (nn == "mixed")
      neighborhood = NeighborhoodType::mixed;
    else if (nn == "imixed")
        neighborhood = NeighborhoodType::imixed;
    else if (nn == "smixed")
            neighborhood = NeighborhoodType::smixed_i;
    else
      cout << "Neighborhood operator " <<
          neighborhood << "not recognized"<<endl;
  } else if (_nseq) {
    string nn = _nseq.Get();
    if (nn == "intensity")
      neighborhood = NeighborhoodType::sequential_i;
    else if (nn == "aperture")
      neighborhood = NeighborhoodType::sequential_a;
  } else if (_nprob) {
    prob_intensity = _nprob.Get();
    neighborhood = NeighborhoodType::sequential_p;
  }

  if (_perturbation_type) {
    string nn = _perturbation_type.Get();
    if (nn == "intensity")
      perturbation_type= PerturbationType::p_intensity;
    else if (nn == "aperture")
      perturbation_type = PerturbationType::p_aperture;
    else
      perturbation_type = PerturbationType::p_mixed;
  }
  if (_perturbation_size) perturbation_size = _perturbation_size.Get();

  if (_targeted_search) {
     targeted_search = true;
     target_type = LSTargetType::target_friends;
  }


  // Archivos del problem
  if (_file) file=_file.Get();
  if (_file2) file2=_file2.Get();
  if (_file3) file3=strdup(_file3.Get().c_str()); //intensidades de partida
  if (_path) path=_path.Get();

  int aux=chdir(path.c_str());

  if (_convergence_file) convergence_file=_convergence_file.Get();
  if(convergence_file=="default"){
     string mkdir = "mkdir output";
     int aux=system(mkdir.c_str());
     string base_name= string("output/") + basename(file.c_str()) + "_" + basename(file2.c_str()) + "_"
                        + strategy+"_"+to_string(maxtime)+"_"+to_string(maxeval)+"_"+to_string(neighborhood)
                        + "_"+to_string(initial_setup)+"_"+to_string(perturbation_type)+"_"+to_string(perturbation_size)
                        + "_"+to_string(targeted_search)+"_"+to_string(initial_intensity)+"_"+to_string(max_apertures)
                        + "_"+to_string(step_intensity)+"_"+to_string(max_intensity)+"_"+to_string(ls_type)+"_"+to_string(tabu_size)
						+ "_"+to_string(ibo_evals)+"_"+to_string(seed);
     convergence_file = base_name + ".conv";
     output_file = base_name+".out";
     json_file = base_name+".json";
  }


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
  } else if(strategy=="ibo+dao") {
    cout << "##******** IMRT-Solver (IBO + DAO Local Search) *********"
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
  vector<double> Zmax={65,60,0};

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
  if (initial_setup==StationSetup::open_all_max)
    cout << "open max intensity" << endl;
  else if (initial_setup==StationSetup::open_all_min)
    cout << "open min intensity" << endl;
  else if (initial_setup==StationSetup::closed_all_max)
    cout << "closed max intensity" << endl;
  else if (initial_setup==StationSetup::closed_all_min)
    cout << "closed min intensity" << endl;
  else if (initial_setup==StationSetup::rand_all_rand)
    cout << "random" << endl;
  else if (initial_setup==StationSetup::rand_int)
    cout << "random" << endl;
  else if (initial_setup==StationSetup::open_min_min)
    cout << "open min min intensity" << endl;
  else if (initial_setup==StationSetup::open_min_k)
    cout << "open min k="<< initial_intensity << " intensity" << endl;
  else if (initial_setup==StationSetup::manual_all_manual)
    cout << "manual intensity" << endl;
  cout << "##   Apertures: " << max_apertures << endl;
  cout << "##   Initial intensity: " << initial_intensity << endl;
  cout << "##   Max intensity: " << max_intensity << endl;
  cout << "##   Step intensity: " << step_intensity << endl;

  if (ls_type == LSType::first)
    cout << "##   Local search: first improvement"  << endl;
  if (ls_type==LSType::best)
    cout << "##   Local search: best improvement"  << endl;

  cout << "##   Tabu list size: " << tabu_size << endl;

  if (neighborhood == NeighborhoodType::intensity)
    cout << "##   Neighborhood: intensity" << endl;
  else if (neighborhood == NeighborhoodType::aperture)
    cout << "##   Neighborhood: aperture" << endl;
  else if (neighborhood == NeighborhoodType::mixed)
    cout << "##   Neighborhood: mixed" << endl;
  else if (neighborhood == NeighborhoodType::sequential_i)
    cout << "##   Neighborhood: sequential intensity first" << endl;
  else if (neighborhood == NeighborhoodType::sequential_a)
    cout << "##   Neighborhood: sequential aperture first" << endl;
  else if (neighborhood == NeighborhoodType::sequential_p) {
    cout << "##   Neighborhood: sequential probabilistic: " <<
            prob_intensity << endl;
  }

  if (targeted_search)
    cout << "##   Targeted search: yes (friends)"  << endl;
  else
    cout << "##   Targeted search: no"  << endl;

  if (perturbation_type == PerturbationType::p_intensity)
    cout << "##   Perturbation: intensity " << endl;
  else if (perturbation_type == PerturbationType::p_aperture)
    cout << "##   Perturbation: aperture " << endl;
  else if (perturbation_type == PerturbationType::p_mixed)
    cout << "##   Perturbation: mixed " << endl;

  cout << "##   Perturbation size: " << perturbation_size << endl;

  cout << "##" << endl << "## Colimator configuration: "<< endl;
  cout << "##   Stations: " << collimator.getNbAngles() << endl;
  cout << "##   Angles: ";
  for (int i=0; i<collimator.getNbAngles();i++)
    cout << collimator.getAngle(i) << " ";
  cout << endl;


  cout << "##" << endl << "## Instance information: "<< endl;
  cout << "##   Volumes: " << volumes.size() << endl;

  cout << "##" << endl << "## Output files: " << endl;
  if (convergence_file!="")
    cout << "##   Convergence file: " << convergence_file << endl;

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




  ILS* ils;
  double cost;
  if (strategy=="dao_ls") {
    ils = new ApertureILS(bsize, vsize, prob_intensity, step_intensity);
    cost = ils->iteratedLocalSearch(P, maxtime, maxeval, ls_type, neighborhood, target_type, perturbation_type, perturbation_size, tabu_size, convergence_file);
  }else if(strategy=="ibo_ls"){

	  ils = new IntensityILS2();
	  cost = ils->iteratedLocalSearch(P, maxtime, maxeval, ls_type, neighborhood, target_type, perturbation_type, perturbation_size, tabu_size, convergence_file);

  }else if(strategy=="ibo+dao"){
    ils = new IntensityILS2();
	  cost = ils->iteratedLocalSearch(P, maxtime, ibo_evals, ls_type, neighborhood, target_type, perturbation_type, perturbation_size, tabu_size, convergence_file);
	  P.generateApertures();
	  int evals=ils->used_evaluations;
	  std::clock_t begin=ils->time_begin;

	  ils = new ApertureILS(bsize, vsize, prob_intensity, step_intensity);
    if(neighborhood == NeighborhoodType::imixed) neighborhood=NeighborhoodType::mixed;
	  cost = ils->iteratedLocalSearch(P, maxtime, maxeval, ls_type, neighborhood, target_type, perturbation_type,
			  perturbation_size, tabu_size, convergence_file, evals, begin);

  }




  cout << "##**************************************************************************"
       << endl;
  cout << "##******************************* RESULTS **********************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;
  cout << "##"<<endl;
  cout << "## Best solution found: " <<  cost << endl;


	cout << endl;
	for(int i=0;i<5;i++)
		P.printIntensity(i);

	cout << endl;

  ofstream o_file,j_file;
  if(convergence_file!=""){
    o_file.open (output_file.c_str(), ios::out);
    j_file.open (json_file.c_str(), ios::out);
  }


  cout <<  cost << " ";
  if(convergence_file!="") o_file  <<  cost << " ";

  std::clock_t time_end = clock();
  double used_time = double(time_end - ils->time_begin) / CLOCKS_PER_SEC;

  cout <<  used_time << " ";
  if(convergence_file!="") o_file  <<  used_time  << " ";

  if(convergence_file!="") o_file << ils->used_evaluations << " ";

  const list<Station*> stations=P.get_stations();

  int tot_alpha=0;
  for(auto s:stations){
    string str=(strategy=="ibo+dao")? "dao_ls":strategy;
    int alpha=s->get_sum_alpha(str);
    cout << alpha << " " ;
    //if(convergence_file!="") o_file  << alpha << " " ;
    tot_alpha+=alpha;
  }
  cout << tot_alpha << " ";
  if(convergence_file!="") o_file << tot_alpha << " ";

  int nb_apertures=0;
  for(auto s:stations){
    string str=(strategy=="ibo+dao")? "dao_ls":strategy;
    int ap=s->get_nb_apertures(str);
    cout << ap << " " ;
    //if(convergence_file!="") o_file << ap << " " ;
    nb_apertures+=ap;
  }
  cout << nb_apertures << endl;
  if(convergence_file!="") o_file << nb_apertures << endl;
  if(convergence_file!="") o_file.close();

  set<int> l = get_angles(file, 5);
  if(_plot){
	  std::stringstream ss;
	  ss << "python plotter/plot.py " << *l.begin()/5 << "_" << strategy << "_" << initial_intensity << "_" << initial_setup << "_" << maxtime << "_" << seed <<".pdf";
	  std::string s = ss.str();
	  system(s.c_str());
  }

  if(convergence_file!=""){
    bool flag=false;
    j_file << "{" << endl;
    for(auto s:stations){
      if (flag) j_file << "," << endl;
      flag=true;
      j_file <<  s->toStringIntensities() ;
    }
    j_file << endl << "}" << endl;
    j_file.close();
  }
  return 0;


	/*cout << "********   Summary of the results    *********"<< endl;
	best_eval = F.eval(P,w,Zmin,Zmax);
	cout << "Final solution: " << best_eval << endl << endl;

  F.generate_voxel_dose_functions ();*/





	return 0;

}
