/*
 * main.cpp
 *
 *  Created on:
 *      Author: iaraya, leslie
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
#include "IntensityILS2.h"
#include "MixedILS.h"
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


vector<Volume> createVolumes (string organ_filename, Collimator& collimator, int max_voxels){
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
    volumes.push_back(Volume(collimator, organ_files[i], max_voxels));

  return(volumes);
}


int main(int argc, char** argv){

  int seed = time(NULL);

  // Budget and execution variables
  int maxiter = 5000;
  int ibo_evals = 2500;
  int maxtime = 0;
  int maxeval = 0;

  // Aperture and intensity configuration
  StationSetup initial_setup = StationSetup::open_min_k;;
  int max_apertures = 5;
  int I_0 = 5;
  int max_intensity = 28;
  int step_intensity = 1;

  // Type of local search
  string strategy = "dao_ls";
  bool continuous = true;
  double prob_intensity = 0.2;
  int tabu_size = 0;

  // Target beamlet heuristic
  bool targeted_search = false;
  double vsize = 0.002;
  double min_improvement=0.05;


  // Files
  string path = ".";
  string file = "data/testinstance_0_70_140_210_280.txt";
  string file2 = "data/test_instance_coordinates.txt";
  char* file3 = NULL;
  int max_voxels=100000;


  args::ArgumentParser parser("********* IMRT-Solver (Aperture solver) *********",
                             "Example.\n../AS -s ibo_ls --maxeval=4000 --ls_sequential=intensity --setup=open_min --seed=2 --ls=first  --max-intensity=20");

  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<int>    _seed     (parser, "int", "Seed  (" +
                                     to_string(seed)+")", {"seed"});

  // Search strategies
  args::Group strat (parser, "Strategy options:");
  args::ValueFlag<int>    _tabu_size (strat, "int",
                                    "Tabu list size(" +
                                     to_string(tabu_size)+")", {"tabu-size"});


  // Execution parameters
  args::Group budget (parser, "Budget options:");
  args::ValueFlag<int>    _maxtime  (budget, "int",
                                    "Maximum time in seconds (" +
                                     to_string(maxtime)+")", {"maxtime"});
  args::ValueFlag<int>    _maxeval  (budget, "int",
                                    "Number of evaluations (" +
                                     to_string(maxiter)+ ")", {"maxeval"});


  // Initial collimator setup (initial solution: aperture, intensity)
  args::Group isetup (parser, "Initial collimator setup:");
  args::ValueFlag<int> _I_0 (isetup, "int",
                                  "Initial intensities in the fluence map matrices " + to_string(I_0), {"setup"});
  args::ValueFlag<int> _max_apertures     (isetup, "int",
                                          "Number of apertures per angle (station) (" +
                                           to_string(max_apertures)+")", {"max-apertures"});

  // Target beamlet heuristic parameters
  args::Group heur (parser, "Parameters for the ordered neighbourghood (targeted):");
  args::ValueFlag<double>    _vsize    (heur, "double",
                                    "Percentage of considered worst voxels (" +
                                     to_string(vsize)+")", {"vsize"});
 args::ValueFlag<double>    _min_improvement    (heur, "double",
                                   "Minimum beamlet improvement estimation(" +
                                    to_string(min_improvement)+")", {"min_impr"});

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
  args::ValueFlag<string> _bac  (io_opt, "string",
                                 "Beam angle configuration ", {"bac"});     
  args::ValueFlag<int> _max_voxels  (io_opt, "int",
                                 "Maximum number of voxels per organ ", {"max_voxels"});                        

  // Output file parameters
  args::ValueFlag<string> _convergence_file (io_opt, "string",
                                 "File to output convergence", {"convergence"});

	try { parser.ParseCLI(argc, argv); }
	catch (args::Help&) {  std::cout << parser; return 0; }
	catch (args::ParseError& e)  { std::cerr << e.what() << std::endl;  std::cerr << parser; return 1; }
	catch (args::ValidationError& e) { std::cerr << e.what() << std::endl; std::cerr << parser; return 1; }

  //setting parameters

	if(_maxtime) maxtime = _maxtime.Get();
  if(_maxeval) maxeval = _maxeval.Get();
	if(_seed) seed = _seed.Get();

  // parameters for the ordered neighbourhood (targeted)
	if(_vsize) vsize = _vsize.Get();
  if(_min_improvement) min_improvement=_min_improvement.Get();

  if(_max_apertures) max_apertures=_max_apertures.Get();
  if(_I_0) I_0=_I_0.Get();

  // Local search strategy
  if(_tabu_size) tabu_size = _tabu_size.Get();
 
  // Input files
  if (_file) file=_file.Get();
  if (_file2) file2=_file2.Get();
  if (_file3) file3=strdup(_file3.Get().c_str()); //intensidades de partida
  if (_path) path=_path.Get();

  int aux=chdir(path.c_str());

  string base_name="output";
  if (_convergence_file) base_name=string("output/") + _convergence_file.Get();
  string convergence_file = base_name + "/" + basename(file.c_str()) + "_" + to_string(seed) + ".conv";

  cout << "##**************************************************************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;

  cout << "##******** IMRT-Solver (DAO-ILS) *********" << endl;

  cout << "##**************************************************************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;

  // Iniciar generador de numeros aleatorios
  srand(seed);

  vector<double> w={1,1,5};
  vector<double> Zmin={0,0,76};
  vector<double> Zmax={65,65,76};

  IntensityILS2::vsize=vsize;
  IntensityILS2::min_improvement=min_improvement;

  if(_max_voxels) max_voxels=_max_voxels.Get();

  // Create colimator object and volumes
  Collimator collimator (file2, get_angles(file, 5));
  vector<Volume> volumes = createVolumes (file, collimator, max_voxels);

  // Create an initial plan
  cout << "creating the initial treatment plan" << endl;

  vector<int> bac;

	if(_bac){
		std::istringstream bac_stream(_bac.Get());
    int a;
    while (bac_stream >> a ) bac.push_back(a);
	}

  Plan P (w, Zmin, Zmax, collimator, volumes, max_apertures,
          max_intensity, I_0, step_intensity,
          initial_setup, file3, bac);

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
  cout << "##   Apertures: " << max_apertures << endl;
  cout << "##   Initial intensity: " << I_0 << endl;
  cout << "##   Max intensity: " << max_intensity << endl;
  cout << "##   Step intensity: " << step_intensity << endl;
  cout << "##   Local search: first improvement"  << endl;
  cout << "##   Tabu list size: " << tabu_size << endl;
  cout << "##   Targeted search: yes (friends)"  << endl;
  cout << "##   Perturbation: mixed " << endl;
  cout << "##   Perturbation size: 3" << endl;

  cout << "##" << endl << "## Colimator configuration: "<< endl;
  cout << "##   Stations: " << collimator.getNbAngles() << endl;
  cout << "##   Angles: ";
  for (int i=0; i<collimator.getNbAngles();i++)
    cout << collimator.getAngle(i) << " ";
  cout << endl;


  cout << "##" << endl << "## Instance information: "<< endl;
  cout << "##   Volumes: " << volumes.size() << endl;

  cout << "##" << endl << "## Output files: " << endl;

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
  for(int i=0;i<P.getNStations();i++)
    P.printIntensity(i);
  cout << endl;




  double cost;
  double used_evaluations = 0;
  std::clock_t begin_time = clock();


  MixedILS mixed_ils(0, 0, 0, 1 /*step_intensity*/);
  cost = mixed_ils.iteratedLocalSearch(P, maxtime, maxeval, LSType::first, LSType::first, 
          false /*continuous*/, NeighborhoodType::sequential_a, NeighborhoodType::sequential_i, LSTargetType::target_friends,
          LSTargetType::target_none, PerturbationType::p_mixed, 3 /*perturbation_size*/,
          tabu_size, convergence_file);
  used_evaluations = mixed_ils.total_evals;

  cout << "##**************************************************************************"
       << endl;
  cout << "##******************************* RESULTS **********************************"
       << endl;
  cout << "##**************************************************************************"
       << endl;
  cout << "##"<<endl;
  cout << "## Best solution found: " <<  cost << endl; //<< " "<< P.eval() << endl;

	cout << endl;
	for(int i=0;i<P.getNStations();i++)
		P.printIntensity(i, true);

	cout << endl;

  cout <<  cost << " ";

  std::clock_t time_end = clock();
  double used_time = double(time_end - begin_time) / CLOCKS_PER_SEC;

  cout <<  used_time << " " << used_evaluations << " ";

  const list<Station*> stations=P.get_stations();

  int tot_alpha=0;
  for(auto s:stations){
    string str="dao_ls";
    int alpha=s->get_sum_alpha(str);
    cout << alpha << " " ;
    tot_alpha+=alpha;
  }
  cout << tot_alpha << " ";


  int nb_apertures=0;
  for(auto s:stations){
    string str="dao_ls";
    int ap=s->get_nb_apertures(str);
    cout << ap << " " ;
    nb_apertures+=ap;
  }
  cout << nb_apertures << endl;

  return 0;


	/*cout << "********   Summary of the results    *********"<< endl;
	best_eval = F.eval(P,w,Zmin,Zmax);
	cout << "Final solution: " << best_eval << endl << endl;

  F.generate_voxel_dose_functions ();*/


	return 0;

}
