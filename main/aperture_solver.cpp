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

/*
double doClose(int beam, int a, Station& station, double c_eval, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax, EvaluationFunction& F, Plan& P) {
  double aux_eval, t_eval, f_eval;
  list<pair<int, double> > aux_diff;

  aux_diff = station.closeBeamlet(beam, a, true);
  if (aux_diff.size()>=1) {
    c_eval = t_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    cout << "  Closing left eval: " << c_eval << " list: ";
    for (list<pair<int,double>>::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
    cout << endl;
    aux_diff = station.undoLast();
    aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
  }

  aux_diff = station.closeBeamlet(beam, a, false);
  if (aux_diff.size()>=1) {
    f_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    cout << "  Closing right eval: " << f_eval << " list: ";
    for (list<pair<int,double>>::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
    //cout << endl;

    if ( f_eval > t_eval ) {
      aux_diff = station.undoLast();
      aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
      aux_diff = station.closeBeamlet(beam, a, true);
      c_eval   = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    }else{
      c_eval=f_eval;
    }
  }
  return(c_eval);
}

double doOpen(int beam, int a, Station& station, double c_eval, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax, EvaluationFunction& F, Plan& P) {
  double aux_eval=0;
  list<pair<int, double> > aux_diff;

  aux_diff = station.openBeamlet(beam, a);
  if(aux_diff.size() <1) return(c_eval);

  aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);

  cout << "  Opening eval: " << aux_eval << " size: " << aux_diff.size() << " list: ";
  for (list<pair<int,double>>::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
  //cout << endl;

  return(aux_eval);
}

double searchFirstAperture(int beam, Station& station, double c_eval, vector<double>& w,
                           vector<double>& Zmin, vector<double>& Zmax, EvaluationFunction& F,
                           bool open_beamlet, Plan & P, double temperature) {
  double current=c_eval, aux_eval=0, local_best=-1;
  list<pair<int, double> > aux_diff;
  int a, last_a, local_a;
  bool flag=true;
  vector<int> a_list;

  if (!open_beamlet) a_list = station.getOpen(beam);
  else a_list = station.getClosed(beam);

  if (a_list.size()<1) {
    //tabu_list.push_back(make_pair(beam, make_pair(&station, open_beamlet)));
    return(c_eval);
  }

  a = last_a = (int)rand() % a_list.size();

  // Check every aperture
  while(flag) {
    if (!open_beamlet) {
      aux_eval = doClose(beam, a_list[a], station, c_eval, w, Zmin, Zmax, F, P);
    } else {
      aux_eval = doOpen(beam, a_list[a], station, c_eval, w, Zmin, Zmax, F, P);
    }

    if (local_best<0 || aux_eval < local_best){
      local_best = aux_eval;
      local_a = a;
    }

    if (current > aux_eval) {
      current = aux_eval;
      last_a=a;
      return(current);
    } else {
      aux_diff = station.undoLast();
      if (aux_diff.size()>0)
        aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    }

    if (a==(a_list.size()-1))  a=0;
    else  a+=1;
    if (last_a==a) flag=false;
  }

  // SA acceptance criterion
  double p = exp((double)-(local_best-c_eval)/temperature);
  double r = ((double)rand() / (RAND_MAX));
  cout << ", eval: " << p << " " << r ;
  if (abs(current-c_eval)<0.001 & abs(local_best-c_eval)>0.001 &  p > r) {
    cout << ", accept worst ";
    if (!open_beamlet) {
      current = doClose(beam, a_list[local_a], station, c_eval, w, Zmin, Zmax, F, P);
    } else {
      current = doOpen(beam, a_list[local_a], station, c_eval, w, Zmin, Zmax, F, P);
    }
  }

  //if (abs(current-c_eval)<0.001)
    //tabu_list.push_back(make_pair(make_pair(beam,&station), make_pair(a_list[local_a], open_beamlet)));

  return(current);
}

double searchFirstIntensity(int beam, Station& station, double c_eval, vector<double>& w,
                            vector<double>& Zmin, vector<double>& Zmax, EvaluationFunction& F,
                            bool open_beamlet, Plan & P, double temperature) {
  double current = c_eval, aux_eval=0, local_best=-1;
  list<pair<int, double> > aux_diff;
  int a= (int) rand()% station.getNbApertures();
  int last_a=a, local_a;
  bool flag=true;

  // Check every aperture
  while(flag) {
    // Left op
    if (open_beamlet) {
      aux_diff = station.modifyIntensityAperture(a, 1);
      aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);

      //cout << endl<<"  Increasing intensity aperture: " << a <<" eval: " << aux_eval << " list: ";
      //for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      //cout << endl;
    } else {
      aux_diff = station.modifyIntensityAperture(a, -1);
      aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);

      //cout << endl<< "  Reducing intensity aperture: " << a <<" eval: " << aux_eval <<  " list: ";
      //for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      //cout << endl;
    }

    if (local_best<0 || local_best> aux_eval){
      local_best=aux_eval;
      local_a=a;
    }

    if (current > aux_eval) {
      current =aux_eval;
      last_a = a;
      return(current);
    } else {
      aux_diff = station.undoLast();
      aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
      if (last_a==a)  flag=false;
    }

    if (a==(station.getNbApertures()-1)) a=0;
    else a+=1;
  }

  double p = exp((double)-(local_best-c_eval)/temperature);
  double r = ((double)rand() / (RAND_MAX));
  cout << ", p " << p << ", r: " << r;
  if (abs(current-c_eval) < 0.001 && abs(local_best-c_eval)>0.001 && p > r)  {
    cout << ", accept worst ";
    if (!open_beamlet) {
      aux_diff = station.modifyIntensityAperture(local_a, 1);
      current = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    } else{
      aux_diff = station.modifyIntensityAperture(local_a, -1);
      current = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    }
  }

  //if (abs(current-c_eval)<0.001)
  //  tabu_list.push_back(make_pair(make_pair(beam,&station),open_beamlet));

  return(current);
}

double coolDownTemperature(double t, double alpha, int iteration) {
  double new_t;
  //exponential
  //new_t=t*pow(alpha,iteration);
  //logarithmic
  new_t=t/(1+log(iteration+1));
  return(new_t);

}

bool isModifiable(int beamlet, Station* station, bool open_flag, bool intensity_flag, bool aperture_flag){
  if (!open_flag) {
    if (station->anyOpen(beamlet) && aperture_flag)
      return(true);
    if (station->canReduceIntensity(beamlet) && intensity_flag)
      return(true);
  } else {
    if (station->anyClosed(beamlet) && aperture_flag)
      return(true);
    if (station->canIncreaseIntensity(beamlet) && intensity_flag)
      return(true);
  }
  return(false);
}

pair <int, pair<Station*, bool> > getLSBeamlet(Plan& P, EvaluationFunction& F, int bsize, int vsize, bool check_intensity, bool check_aperture) {
  auto sb = F.best_beamlets(P, bsize, vsize);
  auto it = sb.begin();
  std::advance(it,rand()%sb.size());
  Station * s = it->second.first;
  int beamlet = it->second.second;
  bool sign = it->first.second;
  int i=0;
  while (!isModifiable(beamlet, s, !sign, check_intensity, check_aperture)) {
    //cout << "Rejecting " << beamlet << endl;
    if (i==sb.size()) return(make_pair(-1, make_pair(s,false)));
    if (it==sb.end()) it = sb.begin();
    else it++;
    i++;
  }
  s = it->second.first;
  beamlet = it->second.second;
  sign = it->first.second;
  return(make_pair(beamlet, make_pair(s,sign)));
}*/


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
	     ils = new ApertureILS(bsize, vsize, search_intensity, search_aperture, prob_intensity, step_intensity, initial_temperature, alphaT, acceptance);
  else if(strategy=="ibo_ls")
      ils = new IntensityILS(bsize, vsize, maxdelta, maxratio, alpha, beta);

  ils->search(P, maxtime, maxiter);
  
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
