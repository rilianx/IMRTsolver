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
#include "args.hxx"


using namespace std;
using namespace imrt;

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

double searchFirstAperture(int beam, Station& station, double c_eval, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax, EvaluationFunction& F, bool open_beamlet, Plan & P, double temperature) {
  double current=c_eval, aux_eval=0, local_best=-1;
  list<pair<int, double> > aux_diff;
  int a= (int)rand() % station.getNbApertures();
  int last_a = a, local_a;
  bool flag=true; 
  
  // Check every aperture 
  while(flag) {
    if (!open_beamlet) {
      // Left op
      aux_eval = doClose(beam, a, station, c_eval, w, Zmin, Zmax, F, P); 
    } else {
      aux_eval = doOpen(beam, a, station, c_eval, w, Zmin, Zmax, F, P);
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
        
    if (a==(station.getNbApertures()-1)) 
      a=0;
    else 
      a+=1;
    
    if (last_a==a) flag=false;
  }
  
  // SA acceptance criterion
  if (current==c_eval && abs(local_best-c_eval)>0.001 &&  exp((double)-(local_best-c_eval)/temperature) > ((double)rand() / (RAND_MAX))) {
    cout << ", accept worst ";
    if (!open_beamlet) {
      current = doClose(beam, local_a, station, c_eval, w, Zmin, Zmax, F, P); 
    } else {
      current = doOpen(beam, local_a, station, c_eval, w, Zmin, Zmax, F, P);
    }
  }
  
  return(current);
}

double searchFirstIntensity(int beam, Station& station, double c_eval, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax, EvaluationFunction& F, bool open_beamlet, Plan & P, double temperature) {
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
      /*dev = F.eval(P, w, Zmin, Zmax);
      cout << endl<<"  Increasing intensity aperture: " << a <<" eval: " << aux_eval << " real: "<< dev << " list: ";
      for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      cout << endl;*/
    } else {
      aux_diff = station.modifyIntensityAperture(a, -1);
      aux_eval = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
      /*dev = F.eval(P, w, Zmin, Zmax);
      cout << endl<< "  Reducing intensity aperture: " << a <<" eval: " << aux_eval << " real: "<< dev << " list: ";
      for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      cout << endl;*/
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
  
  if (current==c_eval && abs(local_best-c_eval)>0.001 && exp((double)-(local_best-c_eval)/temperature) > ((double)rand() / (RAND_MAX)))  {
    cout << ", accept worst ";
    if (!open_beamlet) {
      aux_diff = station.modifyIntensityAperture(local_a, 1);
      current = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    } else{
      aux_diff = station.modifyIntensityAperture(local_a, -1);
      current = F.incremental_eval(station, w, Zmin, Zmax, aux_diff);
    }
  }
  
  return(current);
}

double coolDownTemperature(double t, double alpha, int iteration) {
  double new_t;
  //exponential
  //new_t=t*pow(alpha,iteration);
  //logarithmic
  new_t=t/log(iteration+1);
  return(new_t);
  
}

int main(int argc, char** argv){

  int vsize=20;
  int bsize=5;
  double int0=4.0;
  int maxiter=10000;
  int max_apertures=5;
  double alpha=1.0;
  double beta=1.0;
  double maxdelta=5.0;
  double maxratio=3.0;
  
 
  int initial_intensity=1;
  
  // ls params
  double prob_intensity=0.05;
  double temperature=10;
  double min_temperature=0;
  
  int seed=time(NULL);


	args::ArgumentParser parser("********* IMRT-Solver (Aperture solver) *********", "An IMRT Solver.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	//args::ValueFlag<string> _format(parser, "string", "Format: (BR, BRw, 1C)", {'f'});
	args::ValueFlag<int> _bsize(parser, "int", "Number of considered beamlets for selection ("+to_string(bsize)+")", {"bsize"});
	args::ValueFlag<int> _vsize(parser, "int", "Number of considered worst voxels ("+to_string(vsize)+")", {"vsize"});
  args::ValueFlag<int> _int0(parser, "int", "Initial intensity for beams  ("+to_string(int0)+")", {"int0"});
  args::ValueFlag<int> _max_apertures(parser, "int", "Initial intensity for the station  ("+to_string(max_apertures)+")", {"max_ap"});
  args::ValueFlag<int> _maxdelta(parser, "int", "Max delta  ("+to_string(maxdelta)+")", {"maxdelta"});
  args::ValueFlag<int> _maxratio(parser, "int", "Max ratio  ("+to_string(maxratio)+")", {"maxratio"});
  args::ValueFlag<double> _alpha(parser, "double", "Initial temperature for intensities  ("+to_string(alpha)+")", {"alpha"});
  args::ValueFlag<double> _beta(parser, "double", "Initial temperature for ratio  ("+to_string(beta)+")", {"beta"});
  args::ValueFlag<int> _maxiter(parser, "int", "Number of iterations ("+to_string(maxiter)+")", {"max_iter"});
  args::ValueFlag<int> _seed(parser, "int", "Seed  ("+to_string(seed)+")", {"seed"});
  args::ValueFlag<int> _initial_intensity(parser, "int", "Initial value aperture intensity  ("+to_string(initial_intensity)+")", {"initial_intensity"});
  args::Flag open_setup(parser, "open_setup", "Initialize apertures as open", {"open_setup"});
  args::Group ls_type (parser, "Local search type:", args::Group::Validators::Xor);
  args::Flag ls_apertures(ls_type, "ls_apertures", "Apply aperture local search", {"ls_apertures"});
  args::Flag ls_intensity(ls_type, "ls_intensity", "Apply intensity local search", {"ls_intensity"});
  args::Flag ls_both (ls_type, "ls_both", "Apply intensity and aperture local search based in probabilities", {"ls_both"});
  args::ValueFlag<double> _prob_intensity(parser, "double", "Probability to search over intensity  ("+to_string(prob_intensity)+")", {"prob_intensity"});
  args::ValueFlag<double> _temperature(parser, "double", "Temperature for acceptance criterion  ("+to_string(temperature)+")", {"temperature"});
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
  if(_maxdelta) maxdelta=_maxdelta.Get();
  if(_maxratio) maxratio=_maxratio.Get();
  if(_alpha) alpha=_alpha.Get();
  if(_beta) beta=_beta.Get();
  if(_int0) int0=_int0.Get();
  if(_maxiter) maxiter=_maxiter.Get();
  if(_max_apertures) max_apertures=_max_apertures.Get();
  if(_seed) seed=_seed.Get();
  if(_initial_intensity) initial_intensity=_initial_intensity.Get();
  if(_prob_intensity) prob_intensity=_prob_intensity.Get();
  if(_temperature) temperature=_temperature.Get();

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
	  station = new Station(collimator,volumes, i*70, max_apertures, initial_intensity, open_setup);
	  station->generateIntensity();
	  stations[i]=station;
  }
  
  vector<double> w={1,1,1};
  vector<double> Zmin={0,0,76};
  vector<double> Zmax={65,60,1000};
  
  cout << "************************************************"<< endl;
  cout << "************************************************"<< endl;
  cout << "******** IMRT-Solver (Aperture solver) *********"<< endl;
  cout << "Iterations: " << maxiter << endl;
  cout << "Seed: " << seed << endl;
  cout << "Temperature: " << temperature << endl;
  if (ls_apertures)
    cout << "Searching: aperture pattern" << endl;
  if (ls_intensity)
    cout << "Searching: intensity" << endl;
  if (ls_both){
    cout << "Searching: intensity and aperture pattern" << endl;
    cout << "Probability intensity ls: " << prob_intensity << endl;
  }
  
  cout << endl << "Colimator configuration: "<< endl;
  cout << "  Stations: " << stations.size() << endl;
  cout << "  Angles: ";
  for (int i=0; i<stations.size();i++) cout << stations[i]->getAngle(); " ";
  cout << endl;
  cout << "  Max apertures: " << max_apertures << endl;
  cout << "  Initial intensity: " << initial_intensity << endl;
  cout << "  Open initial setup: " << open_setup << endl;

  
  cout << endl << "Instance information: "<< endl;
  cout << "  Volumes: " << volumes.size() << endl;

  cout << "************************************************"<< endl<< endl;
  cout << "************************************************"<< endl;
  cout << "************************************************"<< endl;
  
  EvaluationFunction F(volumes);
  Plan P(F, w, Zmin, Zmax);
  for(int i=0;i<5;i++)
	  P.add_station(*stations[i]);
	double best_eval=P.doEval();
	
	cout << "Initial solution: " << best_eval << endl;
	
	//From here 
	/*auto sb=F.best_beamlets(P, bsize, vsize);
	auto it=sb.begin();
	Station*s = it->second.first; 
	int beamlet=it->second.second;
	bool sign=it->first.second;*/


	for (int i=0;i<maxiter;i++) {
		auto sb=F.best_beamlets(P, bsize, vsize);
		auto it=sb.begin();
		std::advance(it,rand()%sb.size());

		Station*s = it->second.first; 
		int beamlet=it->second.second;
		bool sign=it->first.second; //impact in F (+ or -)

    cout << "Iteration " << (i+1) << ", best: " << best_eval << ", beamlet: " << beamlet  << ", station: " << s->getAngle() << ", d: " << sign << " t: " << temperature;
    if (ls_apertures) {
      best_eval = searchFirstAperture(beamlet, *s, best_eval, w, Zmin, Zmax, F, !sign, P, temperature);
      cout << ", ls aperture, found: " << best_eval;
    } else if (ls_intensity) {
      best_eval = searchFirstIntensity(beamlet, *s, best_eval, w, Zmin, Zmax, F, !sign, P, temperature);
      cout << ", ls intensity, found: " << best_eval;
    } else if (ls_both) {
      if (((double) rand() / (RAND_MAX)) <= prob_intensity){
        best_eval = searchFirstIntensity(beamlet, *s, best_eval, w, Zmin, Zmax, F, !sign, P, temperature);
        cout << ", ls intensity, found: " << best_eval;
      } else {
	      best_eval = searchFirstAperture(beamlet, *s, best_eval, w, Zmin, Zmax, F, !sign, P, temperature);
        cout << ", ls aperture, found: " << best_eval;
      }
    }
    cout  <<endl;
    //s->printIntensity(false);
    temperature=coolDownTemperature(temperature, 1, i+1); 
	}

	/*cout << endl;
	for(int i=0;i<5;i++){
		//stations[i]->printIntensity();
		stations[i]->printIntensity(false);
        //cout << "nb_apertures:" << stations[i]->int2nb.size() << endl;
    }
	cout << endl;*/

	
	cout << "********   Summary of the results    *********"<< endl;
	best_eval = F.eval(P,w,Zmin,Zmax);
	cout << "Final solution: " << best_eval << endl << endl;

  F.generate_voxel_dose_functions ();
  system("python plotter/plot.py");

	return 0;

}
