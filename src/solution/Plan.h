/*
 * Schedule.h
 *
 *  Created on: 7 may. 2018
 *      Author: iaraya
 */

#ifndef PLAN_H_
#define PLAN_H_

#include "EvaluationFunction.h"
#include "Station.h"

namespace imrt {


/* An IMRT plan
 *
 * It consists in a list of stations, each station corresponds to an
 * angle, and a set of apertures (or a matrix of intensities)
 */

class Plan {
public:
	Plan(EvaluationFunction &ev);

  Plan(EvaluationFunction &ev, vector<double> w, vector<double> Zmin, vector<double> Zmax);

  Plan(vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& collimator, 
       vector<Volume>& volumes, int max_apertures, int max_intensity, int initial_intensity, 
       int step_intensity=2, int setup=6, char* file=NULL);

  Plan(const Plan &p);

	virtual ~Plan() {};

	void newCopy(Plan& p);

	// Adds a new station to the plan
	void add_station(Station& s);

	double eval(vector<double>& w, vector<double>& Zmin, vector<double>& Zmax);

	double eval();

	double incremental_eval (Station& station, int i, int j, double intensity){
		list< pair< int, double > > diff;
		diff.push_back(make_pair(station.pos2beam[make_pair(i,j)],intensity));
		return incremental_eval (station, diff);
	}

	double incremental_eval (Station& station, list< pair< int, double > >& diff);

	double get_delta_eval (Station& s, int i, int j, double intensity, int n_voxels=999999){
    return ev.get_delta_eval(s.getAngle(),
    s.pos2beam.at(make_pair(i,j)), intensity, w, Zmin, Zmax, n_voxels);
  }
	
	double get_delta_eval (Station& s, int b, double intensity, int n_voxels=999999){
	  return ev.get_delta_eval(s.getAngle(), b, intensity, w, Zmin, Zmax, n_voxels);
	}
	
	double get_delta_eval (Station& s, list< pair< int, double > >& diff, int n_voxels=999999){
	  return ev.get_delta_eval(diff, s.getAngle(), w, Zmin, Zmax, n_voxels);
	};

	// This function assumes that there are no changes made without evaluation
	// performed with eval or incrementalEval
	double getEvaluation();

	const list<Station*>& get_stations() const;

	void write_open_beamlets();

	set < pair< pair<double,bool>, pair<Station*, int> >,
       std::greater < pair< pair<double,bool>, pair<Station*, int> > > >
	  best_beamlets(int n, int nv, int mode=0);

	virtual pair<bool, pair<Station*, int>> getLSBeamlet(int bsize, int vsize){
		  auto sb=ev.best_beamlets(*this, bsize, vsize);
		  auto it=sb.begin();
		  std::advance(it,rand()%sb.size());
		  return make_pair(it->first.second, it->second);
	}

	virtual pair<bool, pair<Station*, int>> getBestLSBeamlet(int bsize, int vsize){
	  auto sb=ev.best_beamlets(*this, bsize, vsize);
	  auto it=sb.begin();
	  //std::advance(it,rand()%sb.size());
	  return make_pair(it->first.second, it->second);
	}

  double openBeamlet(int station, int aperture, int beamlet);

  double closeBeamlet(int station, int aperture, int beamlet, int side);

  double modifyIntensityAperture (int station, int aperture, int delta);

	void undoLast();

  void clearLast();

	//Lepi's version
	void undoLast2();

  void printSolution ();

  void printIntensity(int station);

	void writeIntensities(string file, int n);

	EvaluationFunction* getEvaluationFunction(){
		return &ev;
	}

	Station* get_station(int i){
	    list<Station*>::iterator s= stations.begin();
	    advance(s,i);
	    return *s;
	}

  int getNStations () {
     return(stations.size()) ;
  }

  //update Z by increasing the intensity of beamlet (angle,b) in delta_intensity
  //if return_if_unfeasible=true, then it returns when some organ voxel surpasses Zmax
  //return false if some voxel surpasses Zmax
  bool Zupdate(Station* s, int b, double delta_intensity, bool return_if_unfeasible){
	  return ev.Zupdate(s->getAngle(), b, delta_intensity, return_if_unfeasible, Zmax);
  }

  //regresa al savepoint para Z
  void Zrollback(){
  	ev.Zrollback();
  }

  void Zsavepoint(){
  	ev.Zsavepoint();
  }

  //almacena en sorted_beamlets los beamlet (station, b) ordenados segun v/c
  void get_vc_sorted_beamlets(Plan& p, multimap < double, pair<Station*, int>, MagnitudeCompare >& sorted_beamlets){
  	ev.get_vc_sorted_beamlets(p, Zmin, Zmax, sorted_beamlets);
  }


  map<int, Station*> angle2station;
private:
	//The list of stations
	//list<Station> real_stations;
	list<Station*> stations;

  int n_stations;

  Station* last_changed;
	list< pair< int, double > > last_diff;

	EvaluationFunction& ev;

	vector<double> w;
	vector<double> Zmin;
	vector<double> Zmax;

	double evaluation_fx;
};

} /* namespace imrt */

#endif /* PLAN_H_ */
