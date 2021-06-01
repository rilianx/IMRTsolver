/*
 * Schedule.h
 *
 *  Created on: 7 may. 2018
 *      Author: iaraya
 */

#ifndef PLAN_H_
#define PLAN_H_

#include "Station.h"
#include "Evaluator.h"
#include <float.h>

namespace imrt {


/* An IMRT plan
 *
 * It consists in a list of stations, each station corresponds to an
 * angle, and a set of apertures (or a matrix of intensities)
 */

class Plan {
public:

  Plan(Evaluator& ev, Collimator& collimator, vector<Volume>& volumes, int max_apertures,
       int max_intensity, int initial_intensity, int step_intensity=2,
       StationSetup setup = open_all_min, istringstream* fm_stream=NULL, vector<int> bac = vector<int>());

  Plan(const Plan &p);

	virtual ~Plan() {};

	void newCopy(Plan& p);

	// Adds a new station to the plan
	void add_station(Station& s);


	 // Generate apertures from the intensity matrices
	void generateApertures();

	const list<Station*>& get_stations() const ;

	void write_open_beamlets();

  double openBeamlet(int station, int aperture, int beamlet,
		bool delta_eval=false, list<pair<int, double> >* diff=NULL);

  double closeBeamlet(int station, int aperture, int beamlet, int side,
		 bool delta_eval=false, list<pair<int, double> >* diff=NULL);

  double openRow(int station, int aperture, int row, bool side,
		 bool delta_eval=false, list<pair<int, double> >* diff=NULL);
  
  list<pair<int, double> > getOpenRowDiff(int station, int aperture, int row, bool side);

  double closeRow(int station, int aperture, int row, bool side,
		 bool delta_eval=false, list<pair<int, double> >* diff=NULL);
  
  list<pair<int, double> > getCloseRowDiff(int station, int aperture, int row, bool side);

  double modifyIntensityAperture (int station, int aperture, int delta,
		bool delta_eval=false, list<pair<int, double> >* diff=NULL);
  
  list<pair<int, double> > getModifyIntensityDiff(int station, int aperture, int delta);

	void undoLast();

  void clearLast();

  //Lepi's version
  void undoLast2();

  void printSolution ();

  string toStringApertures ();

  string toStringIntensities ();

  void printIntensity(int station, bool vector_form=false);

  void writeIntensities(string file, int n);



	Station* get_station(int i) const{
	    list<Station*>::const_iterator s= stations.begin();
	    advance(s,i);
	    return *s;
	}

  int getNStations () {
     return(stations.size()) ;
  }

  Evaluator& get_evaluator(){return ev;}


  map<int, Station*> angle2station;
private:
	//The list of stations
	//list<Station> real_stations;
	list<Station*> stations;

  int n_stations;

  Station* last_changed;
	list< pair< int, double > > last_diff;



	Evaluator& ev;
	
};

} /* namespace imrt */

#endif /* PLAN_H_ */
