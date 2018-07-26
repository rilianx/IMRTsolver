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
	Plan(EvaluationFunction &ev) : ev(ev) {};
	virtual ~Plan() {};

	// Adds a new station to the plan
	void add_station(Station& s){
		stations.push_back(&s);
	}

	double eval(vector<double>& w, vector<double>& Zmin, vector<double>& Zmax){
		double eval=ev.eval(*this,w,Zmin,Zmax);
		ev.generate_voxel_dose_functions ();
		return eval;
	}

	const list<Station*>& get_stations() const{
		return stations;
	}

	void print_open_beamlets(){
		ofstream myfile;
		myfile.open ("openbeamlets.txt");

		myfile << "Angles\t";
		for(auto s:stations)
			myfile << s->getAngle() << "\t";
		myfile << endl;
		int k=0;
		for(auto s:stations){
			set<int> open_beamlets;
			myfile << endl << "Station Angle\t" << s->getAngle() << endl;
			for(int i=0; i<s->getNbApertures(); i++){
				myfile << "Aperture\t" << i << endl;
				myfile << "Intensity\t" << s->intensity[i] << endl;

				myfile << "OpenBeamlets\t" ;
				bool first=true;
				for(auto beam:s->open_beamlets(i)){
					if(!first) myfile << "\t";
					first=false;
					myfile << k+beam;
				}
				myfile << endl;
			}
			k+=s->getNbBeamlets();
		}
		myfile.close();

	}

private:
	//The list of stations
	list<Station*> stations;

	EvaluationFunction& ev;
};

} /* namespace imrt */

#endif /* PLAN_H_ */
