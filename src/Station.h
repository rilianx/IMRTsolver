/*
 * Station.h
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include <vector>
#include <set>
#include <list>
#include <iostream>

#ifndef STATION_H_
#define STATION_H_

using namespace std;

namespace imrt {


/* An IMRT station
 *
 * A station consists in a beam (matrix of beamlets), an angle, an aperture
 * and an intensisty
 */

class Station {
public:


	/* The Station constructor
	 *
	 * @nb_beamrow number of rows of the beamlet matrix
	 * @nb_beamcol number of columns of the beamlet matrix
	 */
	Station(int nb_beamrow, int nb_beamcol) : nb_beamrow(nb_beamrow), nb_beamcol(nb_beamcol),
	X(nb_beamrow), angle(0.0), beamlets(nb_beamrow*nb_beamcol), intensity(1){
		for(int i=0; i<nb_beamrow*nb_beamcol; i++)
			beamlets[i]=false;

		for(int i=0; i< X.size(); i++)
			X[i]= make_pair(nb_beamcol/2,nb_beamcol/2);
	}

	virtual ~Station();

	// The aperture is set, and the corresponding beamlets are open/closed
	bool set_aperture(int i, int ini, int fin){
		if(ini>nb_beamcol/2 || fin<nb_beamcol/2) return false;
		if(ini==X[i].first && fin==X[i].second) return true;

		X[i]= make_pair(ini,fin);

		for(int b=i*nb_beamcol ; b<(i+1)*nb_beamcol; b++){
			if(b>=i*nb_beamcol+ini && b < i*nb_beamcol+fin){
				if(beamlets[b]==false) switch_beamlet(b);
			}else
				if(beamlets[b]==true) switch_beamlet(b);
		}

		return true;
	}

	// Close all the beamlets
	void close_aperture(){
		for(int i=0; i<nb_beamrow; i++)
			set_aperture(i, nb_beamcol/2,nb_beamcol/2);
	}

	// Turn on/off the beamlet b
	void switch_beamlet(int b){
		if(switched_lets.find(b)==switched_lets.end())
			switched_lets.insert(b);
		else
			switched_lets.erase(b);

		beamlets[b] = !beamlets[b];
	}

	void set_angle(double ang) {angle=ang;}

	double get_angle() const {return angle;}

	void set_intensity(int intens) {intensity=intens;}

	int get_intensity() const {return intensity;}

	// set of beamlets switched from the last evaluation
	set < int > switched_lets;

	// the beamlets in a vector representation required for the evaluation function
	vector<bool> beamlets;


private:
	double angle;

	//a beamlet row is represented by two integer values delimiting its open beamlets
	vector<pair <int,int> > X;

	int intensity;

	int nb_beamrow;
	int nb_beamcol;


};

} /* namespace imrt */

#endif /* STATION_H_ */
