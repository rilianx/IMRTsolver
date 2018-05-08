/*
 * Station.h
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include <vector>

#ifndef STATION_H_
#define STATION_H_

using namespace std;

namespace imrt {



class Station {
public:
	Station(int nb_beamrow, int nb_beamcol) : nb_beamrow(nb_beamrow), nb_beamcol(nb_beamcol),
	X(nb_beamrow), angle(0.0), beamlets(nb_beamrow*nb_beamcol){
		for(int i=0; i<nb_beamrow*nb_beamcol; i++)
			beamlets[i]=false;
	}

	virtual ~Station();

	/*
	 * The aperture is set, and the corresponding beamlets are open/closed
	 */
	void set_aperture(int i, int mid, int ext){
		X[i]= make_pair(mid,ext);

		int j= (int)((mid-ext)/2.0+0.5);
		for(int b=i*nb_beamcol ; b<(i+1)*nb_beamcol; b++){
			if(b>=i*nb_beamcol+j && b < i*nb_beamcol+j + ext){
				if(beamlets[b]==false){
					last_changes.push_back( make_pair(b,true) );
					beamlets[b]=true;
				}
			}else if(beamlets[b]==true){
				last_changes.push_back( make_pair(b,false) );
				beamlets[b]=false;
			}
		}
	}

	double get_angle() {return angle;}

	/* save the last changes (turn on/off) of the beamlets for incremental evaluation */
	list < pair<int, bool> > last_changes;

	vector<bool> beamlets;

private:
	double angle;
	vector<pair <int,int> > X; // (e.g., 2*midpoint and extension of the aperture)

	int nb_beamrow;
	int nb_beamcol;





};

} /* namespace imrt */

#endif /* STATION_H_ */
