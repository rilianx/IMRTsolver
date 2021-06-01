
/*
 * EvaluationFunction.cpp
 *
 *  Created on: 31 may. 2021
 *      Author: iaraya
 */

#include "EvaluatorF.h"
#include <string>
#include <functional>

namespace imrt {

EvaluatorF* EvaluatorF::instance = NULL;

int EvaluatorF::n_evaluations=0;

double EvaluatorF::eval(const Plan& p){
	voxels.clear();
    cout << 1 << endl;
	z_structure.generate_Z(p);
	F=0.0;
	for(int o=0; o<Z.size(); o++){
		double pen=0.0;
		for(int k=0; k<Z[o].size(); k++){
			if(Z[o][k] < Zmin[o] )
				 pen += w[o] * ( pow(Zmin[o]-Z[o][k], 2) );

			if(Z[o][k] > Zmax[o] )
				 pen += w[o] * ( pow(Z[o][k]-Zmax[o], 2) );

			update_sorted_voxels(o, k);
		}
		F+=pen/Z[o].size();
	}
    
	n_evaluations++;
	return F;
}

double EvaluatorF::get_delta_eval(list< pair< int, double > >& changes, double angle) const{

  std::vector<list < pair<int, double>>>& deltaZ = z_structure.compute_deltaZ(changes,angle);

  double delta_F=0.0;

  for(int o=0; o<deltaZ.size(); o++){
      double pen=0.0;
    for(auto kk:deltaZ[o]){
        int k=kk.first;
        double delta = kk.second;
        //with the change in the dose of a voxel we can incrementally modify the value of F
        if(Z[o][k] < Zmin[o] && Z[o][k] + delta < Zmin[o]) //update the penalty
            pen += w[o]*delta*(delta+2*(Z[o][k]-Zmin[o]));
        else if(Z[o][k] < Zmin[o]) //the penalty disappears
            pen -=  w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
        else if(Z[o][k] + delta < Zmin[o]) //the penalty appears
            pen +=  w[o] * ( pow(Zmin[o]-(Z[o][k]+delta), 2) );

        if(Z[o][k] > Zmax[o] && Z[o][k] + delta > Zmax[o]) //update the penalty
            pen += w[o]*delta*(delta+2*(-Zmax[o] + Z[o][k]));
        else if(Z[o][k] > Zmax[o]) //the penalty disappears
            pen -=  w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
        else if(Z[o][k] + delta > Zmax[o]) //the penalty appears
            pen +=  w[o] * ( pow(Z[o][k]+delta - Zmax[o], 2) );

    }
    delta_F += pen/(double) Z[o].size();
  }

  return delta_F;
}


double EvaluatorF::incremental_eval(list< pair< int, double > >& changes, double angle){

	
	double deltaF = get_delta_eval(changes, angle);
    std::vector<list < pair<int, double>>>& deltaZ = z_structure.get_deltaZ();

    z_structure.apply_deltaZ(changes, angle, deltaZ);

    for(int o=0; o<deltaZ.size(); o++){
    for(auto kk:deltaZ[o]){
        int k=kk.first;
        double delta = kk.second;
		update_sorted_voxels(o, k);
	  }
    }

    
    F+=deltaF;
    //cout << F << endl;
    return F;
}

//Return the beamlets sorted by impact on F taking into account the nv worst voxels.
//Each returned beamlet is a pair eval,(station, beamlet)
multimap < double, pair<int, int>, MagnitudeCompare2 > EvaluatorF::sorted_beamlets(const Plan& p, double vsize){
	int tot_voxels=0;
	for(int o=0; o<Z.size(); o++) tot_voxels+=Z[o].size();

	int nv = (double) tot_voxels*vsize;

	multimap < double, pair<int, int>, MagnitudeCompare2 > bestb;
	double max_ev=0.0;
	int k=0;
	for(auto s:p.get_stations()){
		for(int b=0; b<s->getNbBeamlets(); b++){
			double ev=0; int i=0;

            for (auto voxel : voxels){
				  double d = voxel.first;
				  int o = voxel.second.first;
				  int k= voxel.second.second;

				  const Matrix&  Dep = s->getDepositionMatrix(o);
				  ev += d * Dep(k,b);
				  i++; if(i==nv) break;
			 }

			 bestb.insert(make_pair( ev, make_pair(k,b)));

		}
		k++;
	}
	return bestb;
    }


    void EvaluatorF::update_sorted_voxels(int o, int k){
		voxels.erase(make_pair(D[o][k],make_pair(o,k)));
        D[o][k] = 0.0;

		if(Zmin[o]>0 && Z[o][k] < Zmin[o]) D[o][k]=-w[o]*(Z[o][k]-Zmin[o])/Z[o].size();

		else if(Z[o][k] >= Zmax[o]) D[o][k]=w[o]*(Z[o][k]-Zmax[o])/Z[o].size();

		if(D[o][k]!=0.0)
			voxels.insert(make_pair(D[o][k],make_pair(o,k)));
    }

}