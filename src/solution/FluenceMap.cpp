/*
 * FluenceMap.cpp
 *
 *  Created on: 1 june 2021
 *      Author: iaraya
 */

#include "Evaluator.h"
#include <string>
#include <functional>

namespace imrt {
	FluenceMap::FluenceMap(vector<Volume>& volumes, const Collimator& collimator) :
	volumes(volumes), nb_organs(volumes.size()), nb_voxels(volumes.size()) {
	    for(int i=0; i<nb_organs; i++){
		    nb_voxels[i]=volumes[i].getNbVoxels();
		    FM.insert(FM.end(), vector<double>(nb_voxels[i]));
	    }

        deltaFM.resize(nb_organs);

        create_voxel2beamlet_list(volumes, collimator);

    }
    
	FluenceMap::~FluenceMap() { }

	void FluenceMap::create_voxel2beamlet_list(vector<Volume>& volumes, const Collimator& collimator){
	    cout << "creating maps: voxel2beamlet and beamlet2voxel" << endl;
	    for (int i=0;i<collimator.getNbAngles();i++){
            for(int o=0; o<nb_organs; o++){
                int angle = collimator.getAngle(i);
                if( !volumes[o].valid_angle(angle)) continue;
                
                const Matrix&  D = volumes[o].getDepositionMatrix(angle);
                for(int k=0; k<nb_voxels[o]; k++){
                    for(int b=0; b<collimator.getNangleBeamlets(angle); b++)
                        if(D(k,b) > 0.0) {
                            voxel2beamlet_list[angle][make_pair(o,k)].push_back(b);
                            beamlet2voxel_list[angle][b].insert(make_pair(-abs(D(k,b)), make_pair(o,k)));
                        }
                }
            }
	    }
    }

	// Generate the dose distribution matrices Z for each organ
	void FluenceMap::computeFM(const Plan& p){
        const list<Station*>& stations=p.get_stations();

	    for(int o=0; o<nb_organs; o++)
	 	    std::fill(FM[o].begin(), FM[o].end(), 0.0);

        //we update the dose distribution matrices Z with the dose delivered by the station
        for(int o=0; o<nb_organs; o++){
            for(int k=0; k<nb_voxels[o]; k++){
                double dose=0.0;
                for (auto s: p.get_stations()){
                    int angle = s->getAngle();
                    list<int>& beamlets = voxel2beamlet_list[angle][make_pair(o,k)];
                    for(auto b:beamlets){
                        const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
                        dose += Dep(k,b)*p.angle2station.at(angle)->getIntensity(b);
                    }
                }
                FM[o][k] += dose;
            }
        }
    }

	//Get deltaZ related to the list of changes
	std::vector<list < pair<int, double>>>& FluenceMap::compute_deltaFM(list< pair< int, double > >& changes, double angle) const{
        for(int o=0; o<nb_organs; o++){
            deltaFM[o].clear();
            const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
            for(int k=0; k<nb_voxels[o]; k++){
                double delta = 0;
                for (auto let:changes){
                    int b=let.first;
                    if(Dep(k,b)==0.0) continue;
                    delta+= Dep(k,b)*let.second;
                }
                if(delta==0.0) continue;
                deltaFM[o].push_back(make_pair(k,delta));	
            }
        }

        return deltaFM;
    }

	void FluenceMap::updateFM(list< pair< int, double > >& changes, double angle, std::vector<list < pair<int, double>>>& deltaFM){
        for(int o=0; o<deltaFM.size(); o++){
            for(auto kk:deltaFM[o]){
                int k=kk.first;
                double delta = kk.second;
		        FM[o][k] += delta;
	        }
	    }
    }
}