/*
 * ACS.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "ACS.h"

namespace imrt {


  ACS::ACS(vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
           vector<Volume>& volumes, int max_apertures, int max_intensity, int initial_intensity,
           int step_intensity, int _n_ants, double _initial_pheromone, double _alpha, double _beta):
           ACO(w, Zmin, Zmax, _collimator, volumes, max_apertures, max_intensity, initial_intensity,
               step_intensity, _n_ants, _initial_pheromone, _alpha, _beta){
//  cout << "per:" << perturbation_size << endl;
             q0=0.5;
             calculateProbability();
             printProbability();
           };
  
  void ACS::generateTours() {
    double aux;
    bool flag;
    for (int ant=0; ant < n_ants; ant++) {
      for (int i = 0; i < collimator.getNbAngles(); i++) {
        for (int j = 0; j < collimator.getXdim(); j++) {
          for (int ap=0; ap < max_apertures; ap++) {
            aux = rand() / double(RAND_MAX);
            flag = false;
            for (int r = 0; r < ref_size && !flag; r++) {
              if (probability[i](j,r) >= aux) {
                ants[ant]->get_station(i)->setApertureShape(ap,j, reference[r].first, reference[r].second);
                flag = true;
              }
            }
          }
        }
        ants[ant]->get_station(i)->generateIntensity();
      }
      
    }
  };
  
  void ACS::calculateProbability (){
    double sum, aux;
    for (int i = 0; i < collimator.getNbAngles(); i++) {
      for (int j = 0; j < collimator.getXdim(); j++) {
        sum = 0;
        for (int s = 0; s < ref_size; s++) {
          sum = sum + pow(pheromone[i](j,s),alpha);
        }
        aux = 0;
        for (int s = 0; s < ref_size; s++) {
          probability[i](j,s) =  pow(pheromone[i](j,s),alpha);
          probability[i](j,s) = aux + pow(pheromone[i](j,s),alpha) / sum;;
          aux  = probability[i](j,s);
        }
      }
    }
    
  };
  
  void ACS::printProbability() {
    cout << "Probabilities" << endl;
    for (int i = 0; i < collimator.getNbAngles(); i++) {
      for (int j = 0; j < collimator.getXdim(); j++) {
        for (int s = 0; s < ref_size; s++) {
          cout << probability[i](j,s) <<  " ";
        }
        cout << endl;
      }
    }
    
  };


}



