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
           int step_intensity, int _n_ants, double _initial_pheromone, double _alpha, double _beta, double _rho):
           ACO(w, Zmin, Zmax, _collimator, volumes, max_apertures, max_intensity, initial_intensity,
               step_intensity, _n_ants, _initial_pheromone, _alpha, _beta, _rho){
             q0=0.5;
             calculateProbability();
             //printProbability();
           };
  
  void ACS::generateTours() {
    double best_eval=DBL_MAX;
    for (int a=0; a < n_ants; a++) {
      ants[a]->generateTour(probability, iprobability);
      if (ants[a]->getEvaluation() < best_eval) {
        best_eval = ants[a]->getEvaluation();
        iteration_best = ants[a];
      }
    }
  };

  
  void ACS::calculateProbability (){
    double sum, aux;
    for (int s = 0; s < collimator.getNbAngles(); s++) {
      for (int a = 0; a < max_apertures; a++ ){
        // Intensity probability
        sum = 0;
        for (int i=0; i< collimator.getIntensityLevelSize(); i++) {
          sum = sum + pow(ipheromone[s](a,i),alpha)*pow(iheuristic[s](a,i),beta);
        }
        aux=0;
        for (int i=0; i< collimator.getIntensityLevelSize(); i++) {
          iprobability[s](a,i) = aux + (pow(ipheromone[s](a,i),alpha)*pow(iheuristic[s](a,i),beta)/sum);
          aux = iprobability[s](a,i);
        }
        // Aperture probability
        for (int j = 0; j < collimator.getXdim(); j++) {
          sum = 0;
          for (int r = 0; r < collimator.getReferenceSize(); r++) {
            sum = sum + pow(pheromone[s][a](j,r),alpha)*pow(heuristic[s][a](j,r),beta);
          }
          aux = 0;
          for (int r = 0; r < collimator.getReferenceSize(); r++) {
            probability[s][a](j,r) = aux + pow(pheromone[s][a](j,r),alpha)*pow(heuristic[s][a](j,r),beta) / sum;;
            aux  = probability[s][a](j,r);
          }
        }
      }
    }
  };
  
  void ACS::printProbability() {
    cout << "Intensity probabilities:"<< endl;
    for (int s = 0; s < collimator.getNbAngles(); s++) {
      for (int a = 0; a < max_apertures; a++) {
          for (int i = 0; i <  collimator.getIntensityLevelSize();i++) {
            cout << iprobability[s](a,i) << " ";
        }
        cout << endl;
      }
    }
    cout << "Aperture probabilities:" << endl;
    for (int s = 0; s < collimator.getNbAngles(); s++) {
      for (int a = 0; a < max_apertures; a++) {
        for (int j = 0; j < collimator.getXdim(); j++) {
          for (int r = 0; r < collimator.getReferenceSize(); r++) {
            cout << probability[s][a](j,r) <<  " ";
          }
          cout << endl;
        }
      }
    }
  };
  
  
  void ACS::search(int max_eval) {
    int current_eval = 0;
    int iteration = 1;
    while ( (current_eval+ n_ants) <= max_eval) {
      generateTours();
      
      ils->simpleLocalSearch (*(iteration_best->getPlan()), false);
      iteration_best->generateReference();
      
      current_eval = current_eval + n_ants;
      if (global_best->getEvaluation() > iteration_best->getEvaluation())
        global_best->copy(*iteration_best);
      cout << "Iteration: " << iteration << " n_evals: "<< current_eval <<
          " iter_best: " << iteration_best->getEvaluation() <<
          " global_best: " << global_best->getEvaluation()<<endl;
      //iteration_best->printIntensities();
      evaporation();
      globalDeposit(global_best);
      iteration++;
      //printPheromone();
    }
     
  };
  


}



