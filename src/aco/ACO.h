/*
 *  ACO.h
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#ifndef ACO_H_
#define ACO_H_

#include <ctime>
#include <algorithm> 
#include "Plan.h"
#include "Ant.h"
#include "Matrix.h"
#include <math.h>

namespace imrt {

class ACO {
private:
  std::clock_t time_begin;

public:
  int n_ants;
  
  double initial_pheromone;
  map <int, Matrix > pheromone;
  map <int, Matrix> probability;
  
  double alpha;
  double beta;
  
  vector < Ant* > ants;
  Ant * best;
  
  Collimator & collimator;
  int max_apertures;
  double max_no_improvement = 100;
  
  static const int ACCEPT_NONE = 0;
  static const int ACCEPT_SA = 1;

  ACO (vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
       vector<Volume>& volumes, int _max_apertures, int max_intensity, int initial_intensity,
       int step_intensity, int _n_ants, double _initial_pheromone, double _alpha, double _beta): collimator(_collimator), max_apertures(_max_apertures), n_ants(_n_ants), initial_pheromone(_initial_pheromone), alpha(_alpha), beta(_beta) {
  
    Ant * aux;
    initializePheromone();
    printPheromone ();
    
    for (int i=0; i < n_ants; i++) {
      aux = new Ant(w, Zmin, Zmax, collimator, volumes, max_apertures, max_intensity, initial_intensity, step_intensity);
      ants.push_back(aux);
    }
    
  };

  virtual ~ACO() { };
    
  void initializePheromone () {
    
    pair <int,int> aux, raux;
    for (int i = 0; i < collimator.getNbAngles(); i++) {
      pheromone[i] = Matrix(collimator.getXdim(), collimator.getReferenceSize());
      probability[i] = Matrix(collimator.getXdim(), collimator.getReferenceSize());
      for (int j = 0; j < collimator.getXdim(); j++) {
        aux = collimator.getActiveRange(j, collimator.getAngle(i));
        for (int s = 0; s < collimator.getReferenceSize(); s++) {
          raux = collimator.getReference(s);
          if (raux.first >= aux.first &&
              raux.second <= aux.second) {
             pheromone[i](j,s) = initial_pheromone;
          } else {
             pheromone[i](j,s) = 0;
          }
        }
      }
    }
  };
  
  void printPheromone () {
    cout << "Pheromone matrix" << endl;
    for (int i = 0; i < collimator.getNbAngles(); i++) {
      cout << "Angle " << i << endl;
      for (int j = 0; j < collimator.getXdim(); j++) {
        for (int s = 0; s <  collimator.getReferenceSize(); s++) {
          cout << pheromone[i](j,s) << " ";
        }
        cout << endl;
      }
    }
  };
  
  void printAnts(){
    cout << "ALL ANTS" << endl;
    for (int ant = 0; ant < n_ants; ant++) {
      cout << "ANT " << ant << endl;
      for (int i = 0; i < collimator.getNbAngles(); i++) {
        ants[ant]->printIntensity(i);
      }
    }
  };
  
  virtual void generateTours () = 0;
  
  virtual void calculateProbability ()=0;
  
};

}

#endif /* ACO_H_ */
