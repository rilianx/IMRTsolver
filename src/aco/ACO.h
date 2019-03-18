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
  map < int, pair <int,int> > reference;
  map <int, Matrix> probability;
  int ref_size;
  
  double alpha;
  double beta;
  
  vector < Plan* > ants;
  Plan * best;
  
  
  Collimator & collimator;
  int max_apertures;
  double max_no_improvement = 100;
  
  static const int ACCEPT_NONE = 0;
  static const int ACCEPT_SA = 1;

  ACO (vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
       vector<Volume>& volumes, int _max_apertures, int max_intensity, int initial_intensity,
       int step_intensity, int _n_ants, double _initial_pheromone, double _alpha, double _beta): collimator(_collimator), max_apertures(_max_apertures), n_ants(_n_ants), initial_pheromone(_initial_pheromone), alpha(_alpha), beta(_beta) {
  
    generateReference ();
    printReference ();
    initializePheromone();
    printPheromone ();
    
    best = new Plan(w, Zmin, Zmax, collimator, volumes, max_apertures, max_intensity, initial_intensity, step_intensity, -1, Station::OPEN_MIN_SETUP);
    for (int i=0; i < n_ants; i++) {
      ants.push_back(new Plan(w, Zmin, Zmax, collimator, volumes, max_apertures, max_intensity, initial_intensity, step_intensity, -1, Station::CLOSED_MIN_SETUP));
    }
    
  };

  virtual ~ACO() { };
  
  void generateReference () {
    int i = 0;
    reference[0] = make_pair(-1,-1);
    i = 1;
    for (int s = 1; s < collimator.getYdim(); s++) {
      for (int j = 0; (j+s) < collimator.getYdim(); j++) {
        reference[i] = make_pair(j,j+s);
        i++;
      }
    }
    reference[i] = make_pair(0,collimator.getYdim()-1);
    ref_size = i;
  };
  
  void printReference () {
    cout << "Reference list" << endl;
    for (int i=0; i < ref_size; i++) {
      cout << i << " " << reference[i].first << "," << reference[i].second << endl;
    }
  };
    
  void initializePheromone () {
    
    pair <int,int> aux;
    for (int i = 0; i < collimator.getNbAngles(); i++) {
      pheromone[i] = Matrix(collimator.getXdim(), ref_size);
      probability[i] = Matrix(collimator.getXdim(), ref_size);
      for (int j = 0; j < collimator.getXdim(); j++) {
        aux = collimator.getActiveRange(j, collimator.getAngle(i));
        for (int s = 0; s <  ref_size; s++) {
          if (reference[s].first >= aux.first &&
              reference[s].second <= aux.second) {
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
        for (int s = 0; s <  ref_size; s++) {
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
