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
#include "ApertureILS.h"
#include "Matrix.h"
#include <math.h>
#include <float.h>

using namespace std;

namespace imrt {

class ACO {
private:
  std::clock_t time_begin;
  

public:
  ApertureILS * ils;
  int n_ants;
  
  double initial_pheromone;
  vector <map<int, Matrix> > pheromone;
  vector <map<int, Matrix> > heuristic;
  vector <map<int, Matrix>> probability;
  vector <Matrix> ipheromone;
  vector <Matrix> iheuristic;
  vector <Matrix> iprobability;
  
  double alpha;
  double beta;
  double rho;
  
  vector < Ant* > ants;
  Ant* iteration_best;
  Ant * global_best;
  
  Collimator & collimator;
  int max_apertures;
  double max_no_improvement = 100;
  

  ACO (vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
       vector<Volume>& volumes, int _max_apertures, int max_intensity, int initial_intensity,
       int step_intensity, int _n_ants, double _initial_pheromone, double _alpha, double _beta, 
       double _rho): collimator(_collimator), max_apertures(_max_apertures), n_ants(_n_ants), 
        initial_pheromone(_initial_pheromone), alpha(_alpha), beta(_beta), rho(_rho) {
  
    Ant * aux;
    double initial_best = DBL_MAX;
    
    initializePheromone();
    initializeHeuristic(w,Zmin, Zmax, volumes, max_intensity, initial_intensity, step_intensity);
    
    for (int i=0; i < n_ants; i++) {
      aux = new Ant(w, Zmin, Zmax, collimator, volumes, max_apertures, max_intensity, initial_intensity, step_intensity);
      ants.push_back(aux);
      if (initial_best > aux->getEvaluation()) {
        initial_best = aux->getEvaluation();
        iteration_best = aux;
      }
    }
    evaporation();
    globalDeposit (iteration_best);
    global_best =new Ant (*iteration_best);
    
    
  };

  virtual ~ACO() {
    for (int i=0; i < n_ants; i++) {
       delete ants[i];
    }
  };
  
  void initializeHeuristic (vector<double> w, vector<double> Zmin, vector<double> Zmax, vector<Volume>& volumes, 
                            int max_intensity, int initial_intensity, int step_intensity) {
    int nls=5;
    //Initialize structure
    for (int s=0; s< collimator.getNbAngles();s++) {
      map<int, Matrix> mheuristic;
      iheuristic.push_back(Matrix(max_apertures,  collimator.getIntensityLevelSize()));
      for (int a=0; a< max_apertures; a++) {
        for (int i=0; i< collimator.getIntensityLevelSize(); i++){
          iheuristic[s](a,i) = 0.5;
        }
        mheuristic[a] = Matrix(collimator.getXdim(), collimator.getReferenceSize());
        for (int j = 0; j < collimator.getXdim(); j++) {
          for (int r = 0; r < collimator.getReferenceSize(); r++) {
            mheuristic[a](j,r) =0.5;
          }
        }
      }
      heuristic.push_back(mheuristic);
    }
    
    //Perform local search and generate heuristic information
    Station *st;
    int aux;
    double deposit;
    ils = new ApertureILS(1, 1, true, true, 1, step_intensity, 100, 0.5, false, 0, ILS::ACCEPT_NONE, ApertureILS::FIRST_IMPROVEMENT);
    cout << "Initializing heuristic information with " << nls << " ls executions"<< endl;
    for (int i=0 ;i < nls; i++) {
      Plan *p = new Plan(w, Zmin, Zmax, collimator, volumes, max_apertures, max_intensity, initial_intensity, step_intensity, -1, Station::OPEN_MIN_SETUP);
      ils->simpleLocalSearch (*p, false);
      cout << " Local search "<<i<<": "<< p->getEvaluation() << endl;
      
      // Place deposit in heuristic information
      deposit = 5.0/p->getEvaluation();
      for (int s = 0; s < collimator.getNbAngles(); s++) {
        st = p->get_station(s);
        for (int a = 0; a < max_apertures; a++) {
          // Intensity heuristic 
          aux = collimator.getIntensityLevel(st->getApertureIntensity(a));
          iheuristic[s](a,aux) = iheuristic[s](a,aux) + deposit;
          // Aperture heuristic
          for (int j = 0; j < collimator.getXdim(); j++) {
            aux = collimator.searchReferenceIndex(st->getApertureShape(a, j));
            heuristic[s][a](j,aux) = heuristic[s][a](j,aux) + deposit;
          }
        }
      }
    }
  };
    
  void initializePheromone () {
    pair <int,int> aux, raux;
    for (int s = 0; s < collimator.getNbAngles(); s++) {
      map<int, Matrix> mpheromone;
      map<int, Matrix> mprobability;
      ipheromone.push_back(Matrix(max_apertures,  collimator.getIntensityLevelSize()));
      iprobability.push_back(Matrix(max_apertures,  collimator.getIntensityLevelSize()));
      for (int a=0; a< max_apertures; a++) {
        // Intensity initial pheromones
        for (int i=0; i< collimator.getIntensityLevelSize(); i++){
          ipheromone[s](a,i) = initial_pheromone;
        }
        // Aperture initial pheromones
        mpheromone[a] = Matrix(collimator.getXdim(), collimator.getReferenceSize());
        mprobability[a] = Matrix(collimator.getXdim(), collimator.getReferenceSize());
        for (int j = 0; j < collimator.getXdim(); j++) {
          aux = collimator.getActiveRange(j, collimator.getAngle(s));
          for (int r = 0; r < collimator.getReferenceSize(); r++) {
            raux = collimator.getReference(r);
            if (raux.first >= aux.first &&
                raux.second <= aux.second) {
               mpheromone[a](j,r) = initial_pheromone;
            } else {
               mpheromone[a](j,r) = 0;
            }
          }
        }
      }
      pheromone.push_back(mpheromone);
      probability.push_back(mprobability);
    }
  };
  
  void printPheromone () {
    cout << "Pheromone matrix: " << endl;
    for (int s = 0; s < collimator.getNbAngles(); s++) {
      cout << "Angle " << s << endl;
      for (int a=0; a<max_apertures; a++) {
        cout << "Aperture " << a << endl;
        for (int j = 0; j < collimator.getXdim(); j++) {
          for (int r = 0; r <  collimator.getReferenceSize(); r++) {
            cout << pheromone[s][a](j,r) << " ";
          }
          cout << endl;
        }
      }
    }
  };
  
  void printiPheromone () {
    cout << "intensity pheromone: " <<endl;
    for (int s=0; s<collimator.getNbAngles(); s++) {
      cout << "station "<< s<< endl;
      for (int a=0; a < max_apertures; a++) {
        for (int i=0; i <  collimator.getIntensityLevelSize();i++)
          cout << ipheromone[s](a,i) << " ";
        cout << endl;
      }
    }
  }
  
  void printHeuristic () {
    cout << "Heuristic matrix: " << endl;
    for (int s = 0; s < collimator.getNbAngles(); s++) {
      cout << "Angle " << s << endl;
      for (int a=0; a<max_apertures; a++) {
        cout << "Aperture " << a << endl;
        for (int j = 0; j < collimator.getXdim(); j++) {
          for (int r = 0; r <  collimator.getReferenceSize(); r++) {
            cout << heuristic[s][a](j,r) << " ";
          }
          cout << endl;
        }
      }
    }
  };
  
  void printiHeuristic () {
    cout << "intensity heuristic: " <<endl;
    for (int s=0; s<collimator.getNbAngles(); s++) {
      cout << "station "<< s<< endl;
      for (int a=0; a < max_apertures; a++) {
        for (int i=0; i <  collimator.getIntensityLevelSize();i++)
          cout << iheuristic[s](a,i) << " ";
          cout << endl;
      }
    }
  }
  
  
  void printAnts(){
    cout << "Printing all ants intensity..." << endl;
    for (int ant = 0; ant < n_ants; ant++) {
      cout << "ant # " << ant << endl;
        ants[ant]->printIntensities();
    }
  };
  
  void printAntsReference(){
    cout << "Printing all ants reference ..." << endl;
    for (int ant = 0; ant < n_ants; ant++) {
      cout << "ant # " << ant << endl;
      for (int i = 0; i < collimator.getNbAngles(); i++) {
        ants[ant]->printReferencePlan();
      }
    }
  };
  
  Ant* getBestAnt (){
    Ant * best;
    double best_eval= DBL_MAX;
    for (int a =0; a< n_ants;a++){
      if (best_eval< ants[a]->getEvaluation()) {
        best_eval = ants[a]->getEvaluation();
        best = ants[a];
      }
      return(best);
    }
    
  };
  
  void evaporation () {
    for (int s =0 ; s< collimator.getNbAngles();s++) 
      for (int a=0; a<max_apertures;a++) {
        pheromone[s][a] =  pheromone[s][a] *(1.0-rho);
        ipheromone[s] = ipheromone[s] *(1.0-rho);
      }
  };
  
  void globalDeposit (Ant * a) {
    int raux;
    double deposit = 1/a->getEvaluation();
    for (int s = 0; s < collimator.getNbAngles(); s++) {
      for (int ap = 0; ap < a->n_apertures; ap++) {
        // Intensity pheromone deposit
        for (int i =0; i <  collimator.getIntensityLevelSize(); i++) {
          raux = a->getIntensityReferenceAt(s,ap);
          ipheromone[s](ap,raux) = ipheromone[s](ap,raux) + deposit;
        }
        // Aperture pheromone deposit
        for (int j = 0; j < collimator.getXdim(); j++) {
          raux = a->getReferenceAt(s,ap,j);
          (pheromone[s])[ap](j,raux) = (pheromone[s])[ap](j,raux) + deposit;
        }
      }
    }
  };
  
  virtual void generateTours () = 0;
  
  virtual void calculateProbability ()=0;
  
  virtual void search (int max_eval) = 0;
  
};

}

#endif /* ACO_H_ */
