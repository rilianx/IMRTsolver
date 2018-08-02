/*
 * ils.h
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#ifndef ILS_H_
#define ILS_H_

#include "Plan.h"

namespace imrt {

class ILS {
private:
  double max_no_improvement=100;

public:
  int bsize;
  int vsize;
  ILS(int bsize, int vsize): bsize(bsize), vsize(vsize) {};
  
  virtual double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P) = 0;
  virtual bool acceptanceCriterion(double new_eval, double prev_eval)=0;

  virtual pair<bool, pair<Station*, int>> getLSBeamlet(Plan& P) =0;
  //virtual double perturbation(Plan& P)=0;
  //virtual bool perturbate(int no_improvement)=0;
  
  double search(Plan& P, EvaluationFunction& F, int max_iterations) {
    Plan best_plan=Plan(P);
    pair<bool, pair<Station*, int>> target_beam;
    double local_eval, best_eval, aux_eval;
    int no_improvement;
    
    no_improvement = 0;
    for (int s=0;s<max_iterations;s++) {
      target_beam = getLSBeamlet(P);
      aux_eval = localSearch (target_beam, P);
      
      if (aux_eval < best_eval) {
        best_eval=aux_eval;
        //best_plan=Plan(P);
      }
      
      if (aux_eval < local_eval) {
        local_eval = aux_eval;
        no_improvement = 0;
      } else if (acceptanceCriterion(aux_eval, local_eval)) {
        local_eval = aux_eval;
        no_improvement = 0;
      } else {
        P.undoLast();
        no_improvement ++;
      }
      
      //if ( perturbate(no_improvement))
      //  local_eval = perturbation(P);
    }
    
    //P=Plan(best_plan);
    return(best_eval);
  };
  
  void printHeader() {
    cout << "************************************************"<< endl;
    cout << "************************************************"<< endl;
    cout << "****************** IMRT-Solver  ****************"<< endl;
    cout << "Iterations: " << maxiter << endl;
    cout << "Seed: " << seed << endl;
    cout << "Temperature: " << temperature << endl;
    if (ls_apertures)
      cout << "Searching: aperture pattern" << endl;
    if (ls_intensity)
      cout << "Searching: intensity" << endl;
    if (ls_both){
      cout << "Searching: intensity and aperture pattern" << endl;
      cout << "Probability intensity ls: " << prob_intensity << endl;
    }
    
    cout << endl << "Colimator configuration: "<< endl;
    cout << "  Stations: " << stations.size() << endl;
    cout << "  Angles: ";
    for (int i=0; i<stations.size();i++) cout << stations[i]->getAngle() << " ";
    cout << endl;
    cout << "  Max apertures: " << max_apertures << endl;
    cout << "  Initial intensity: " << initial_intensity << endl;
    cout << "  Open initial setup: " << open_setup << endl;
    
    
    cout << endl << "Instance information: "<< endl;
    cout << "  Volumes: " << volumes.size() << endl;
    
    cout << "************************************************"<< endl<< endl;
    cout << "************************************************"<< endl;
    cout << "************************************************"<< endl;
  }
  
};

}

#endif /* ILS_H_ */
