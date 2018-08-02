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
  
  double search(Plan& P, int max_iterations) {
    
    cout << "Staring ILS search." << endl;
    //Plan best_plan=Plan(P);
    pair<bool, pair<Station*, int>> target_beam;
    double local_eval, aux_eval,  best_eval=P.eval();
    int no_improvement;
    
    no_improvement = 0;
    for (int s=0;s<max_iterations;s++) {
      target_beam = getLSBeamlet(P);
      
      cout << "Iteration " << (s+1) << ", best: " << best_eval << ", beamlet: " << target_beam.second.second  << ", station: " << target_beam.second.first->getAngle() << ", +-: " << target_beam.first;
      aux_eval = localSearch (target_beam, P);
      cout << endl;
      
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
  
  
};

}

#endif /* ILS_H_ */
