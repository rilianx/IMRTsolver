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
  
  double search(Plan& current_plan, int max_iterations) {
    
    cout << "Staring ILS search." << endl;
    Plan best_plan (current_plan);
    pair<bool, pair<Station*, int>> target_beam;
    double local_eval, aux_eval,  best_eval=current_plan.eval();
    int no_improvement;
    no_improvement = 0;
    for (int s=0;s<max_iterations;s++) {
      //cout << "ss"<< endl;
      target_beam = getLSBeamlet(current_plan);
      if (target_beam.second.second<0) 
        cout << "ERROR: No beamlet available" << endl;
      
      cout << "Iteration " << (s+1) << ", best: " << best_eval << ", beamlet: " << target_beam.second.second  << 
              ", station: " << target_beam.second.first->getAngle() << ", +-: " << target_beam.first;
      aux_eval = localSearch (target_beam, current_plan);
      cout << endl;
      
      if (aux_eval < best_eval) {
        best_eval=aux_eval;
        best_plan.newCopy(current_plan);
      }
      
      if (aux_eval < local_eval) {
        local_eval = aux_eval;
        no_improvement = 0;
      } else if (aux_eval!= local_eval && acceptanceCriterion(aux_eval, local_eval)) {
        local_eval = aux_eval;
        no_improvement = 0;
      } else {
        current_plan.undoLast();
        no_improvement ++;
      }
      //if ( perturbate(no_improvement))
      //  local_eval = perturbation(P);
    }
    
    aux_eval=best_plan.eval();
    cout << "comparison:" << aux_eval << " vs "<< best_eval<<endl;
    return(best_eval);
  };
  
  
};

}

#endif /* ILS_H_ */
