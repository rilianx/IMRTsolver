/*
 *  ILS.h
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#ifndef ILS_H_
#define ILS_H_

#include <ctime>
#include "Plan.h"

namespace imrt {

class ILS {
private:
  double max_no_improvement=100;
  std::clock_t time_begin;

public:
  int bsize;
  int vsize;
  int acceptance;

  static const int ACCEPT_NONE = 0;
  static const int ACCEPT_SA = 1;

  ILS(int bsize, int vsize, int acceptance=ACCEPT_NONE): bsize(bsize), vsize(vsize), acceptance(acceptance) {
  };

  virtual double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P) = 0;
  virtual bool acceptanceCriterion(double new_eval, double prev_eval)=0;

  virtual pair<bool, pair<Station*, int>> getLSBeamlet(Plan& P){
	  return P.getLSBeamlet(bsize, vsize);
  }
  //virtual double perturbation(Plan& P)=0;
  //virtual bool perturbate(int no_improvement)=0;

  virtual void updateTemperature() {};

  double search(Plan& current_plan, int max_time, int max_iterations) {

    cout << "Staring ILS search." << endl;
    std::clock_t time_end;
    Plan best_plan (current_plan);
    pair<bool, pair<Station*, int>> target_beam;
    double local_eval, aux_eval,  best_eval=current_plan.eval();
    double used_time=0;
    bool flag=true;
    int no_improvement, iteration=1;
    no_improvement = 0;

    local_eval=best_eval;
    //Start time
    time_begin=clock();

    while (flag) {
      //cout << "ss"<< endl;
      target_beam = getLSBeamlet(current_plan);
      if (target_beam.second.second<0)
        cout << "ERROR: No beamlet available" << endl;

      cout << "Iteration: " << iteration << ", time: "<< used_time << ", best: " << best_eval <<
              ", current: " << local_eval  << ", beamlet: " << target_beam.second.second  <<
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
        current_plan.undoLast2();
        no_improvement ++;
      }


      if (acceptance==ACCEPT_SA)
         updateTemperature();
      iteration++;

      // Termination criterion
      time_end=clock();
      used_time=double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) flag=false;
      if (max_iterations!=0 && iteration>=max_iterations) flag=false;

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
