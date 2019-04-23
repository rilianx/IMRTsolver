/*
 *  ILS.h
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#ifndef ILS_H_
#define ILS_H_

#include <ctime>
#include <algorithm>
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

  ILS(int bsize, int vsize, int acceptance=ACCEPT_NONE): bsize(bsize), vsize(vsize), acceptance(acceptance){
  };

  ILS(const ILS & ils ){
    bsize=ils.bsize;
    vsize=ils.vsize;
    acceptance=ils.acceptance;
  };

  virtual ~ILS() { };

  virtual double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P) = 0;
  
  virtual double iLocalSearch(Plan& P,double max_time, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };
  
  virtual double aLocalSearch(Plan& P,  double max_time, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };
  
  virtual double iSLocalSearch(Plan& P, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };
  
  virtual double aSLocalSearch(Plan& P, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };
  
  virtual bool acceptanceCriterion(double new_eval, double prev_eval)=0;

  virtual pair<bool, pair<Station*, int>> getLSBeamlet(Plan& P){
	  return P.getLSBeamlet(bsize, vsize);
  }

  virtual pair<bool, pair<Station*, int>> getBestLSBeamlet(Plan& P){
    return P.getBestLSBeamlet(bsize, vsize);
  }

  virtual double perturbation(Plan& P) {
    return(P.getEvaluation());
  };

  virtual bool perturbate(int no_improvement, int iteration) {
    return(false);
  };

  virtual void undoLast(Plan& p){
	  p.undoLast();
  }

  virtual void updateTemperature() {};

  //  Plan* init_plan;
  double beamTargetedSearch(Plan& current_plan, int max_time, int max_iterations) {

    cout << "## Staring ILS search." << endl;
    std::clock_t time_end;

    time_begin=clock();

  //  if(init_plan) delete init_plan;
  //  init_plan = new Plan(current_plan);


    Plan& best_plan= *new Plan(current_plan);

    pair<bool, pair<Station*, int>> target_beam;
    double local_eval, aux_eval,  best_eval=current_plan.eval();
    double used_time=0;
    bool flag=true;
    int no_improvement, iteration=1, perturbation_iteration=0;
    no_improvement = 0;

    local_eval=best_eval;
    //Start time


    while (flag) {
      //cout << "ss"<< endl;
      target_beam = getLSBeamlet(current_plan);
      //cout << target_beam.second.second << endl;
      while (target_beam.second.second < 0) {
        cout << "NOTE: No beamlet available." << endl;
        local_eval = perturbation(current_plan);
        perturbation_iteration=iteration;
        target_beam = getLSBeamlet(current_plan);
        if (local_eval < best_eval) {
          best_eval=local_eval;
          best_plan.newCopy(current_plan);
        }
      }

      cout << "Iteration: " << iteration << ", eval: " << EvaluationFunction::n_evaluations << ", time: "<< (roundf(used_time * 1000) / 1000)  << ", best: " << best_eval <<
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
        undoLast(current_plan);
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

      if ( perturbate(no_improvement, iteration )) {
        local_eval = perturbation(current_plan);
        perturbation_iteration = iteration;
        no_improvement=no_improvement/2;
      }
    }
    current_plan.newCopy(best_plan);
    aux_eval=current_plan.getEvaluation();
    cout << 1 << endl;
    best_plan.getEvaluationFunction()->generate_voxel_dose_functions();
    cout << 2 << endl;
    return(aux_eval);
  };

  /* Not targeted version having first ils and then als */
  double notTargetedSearch(Plan& current_plan, int max_time, int max_iterations) {

    cout << "## Staring ILS search." << endl;
    std::clock_t time_end;

    //Start time
    time_begin=clock();
    Plan& best_plan= *new Plan(current_plan);

    double local_eval, aux_eval,  best_eval;
    double used_time=0;
    double ls_time=0;
    bool flag=true, impils=true, impals=true;
    int no_improvement, iteration=1, perturbation_iteration=0;
    no_improvement = 0;
    
    best_eval = local_eval = current_plan.getEvaluation();
    while (flag) {
      cout << "Iteration: " << iteration << ", eval: " << EvaluationFunction::n_evaluations << 
        ", time: "<< (roundf(used_time * 1000) / 1000)  << ", best: " << best_eval << ", current: " << local_eval << endl;
      
      // Perturbate if we haven't improved last iteration
      if (!impals) {
        local_eval = perturbation(current_plan);
        cout << "Iteration: " << iteration << ", per: " << local_eval << endl;
      }
      impals = impils = false;

      if (max_time!=0) ls_time= max_time-used_time;

      
      // Intensity ls
      aux_eval =iLocalSearch (current_plan, ls_time, true);
      if (aux_eval < local_eval) {
        local_eval = aux_eval;
        impils = true;
      }
      
      // Termination criterion
      time_end = clock();
      used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) {
        flag = false;
        if (local_eval < best_eval) {
          best_eval = local_eval;
          best_plan.newCopy(current_plan);
        }
        break;
        cout << "## Timed out" << endl;
      }

      if (max_time!=0) ls_time= max_time-used_time;
      
      // Aperture ls
      aux_eval = aLocalSearch (current_plan, ls_time , false);
      if (aux_eval < local_eval) {
        local_eval = aux_eval;
        impals = true;
      }
      
      /*
      cout << endl;
      for(int j=0;j<5;j++)
        current_plan.printIntensity(j);
      cout << endl;
      */
      
      if (local_eval < best_eval) {
        best_eval = local_eval;
        best_plan.newCopy(current_plan);
      }
      iteration++;

      // Termination criterion
      time_end=clock();
      used_time=double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) flag=false;
      if (max_iterations!=0 && iteration>=max_iterations) flag=false;
    }
    current_plan.newCopy(best_plan);
    aux_eval=current_plan.getEvaluation();
    best_plan.getEvaluationFunction()->generate_voxel_dose_functions();
    
    time_end=clock();
    used_time=double(time_end- time_begin) / CLOCKS_PER_SEC;
    cout << "## Total used time: "<< used_time << endl;
    return(aux_eval);
  };

};

}

#endif /* ILS_H_ */
