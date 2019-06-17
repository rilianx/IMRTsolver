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

struct NeighborMove {
  // Type
  // 1: intensity
  // 2: aperture
  int type;
  int station_id;
  int aperture_id;
  // Action
  // 1: open (increase)
  // -1 close left (reduce intensity)
  // -2 close right
  int action;
  int beamlet_id;
};

enum NeighborhoodType {
  intensity = 1,
  aperture = 2,
  mixed = 3,
  sequential = 4
};

enum LSType {
  best = 1,
  first = 2
};

enum LSTarget {
  none = 1,
  beamlet = 2  
};

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

  ILS(int bsize, int vsize, int acceptance=ACCEPT_NONE): bsize(bsize),
                                  vsize(vsize), acceptance(acceptance){
  };

  ILS(const ILS & ils ){
    bsize=ils.bsize;
    vsize=ils.vsize;
    acceptance=ils.acceptance;
  };

  virtual ~ILS() { };
  
  virtual double iLocalSearch(Plan& P,double max_time, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };
  
  virtual double aLocalSearch(Plan& P,  double max_time, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };
  
  virtual bool acceptanceCriterion(double new_eval, double prev_eval) = 0;

  virtual pair<bool, pair<Station*, int>> getLSBeamlet(Plan& P){
	  return P.getLSBeamlet(bsize, vsize);
  };

  virtual pair<bool, pair<Station*, int>> getBestLSBeamlet(Plan& P){
    return P.getBestLSBeamlet(bsize, vsize);
  };

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

  virtual vector <NeighborMove> getNeighborhood(Plan & current_plan,
                                                NeighborhoodType ls_neighborhood,
                                                LSTarget ls_target) = 0;

  virtual double applyMove (Plan & current_plan, NeighborMove move) = 0;

  double iteratedLocalSearch (Plan& current_plan, int max_time, int max_evaluations, 
                              LSType ls_type, NeighborhoodType ls_neighborhood,
                              LSTarget ls_target) {
     int current_iteration = 0;
     double aux_eval = current_plan.getEvaluation();
     double current_eval = current_plan.getEvaluation();
     double used_time = 0;
     int used_evaluations = 0;

     //Start time
     std::clock_t time_end;
     time_begin = clock();

     cout << "Starting iterated local search " << endl;
     while (true) {
       // Apply local search
       aux_eval = localSearch(current_plan, max_time, max_evaluations,  
                  used_evaluations, ls_type, ls_neighborhood, ls_target);
        
       if (aux_eval < current_eval) {
         current_eval = aux_eval;
       }

       current_iteration++;
       cout << "  iteration: " << current_iteration << " ; eval: " << aux_eval << " ; best: " << current_eval << endl;

       // Termination criterion
       time_end = clock();
       used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
       if (max_time!=0 && used_time >= max_time) break;
       if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;

       //Perturbation
       perturbation(current_plan);

     }
     return(current_eval);
  }

  double localSearch (Plan& current_plan, int max_time, int max_evaluations, 
                      int& used_evaluations, LSType ls_type, NeighborhoodType ls_neighborhood, 
                      LSTarget ls_target) {
    bool improvement = true;
    vector <NeighborMove> neighborhood;
    double current_eval = current_plan.getEvaluation();    
    NeighborMove best_move = {0,0,0,0,0}; 
    //Start time
    std::clock_t time_end;
    double used_time;

    cout << "Starting local search" << endl;

    while (improvement) {
      improvement = false;
      //TODO: Hacer enum para ls_neighborhood, ls_target, ls_type
      neighborhood = getNeighborhood(current_plan, ls_neighborhood, ls_target); //TODO: implementar
      cout << " Neighborhood size: " << neighborhood.size() << endl;
      for (NeighborMove move:neighborhood) {
        applyMove(current_plan, move);  //TODO: implementar
        // Check if there is an improvement
        if (current_plan.getEvaluation() < current_eval) {
            cout << "  Improvement" << current_plan.getEvaluation()<< endl;
          improvement = true;
          if (ls_type == LSType::first) {
            // First improvement  
            current_eval = current_plan.getEvaluation();
            break;
          } else {
            // Best improvement
            current_eval = current_plan.getEvaluation();
            best_move = move;
            current_plan.undoLast(); //TODO: ignacio implementame.
          }
        } else {
          current_plan.undoLast();
        }
       
        used_evaluations++;

        // Termination criterion
        time_end = clock();
        used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
        if (max_time!=0 && used_time >= max_time) break;
        if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
      }

      // Apply best improvement move    
      if (ls_type == LSType::best && improvement) {
        applyMove(current_plan, best_move);
        current_eval = current_plan.getEvaluation();
        used_evaluations++;
      }     

      // Termination criterion
      time_end = clock();
      used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) break;
      if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
    }
    return(current_eval);
  };


  /* Targeted version of the local search where a beamlert is identified as 
     promising and local search is directed to it */
  /*double beamTargetedSearch (Plan& current_plan, int max_time, int max_iterations) {

    cout << "## Staring ILS search." << endl;

    //Start time
    std::clock_t time_end;
    time_begin=clock();

    //Best plan found so far
    Plan& best_plan= *new Plan(current_plan);

    pair<bool, pair<Station*, int>> target_beam;
    double local_eval, aux_eval,  best_eval=current_plan.eval();
    double used_time=0;
    // flag is used for the termination criterion
    bool flag = true;
    int no_improvement, iteration=1, perturbation_iteration=0;
    no_improvement = 0;

    local_eval=best_eval;
    
    while (flag) {
      // Get the targeted beamlet for local search
      target_beam = getLSBeamlet(current_plan);
      while (target_beam.second.second < 0) {
        // If there is no beamlet available perturbate
        // until we get one
        cout << "NOTE: No beamlet available." << endl;
        local_eval = perturbation(current_plan);
        perturbation_iteration = iteration;
        target_beam = getLSBeamlet(current_plan);
        if (local_eval < best_eval) {
          best_eval=local_eval;
          best_plan.newCopy(current_plan);
        }
      }

      cout << "Iteration: " << iteration << 
              ", eval: " << 
              EvaluationFunction::n_evaluations << 
              ", time: " << 
              (roundf(used_time * 1000) / 1000)  << 
              ", best: " << best_eval <<
              ", current: " << local_eval  << 
              ", beamlet: " << target_beam.second.second  <<
              ", station: " << target_beam.second.first->getAngle() << 
              ", +-: " << target_beam.first;

      // Apply local search on the targeted beamlet
      aux_eval = localSearch (target_beam, current_plan);
      cout << endl;

      // Check if there is a new global best solution 
      if (aux_eval < best_eval) {
        best_eval=aux_eval;
        best_plan.newCopy(current_plan);
      }

      // Acceptance criterion
      if (aux_eval < local_eval) {
        // Accept improving solution
        local_eval = aux_eval;
        no_improvement = 0;
      } else if (aux_eval!= local_eval && acceptanceCriterion(aux_eval, local_eval)) {
        // Accept non-improving solution
        local_eval = aux_eval;
        no_improvement = 0;
      } else {
        // Not accept solution
        undoLast(current_plan);
        no_improvement ++;
      }

      // Update SA temperature
      if (acceptance==ACCEPT_SA)
         updateTemperature();

      iteration++;

      // Termination criterion
      time_end=clock();
      used_time=double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) flag = false;
      if (max_iterations!=0 && iteration>=max_iterations) flag = false;

      // Check if we should perturbate
      if ( perturbate(no_improvement, iteration )) {
        local_eval = perturbation(current_plan);
        perturbation_iteration = iteration;
        no_improvement=no_improvement/2;
      }
    }

    current_plan.newCopy(best_plan);
    aux_eval=current_plan.getEvaluation();
    best_plan.getEvaluationFunction()->generate_voxel_dose_functions();
    return(aux_eval);
  };*/

  /* Not targeted version having first ils and then als */
 /* double notTargetedSearch(Plan& current_plan, int max_time, int max_iterations) {

    cout << "## Staring ILS search." << endl;
    std::clock_t time_end;

    //Start time
    time_begin=clock();

    //Best plan found so far
    Plan& best_plan= *new Plan(current_plan);

    double local_eval, aux_eval,  best_eval;
    double used_time=0;
    double ls_time=0;
    // flag is used for the termination criterion
    bool flag=true;
    // improvement signals a new best solution found
    bool improvement=true;
    int no_improvement, iteration=1, perturbation_iteration=0;
    no_improvement = 0;
    
    best_eval = local_eval = current_plan.getEvaluation();
    while (flag) {
      cout << "Iteration: " << iteration << ", eval: " << 
              EvaluationFunction::n_evaluations << 
              ", time: "<< (roundf(used_time * 1000) / 1000)  << 
              ", best: " << best_eval << ", current: " << 
              local_eval << endl;
      
      // Perturbate if we haven't improved last iteration
      if (!improvement) {
        local_eval = perturbation(current_plan);
        cout << "Iteration: " << iteration << ", per: " << local_eval << endl;
      }
      
      improvement = false;

      // Track the used time
      if (max_time!=0) ls_time = max_time - used_time;
      
      // Apply intensity ls
      aux_eval = iLocalSearch (current_plan, ls_time, false);
      // Check if there if there is a new solution was found
      if ((local_eval - aux_eval) > 0.00001) {
        local_eval = aux_eval;
      }
      
      // Check if there is a better global best
      if ((best_eval - local_eval) > 0.00001) {
        best_eval = local_eval;
        best_plan.newCopy(current_plan);
      }
      
      //Check print (comment if not needed)
      //current_plan.eval();
      //cout << "returned eval: " << aux_eval << " current local: " << local_eval << " in current plan: " << current_plan.getEvaluation()<< endl;
      //cout << endl;
      //for(int j=0;j<5;j++)
      //current_plan.printIntensity(j);
      //cout << endl;
      
      // Check termination criterion
      time_end = clock();
      used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) {
        flag = false;
        break;
      }
      // Update used time
      if (max_time!=0) ls_time = max_time-used_time;
      
      // Apply aperture ls
      aux_eval = aLocalSearch (current_plan, ls_time , false);
      // Check if there if there is a new solution was found
      if ((local_eval - aux_eval) > 0.00001) {
        local_eval = aux_eval;
        improvement = true;
      }
      
      // Check if there is a better global best
      if ((best_eval - local_eval) > 0.00001) {
        best_eval = local_eval;
        best_plan.newCopy(current_plan);
      }

      iteration++;

      // Check termination criterion
      time_end=clock();
      used_time=double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) flag = false;
      if (max_iterations!=0 && iteration>=max_iterations) flag = false;
    }

    current_plan.newCopy(best_plan);
    aux_eval = current_plan.getEvaluation();
    best_plan.getEvaluationFunction()->generate_voxel_dose_functions();
    
    time_end = clock();
    used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
    cout << "## Total used time: " << used_time << endl;
    return(aux_eval);
  };*/

};

}

#endif /* ILS_H_ */
