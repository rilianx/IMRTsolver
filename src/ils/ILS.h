/*
 *  ILS.h
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#ifndef ILS_H_
#define ILS_H_

#include <ctime>
#include <iostream>
#include <fstream>
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
  sequential_i = 4,
  sequential_a = 5,
  sequential_p = 6,
  imixed = 7,
  smixed_i= 8,
  smixed_a= 9,
};

enum LSType {
  best = 1,
  first = 2
};

enum LSTarget {
  none = 1,
  beamlet = 2  
};

// Basically this ones will be the same as the moves
enum PerturbationType {
  p_intensity = 1,
  p_aperture = 2,
  p_mixed = 3
};

class ILS {
private:
  double max_no_improvement=100;


public:
  int bsize;
  int vsize;
  
  std::clock_t time_begin;

  static const int ACCEPT_NONE = 0;
  static const int ACCEPT_SA = 1;

  ILS(int bsize=1, int vsize=10, int acceptance=ACCEPT_NONE): bsize(bsize),
                                  vsize(vsize){
  };

  ILS(const ILS & ils ){
    bsize=ils.bsize;
    vsize=ils.vsize;
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
  
  virtual pair<bool, pair<Station*, int>> getLSBeamlet(Plan& P){
	  return P.getLSBeamlet(bsize, vsize);
  };

  virtual pair<bool, pair<Station*, int>> getBestLSBeamlet(Plan& P){
    return P.getBestLSBeamlet(bsize, vsize);
  };

  double perturbation(Plan& current_plan, PerturbationType perturbation_type,
                              int perturbation_size) {
    vector <NeighborMove> neighborhood;
    bool is_move = false;
    NeighborhoodType neighborhood_type;
    NeighborMove move;

    // Move-based perturbation
    if (perturbation_type == PerturbationType::p_intensity){
      neighborhood_type = NeighborhoodType::intensity;
      is_move = true;
    } else if (perturbation_type == PerturbationType::p_aperture){
      neighborhood_type = NeighborhoodType::aperture;
      is_move = true;
    } else if (perturbation_type == PerturbationType::p_mixed){
      neighborhood_type = NeighborhoodType::mixed;
      is_move = true;
    }
 
    if (is_move ){
      cout << "Perturbation:" << endl; 
      for (int i=0; i < perturbation_size; i++) {
          neighborhood = getNeighborhood(current_plan, neighborhood_type,
                                         LSTarget::none);
    	  move = neighborhood[i];
    	  cout << "  move type " << move.type << ", s:" <<
               move.station_id << ", a:" << move.aperture_id << 
               ", b:" << move.beamlet_id << ", action:"<< move.action << endl;

		  applyMoveP(current_plan, move);
      }
    }
    return(current_plan.getEvaluation());
  };

  virtual bool perturbate(int no_improvement, int iteration) {
    return(false);
  };

  virtual void undoLast(Plan& p){
	  p.undoLast();
  }

  virtual vector <NeighborMove> getNeighborhood(Plan & current_plan,
                                                NeighborhoodType ls_neighborhood,
                                                LSTarget ls_target) = 0;

  virtual double applyMoveP (Plan & current_plan, NeighborMove move){
	  return applyMove (current_plan, move);
  }

  virtual double applyMove (Plan & current_plan, NeighborMove move) = 0;

  virtual string planToString(Plan & current_plan) = 0;

  int used_evaluations;

  double iteratedLocalSearch (Plan& current_plan, int max_time, int max_evaluations, 
                              LSType ls_type, NeighborhoodType ls_neighborhood,
                              LSTarget ls_target, PerturbationType perturbation_type,
                              int perturbation_size, string convergence_file) {
     int current_iteration = 0;
     double aux_eval = current_plan.getEvaluation();
     double current_eval = current_plan.getEvaluation();
     double used_time = 0;
     used_evaluations = 0;

     //Start time
     std::clock_t time_end;
     time_begin = clock();

     //Create log files
     ofstream c_file;
     ofstream t_file;
     string trajectory_file="";
     if (convergence_file!="") {
       c_file.open (convergence_file.c_str(), ios::out);
       trajectory_file = convergence_file + ".traj";
       t_file.open (trajectory_file.c_str(), ios::out);
       t_file << "evalutions;quality" << endl;
       t_file.close();
     }

     cout << "Starting iterated local search " << endl;
     while (true) {
       // Apply local search
       aux_eval = localSearch(current_plan, max_time, max_evaluations,  
                  used_evaluations, ls_type, ls_neighborhood, ls_target, 
                  trajectory_file);
        
       if (aux_eval < current_eval) {
         current_eval = aux_eval;
       }

       current_iteration++;
       cout << "  iteration: " << current_iteration << " ; eval: " << 
               aux_eval << " ; best: " << current_eval << endl;

       // Print convergence information
       if (convergence_file!="") {
         c_file << used_evaluations << ";" << current_iteration << 
                   ";" << aux_eval << ";" << current_eval << 
                    planToString(current_plan) << endl;
       }

       // Termination criterion
       time_end = clock();
       used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
       if (max_time!=0 && used_time >= max_time) break;
       if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;

       //Perturbation
       perturbation(current_plan, perturbation_type, perturbation_size);

     }

     if (c_file.is_open()) {
       c_file.close();
     }
     return(current_eval);
  };

  //Improvement indicates we should keep the neighborhood!
  NeighborhoodType selectNeighborhood (NeighborhoodType current, 
                                       NeighborhoodType user, 
                                       bool keep_flag) {
     
     if (user == NeighborhoodType::intensity || 
         user == NeighborhoodType::aperture ||
         user == NeighborhoodType::mixed ||
		 user == NeighborhoodType::smixed_i ||
		 user == NeighborhoodType::smixed_a)
       return (user);

     //TODO: add probability parameter
     if (user == NeighborhoodType::sequential_p) {
       if (keep_flag) {
         // If we found an improvement we choose randomly
         float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
         if (r < 0.5) return(NeighborhoodType::intensity);
         else return(NeighborhoodType::aperture);
       } else {
         if (current == NeighborhoodType::aperture)
           return (NeighborhoodType::intensity);
         return(NeighborhoodType::aperture); 
       }
     }

     if (user == NeighborhoodType::sequential_i ||
         user == NeighborhoodType::sequential_a) {
       if (keep_flag) {
         return (current);
       } else {
         if (current == NeighborhoodType::aperture)
           return (NeighborhoodType::intensity);
         return(NeighborhoodType::aperture);  
       }
     }

     if (user == NeighborhoodType::imixed)
    	 return (current);


     cout << "Error: Unknown neighborhood operator" << endl;
     return(NeighborhoodType::intensity);  
  };

  NeighborhoodType selectInitialNeighborhood (NeighborhoodType user) {
     if (user == NeighborhoodType::sequential_i || user == NeighborhoodType::imixed)
       return(NeighborhoodType::intensity);
     else if (user == NeighborhoodType::sequential_a)
       return(NeighborhoodType::aperture);
     else if (user == NeighborhoodType::sequential_p){
       float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
       if (r < 0.5) return(NeighborhoodType::intensity);
       else return(NeighborhoodType::aperture);
     } 
     return(user);
  };

  double localSearch (Plan& current_plan, int max_time, int max_evaluations, 
                      int& used_evaluations, LSType ls_type, 
                      NeighborhoodType ls_neighborhood, 
                      LSTarget ls_target, string trajectory_file) {
    bool improvement = true;
    vector <NeighborMove> neighborhood;
    double current_eval = current_plan.getEvaluation();    
    NeighborMove best_move = {0,0,0,0,0}; 
    NeighborhoodType current_neighborhood;
    int n_neighbor = 1; 

    // the sequential flag indicates that the previous neighborhood was checked
    // unsuccesfully 
    bool sequential_flag = false; 
    bool is_sequential = (ls_neighborhood == NeighborhoodType::sequential_i) ||
                         (ls_neighborhood == NeighborhoodType::sequential_a) ||
                         (ls_neighborhood == NeighborhoodType::sequential_p);

    // Select initial neighborhood
    current_neighborhood = selectInitialNeighborhood(ls_neighborhood);

    //Output file
    ofstream t_file;
    if (trajectory_file!="") {
      t_file.open(trajectory_file.c_str(), ios::app);
    }

    //Start time
    std::clock_t time_end;
    double used_time;   
    cout << "Starting local search" << endl;

    while (improvement) {
      //Select the neighborhood (done for the cases in which sequenced neighborhoods are chosen)
      current_neighborhood = selectNeighborhood (current_neighborhood, 
                                                 ls_neighborhood, 
                                                 improvement && !sequential_flag); // improvement is always true?
      improvement = false;
      //Get the moves in the neighborhood (this is possible because the moves are not that many!)
      neighborhood = getNeighborhood(current_plan, current_neighborhood, ls_target); 
      n_neighbor = 1;//neighbor counter

      cout << " Neighborhood: " << current_neighborhood << "; size: " << 
              neighborhood.size() << endl;

      for (NeighborMove move:neighborhood) {
        //Generate the solution in the neighborhood
        applyMove(current_plan, move);  //TODO: ignacio implementame!.


        // Check if there is an improvement
        if (current_plan.getEvaluation() < (current_eval-0.001)) {
          cout << "  Neighbor: " << n_neighbor  << "(" << move.station_id << 
                  "," << move.aperture_id << "," << move.action << "); Improvement: " <<
                  current_plan.getEvaluation() << endl;
          improvement = true;

          if(ls_neighborhood==NeighborhoodType::smixed_i && move.type==2)
        	  ls_neighborhood==NeighborhoodType::smixed_a;
          else if(ls_neighborhood==NeighborhoodType::smixed_a && move.type==1)
        	  ls_neighborhood==NeighborhoodType::smixed_i;


          if (ls_type == LSType::first) {
            // First improvement  
            current_eval = current_plan.getEvaluation();
            current_plan.clearLast();
            break;
          } else {
            // Best improvement
            current_eval = current_plan.getEvaluation();
            best_move = move;
            current_plan.undoLast(); //TODO: ignacio revisame!.
          }
        } else {
          //No improvement
          current_plan.undoLast();
        }
       
        //Counter updates 
        used_evaluations++;
        n_neighbor++;

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
        current_plan.clearLast();
        used_evaluations++;
      }

      if (improvement) {
         if (t_file.is_open())
            t_file << used_evaluations << ";" << used_time << ";" << current_eval << "\n";
         
      }

      //Check if imixed neighborhood
      if (!improvement && ls_neighborhood==NeighborhoodType::imixed && current_neighborhood==NeighborhoodType::intensity){
    	  ls_neighborhood=NeighborhoodType::mixed;
    	  improvement = true;
      }

      // Check sequential neighborhood
      if (is_sequential) {
        if (!improvement & !sequential_flag) {
          // If we are in sequential mode, when there is no 
          // improvement we allow to check also the next neighborhood
          // Note: this should be coordinated with the select neighborhood 
          // function.
		  sequential_flag = true;
          improvement = true;
        } else {
          sequential_flag = false;
        }
      }

      // Termination criterion
      time_end = clock();
      used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) break;
      if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
    }

    if (trajectory_file!="")
      t_file.close();

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
