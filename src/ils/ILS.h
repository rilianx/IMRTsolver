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
#include "EvaluatorF.h"
#include "EvaluatorGScore.h"

namespace imrt {

struct NeighborMove {
  // Type
  // 1: intensity
  // 2: aperture
  int type;
  int station_id;
  int aperture_id;
  // Action
  // 1: open left (increase)
  // 2: open right (increase)
  // -1 close left (reduce intensity)
  // -2 close right
  // -3 close from closest border
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
  smixed_a= 9
};

enum LSType {
  best = 1,
  first = 2
};

enum LSTargetType {
  target_none = 1,
  target_friends = 2
};

struct LSTarget {
  LSTargetType target_type;
  NeighborMove target_move;
};

// Basically this ones will be the same as the moves
enum PerturbationType {
  p_intensity = 1,
  p_aperture = 2,
  p_mixed = 3,
  none = 4
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

  /*virtual double iLocalSearch(Plan& P,double max_time, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };

  virtual double aLocalSearch(Plan& P,  double max_time, bool verbose=true) {
    cout << "Not implemented "<< endl;
    return 0.0;
  };*/

  double perturbation(Plan& current_plan, Evaluator& evaluator, PerturbationType perturbation_type,
                              int perturbation_size, bool verbose=true) {
    vector <NeighborMove> neighborhood;
    bool is_move = false;
    NeighborhoodType neighborhood_type;
    NeighborMove move;
    LSTarget ls_target = {LSTargetType::target_none, {0,0,0,0,0}};

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
      if(verbose) cout << "ILS perturbation current:" << evaluator.eval(current_plan)  << endl;
      for (int i=0; i < perturbation_size; i++) {
          neighborhood = getNeighborhood(current_plan, neighborhood_type,
                                         ls_target);
    	  move = neighborhood[i];

        //Update delta_eval
        list<pair<int, double> > changes = get_changes_in_fm(current_plan, move); //changes in fluence_map
        double delta_eval = evaluator.get_delta_eval (changes, current_plan.get_station(move.station_id)->getAngle());
       

        //list<pair<int, double> > diff;
        //double delta_eval = get_delta_eval(current_plan, move, diff);
        //Station *s = current_plan.get_station(move.station_id);

        if(changes.empty()){i--; continue;}

        //Update delta_eval
        applyMove(current_plan, move);
        evaluator.incremental_eval(changes, current_plan.get_station(move.station_id)->getAngle());

      if(verbose)
    	  cout << " -move type " << move.type << ", s:" <<
               move.station_id << ", a:" << move.aperture_id <<
               ", r:" << move.beamlet_id << ", action:"<< move.action << endl;
      }
    }
    return(evaluator.get_evaluation());
  };

  virtual bool perturbate(int no_improvement, int iteration) {
    return(false);
  };

  //virtual void undoLast(Plan& p){
	//  p.undoLast();
  //}

  virtual vector <NeighborMove> getNeighborhood(Plan & current_plan,
                                                NeighborhoodType ls_neighborhood,
                                                LSTarget ls_target) = 0;

  //Update delta_eval
  virtual list<pair<int, double> > get_changes_in_fm(Plan &current_plan, NeighborMove move) = 0;
  //virtual double get_delta_eval(Plan &P, NeighborMove move, list<pair<int, double> >& diff) = 0;

  //debería llamar a incremental_eval
  virtual double applyMove (Plan & current_plan, NeighborMove move) = 0;

  virtual string planToString(Plan & current_plan) = 0;

  int used_evaluations;

  double iteratedLocalSearch (Plan& current_plan, Evaluator& evaluator, int max_time, int max_evaluations,
                              LSType ls_type, bool continuous, NeighborhoodType ls_neighborhood,
                              LSTargetType ls_target_type, PerturbationType perturbation_type,
                              int perturbation_size, int tabu_size, string convergence_file,
			      int evaluations=0, std::clock_t begin = clock(), bool verbose=false) {
     int current_iteration = 0;
     double aux_eval = evaluator.get_evaluation();
     double best_eval = aux_eval;
     double used_time = 0;
     used_evaluations = evaluations;
     Plan* best_plan = new Plan(current_plan);

     //Start time
     std::clock_t time_end;
     time_begin = begin;

     //Create log files
     ofstream c_file;
     ofstream t_file;
     string trajectory_file="";
     if (convergence_file!="") {
       c_file.open (convergence_file.c_str(), ios::out);
       trajectory_file = convergence_file + ".traj";
       t_file.open (trajectory_file.c_str(), ios::app);
       if (evaluations == 0) {
         t_file << "evalutions;time;quality" << endl;
         t_file << "0;0"  << ";" << best_eval << "\n";
       }
       t_file.close();
     }

     cout << "Starting iterated local search " << endl;
     while (true) {
       // Apply local search
       if (ls_type == LSType::first) {
         aux_eval = FILocalSearch(current_plan, evaluator, max_time, max_evaluations,
                    used_evaluations, ls_neighborhood, ls_target_type,
		    tabu_size, trajectory_file, continuous, verbose);
       } else {
         aux_eval = BILocalSearch(current_plan, evaluator, max_time, max_evaluations,
                    used_evaluations, ls_neighborhood, ls_target_type,
                    tabu_size, trajectory_file, verbose);
       }

       if (aux_eval < best_eval) {
         
         best_eval = aux_eval;
         delete best_plan;
         best_plan = new Plan(current_plan);
                  
         //cout << best_plan->eval()  << endl;
       }

       current_iteration++;
       if(verbose)
          cout << "ILS  iteration: " << current_iteration << " ; current: " <<
            aux_eval << " ; best: " << best_eval << "; evals:" << used_evaluations << endl;

       // Print convergence information
       if (convergence_file!="") {

         c_file << used_evaluations << ";" << current_iteration <<
                   ";" << aux_eval << ";" << best_eval <<
                    planToString(current_plan) << endl;

       }

       // Termination criterion
       time_end = clock();
       used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
       if (max_time!=0 && used_time >= max_time) break;
       if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
       if (perturbation_size == 0)  break;

       cout << 0 << endl;
       current_plan.newCopy(*best_plan);
       //Perturbation
       perturbation(current_plan, evaluator, perturbation_type, perturbation_size, verbose);

     }

     if (c_file.is_open()) {
       c_file.close();
     }
     cout << 1 << endl;
     current_plan.newCopy(*best_plan);
     cout << 2 << endl;
     delete best_plan;
     return(best_eval);
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
         user == NeighborhoodType::sequential_a ) {
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


  void addUndoTabu (vector<NeighborMove>& tabu_list, NeighborMove move, int tabu_size) {
    NeighborMove tabu_move = move;

    if (tabu_list.size() == tabu_size)
      tabu_list.erase(tabu_list.begin());

    if (move.type == 1) {
      //Intensity move
      if (move.action < 0){
        //Move reduces intensity
        tabu_move.action = 1;
      } else {
        //Move increases intensity
        tabu_move.action = -1;
      }
    } else if (move.type == 2) {
      //Aperture move
      if (move.action < 0){
        //Move closes aperture
        tabu_move.action = 1;
      } else {
        //Move open aperture
        tabu_move.action = -3;
      }
    }
    tabu_list.push_back(tabu_move);
  };


  void addTabu (vector<NeighborMove>& tabu_list, NeighborMove move, int tabu_size) {
    NeighborMove tabu_move = move;

    if (tabu_list.size() == tabu_size)
      tabu_list.erase(tabu_list.begin());

    tabu_list.push_back(tabu_move);
  };

  bool isTabu (NeighborMove move, vector<NeighborMove> tabu_list) {
    for (auto tabu_move : tabu_list) {
      if (tabu_move.beamlet_id != move.beamlet_id) continue;
      if (tabu_move.aperture_id != move.aperture_id) continue;
      if (tabu_move.station_id != move.station_id) continue;
      if (tabu_move.type != move.type) continue;
      if (tabu_move.action < 0 && move.action > 0) continue;
      if (tabu_move.action > 0 && move.action < 0) continue;
      //cout << "CENSORED! "<< move.station_id <<"," << move.aperture_id <<"," << move.beamlet_id << "," << move.action << endl;

      return(true);
    }
    return(false);
  };

  double FILocalSearch (Plan& current_plan, Evaluator& evaluator, int max_time, int max_evaluations,
                        int& used_evaluations, NeighborhoodType ls_neighborhood,
                        LSTargetType ls_target_type, int tabu_size,
                        string trajectory_file, bool continuous, bool verbose=false) {

    EvaluatorGS evalGS(evaluator);
    evalGS.eval(current_plan);

    EvaluatorF evalF(evaluator);
    evalF.eval(current_plan);

    vector <NeighborMove> neighborhood;
    double current_eval = evaluator.get_evaluation();
    NeighborMove best_move = {0,0,0,0,0};
    NeighborhoodType current_neighborhood;
    LSTarget ls_target = {ls_target_type, best_move};

    vector <NeighborMove> tabu_list;
    NeighborMove move;

    int id_neighbor = 0;
    int n_neighbors = 1;
    bool generate_neighborhood = true;

    bool improvement = true;

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

    cout << "Local search" << endl;

    while (improvement) {
      if (generate_neighborhood){
        //Select the neighborhood (done for the cases in which sequenced neighborhoods are chosen)
        current_neighborhood = selectNeighborhood (current_neighborhood,
                                                   ls_neighborhood,
                                                   improvement && !sequential_flag); // improvement is always true?

        //Get the moves in the neighborhood (this is possible because the moves are not that many!)
        neighborhood = getNeighborhood(current_plan, current_neighborhood, ls_target);
        id_neighbor = 0;
	      n_neighbors = 0;
	      generate_neighborhood = false;
      }

      improvement = false;

      if(verbose)
        cout << " -neighborhood: " << current_neighborhood << "; size: " <<
             neighborhood.size() << "; curren best: " << current_eval << endl;;

      while (n_neighbors < neighborhood.size()) {

        move = neighborhood[id_neighbor];

        //Counter updates
        n_neighbors++; id_neighbor++;
        if (id_neighbor >= neighborhood.size()) id_neighbor = 0;

        //Skip neighbor if its marked as tabu
        if (tabu_size > 0 && isTabu(move, tabu_list))  continue;
        
        list<pair<int, double> > changes = get_changes_in_fm(current_plan,move); //changes in fluence_map
        //double delta_eval = current_plan.get_delta_eval (move.station_id, changes);
        double delta_eval = evaluator.get_delta_eval (changes, current_plan.get_station(move.station_id)->getAngle());

        //Update delta_eval (eliminar 3 lineas)
        //list<pair<int, double> > diff;
        //double delta_eval = get_delta_eval(current_plan, move, diff);
        //Station *s = current_plan.get_station(move.station_id);

        if(changes.empty()) continue;

        //Evaluation updates
        used_evaluations++;

        // Check if there is an improvement
        if (delta_eval < -0.001) {
            applyMove(current_plan, move); //NO llama a incremental_eval
            current_eval = evaluator.incremental_eval (changes,  current_plan.get_station(move.station_id)->getAngle());
            
        //these evals compute the evaluation using the already updated z structure
        cout << used_evaluations <<",";
        cout << evalF.incremental_eval() << ",";
        cout << evalGS.incremental_eval() <<endl;

        if(false)
            cout << "    neighbor: " << n_neighbors  << "; "
              << "(" << move.station_id <<   "," << move.beamlet_id << "," << move.action
              << "); improvement: " << evaluator.get_evaluation() << ";" << endl;
              //<< evalGS.incremental_eval(changes,current_plan.get_station(move.station_id)->getAngle()) <<endl;

            improvement = true;

            if (ls_neighborhood==NeighborhoodType::smixed_i && move.type==2) {
              ls_neighborhood=NeighborhoodType::smixed_a;
              generate_neighborhood = true;
            } else {
              if (ls_neighborhood==NeighborhoodType::smixed_a && move.type==1) {
                ls_neighborhood=NeighborhoodType::smixed_i;
                generate_neighborhood = true;
              }
          }

          if (!continuous)
              generate_neighborhood = true;
          
          ls_target.target_type = ls_target_type;
          ls_target.target_move = move;

          // Add the undo movements as tabu
          if (tabu_size > 0)
            addUndoTabu(tabu_list, move, tabu_size);

            n_neighbors = 0;
            break;

        } else {
          //no improvement

          //Update delta_eval (eliminar)
          //s->undoLast();
          //s->diff_undo(diff);


	        // Add the non-improving movement as tabu
          if (tabu_size > 0)
            addTabu(tabu_list, move, tabu_size);
        }

        // Termination criterion
        time_end = clock();
        used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
        if (max_time!=0 && used_time >= max_time) break;
        if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
      }

      if (improvement) {
	      time_end = clock();
        used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
        if (t_file.is_open())
           t_file << used_evaluations << ";" << used_time << ";" << current_eval << "\n";

      }

      //Check if imixed neighborhood
      if (!improvement && ls_neighborhood==NeighborhoodType::imixed &&
	      current_neighborhood==NeighborhoodType::intensity){
    	    ls_neighborhood=NeighborhoodType::mixed;
    	    improvement = true;
	        generate_neighborhood = true;
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
	  generate_neighborhood = true;
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

    if(verbose)
      cout << endl;


    return(current_eval);
  };

  double BILocalSearch (Plan& current_plan, Evaluator& evaluator, int max_time, int max_evaluations,
                      int& used_evaluations,  NeighborhoodType ls_neighborhood,
                      LSTargetType ls_target_type, int tabu_size,
                      string trajectory_file, bool verbose=true) {

    vector <NeighborMove> neighborhood;
    double current_eval = evaluator.get_evaluation();
    NeighborMove best_move = {0,0,0,0,0};
    NeighborhoodType current_neighborhood;
    LSTarget ls_target = {LSTargetType::target_none, best_move};
    vector <NeighborMove> tabu_list;

    int n_neighbors = 0;
    bool improvement = true;

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


    cout << "Local search" << endl;

    while (improvement) {

      improvement = false;

      //Update counter
      n_neighbors = 0;

      //Select the neighborhood (done for the cases in which sequenced neighborhoods are chosen)
      current_neighborhood = selectNeighborhood (current_neighborhood,
                                                   ls_neighborhood,
                                                   improvement && !sequential_flag); // improvement is always true?
      //Get the moves in the neighborhood (this is possible because the moves are not that many!)
      neighborhood = getNeighborhood(current_plan, current_neighborhood, ls_target);

      if(verbose)
        cout << " -neighborhood: " << current_neighborhood << "; size: " <<
        neighborhood.size() << " ";

       for (NeighborMove move:neighborhood) {

        //Counter updates
        n_neighbors++;

        //cout << "  Neighbor: " << n_neighbors  << "(" << move.station_id <<
        //          "," << move.aperture_id << "," <<move.beamlet_id << ","<< move.action << ");" << endl;

        //Skip neighbor if its marked as tabu
        if (tabu_size > 0 && isTabu(move, tabu_list)) continue;

        //-1.0 means that the move is not a valid move
        if(applyMove(current_plan, move) == -1.0)  continue;

        //Evaluation updates
        used_evaluations++;

        // Check if there is an improvement
        if ( evaluator.get_evaluation() < (current_eval-0.001)) {
          if(verbose)
            cout << "; neighbor: " << n_neighbors  << "(" << move.station_id <<
              "," << move.aperture_id << "," << move.action << "); Improvement: " <<
               evaluator.get_evaluation() << endl;

          improvement = true;

          if(ls_neighborhood==NeighborhoodType::smixed_i && move.type==2)
            ls_neighborhood=NeighborhoodType::smixed_a;
          else if(ls_neighborhood==NeighborhoodType::smixed_a && move.type==1)
            ls_neighborhood=NeighborhoodType::smixed_i;

          // Best improvement
          current_eval =  evaluator.get_evaluation();
          best_move = move;
          //current_plan.undoLast(); //TODO: ignacio revisame!.
          ls_target.target_type = ls_target_type;
          ls_target.target_move = move;

        } else {
          //No improvement
          //current_plan.undoLast();
          // Add the non-improving movement as tabu
          if (tabu_size > 0)
            addTabu(tabu_list, move, tabu_size);
        }

        // Termination criterion
        time_end = clock();
        used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
        if (max_time!=0 && used_time >= max_time) break;
        if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
      }

      // Apply best improvement move
      if (improvement) {
        applyMove(current_plan, best_move);
        current_eval =  evaluator.get_evaluation();
        current_plan.clearLast();
        // Add the undo movements as tabu
        if (tabu_size > 0)
          addUndoTabu(tabu_list, best_move, tabu_size);
        //No se agregan evaluaciones en este caso por que
        //no es una solucion nueva...
        //used_evaluations++;
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

};

}

#endif /* ILS_H_ */
