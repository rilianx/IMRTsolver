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

public:


  ILS(vector<Evaluator*>& evaluators, int sf_eval=0, int of_eval=0) : evaluators(evaluators), best_evals(evaluators.size())  {
     SF_evaluator=evaluators[sf_eval];
     OF_evaluator=evaluators[of_eval];
  };


  virtual ~ILS() { };


  double perturbation(Plan& current_plan, PerturbationType perturbation_type,
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
      if(verbose) cout << "ILS perturbation current:" << SF_evaluator->eval(current_plan)  << endl;
      for (int i=0; i < perturbation_size; i++) {
        neighborhood = getNeighborhood(current_plan, neighborhood_type, ls_target);
        move = neighborhood[i];

        //Update delta_eval
        list<pair<int, double> > changes = get_changes_in_fm(current_plan, move); //changes in fluence_map
        if(changes.empty()){i--; continue;}
        //double delta_eval = SF_evaluator->get_delta_eval (changes, current_plan.get_station(move.station_id)->getAngle());

        //Update delta_eval
        applyMove(current_plan, move);
        SF_evaluator->incremental_eval(changes, current_plan.get_station(move.station_id)->getAngle());

        if(verbose)
        cout << " -move type " << move.type << ", s:" <<
              move.station_id << ", a:" << move.aperture_id <<
              ", r:" << move.beamlet_id << ", action:"<< move.action << endl;
      }
    }
    return(SF_evaluator->get_evaluation());
  };

  virtual bool perturbate(int no_improvement, int iteration) {
    return(false);
  };


  virtual vector <NeighborMove> getNeighborhood(Plan & current_plan,
                                                NeighborhoodType ls_neighborhood,
                                                LSTarget ls_target= {LSTargetType::target_none, {0,0,0,0,0}}) = 0;

  //Update delta_eval
  virtual list<pair<int, double> > get_changes_in_fm(Plan &current_plan, NeighborMove move) = 0;
  //virtual double get_delta_eval(Plan &P, NeighborMove move, list<pair<int, double> >& diff) = 0;

  //debería llamar a incremental_eval
  virtual double applyMove (Plan & current_plan, NeighborMove move) = 0;

  virtual string planToString(Plan & current_plan) = 0;

  double iteratedLocalSearch (Plan& current_plan, int max_evaluations,
                              vector<NeighborhoodType>& ls_neighborhoods, int perturbation_size, ofstream& output_stream,
			                        int& used_evaluations, std::clock_t time_begin = clock(), bool verbose=false, 
                              double min_delta_eval=0.001, double alpha=1.0, int switch_patience=5, double pr_first_neigh=1.0);

  double FILocalSearch (Plan& current_plan, int max_evaluations,
                              vector<NeighborhoodType>& ls_neighborhoods, ofstream& output_stream,
			                        int& used_evaluations, std::clock_t time_begin, bool verbose, 
                              double& min_delta_eval, double alpha=1.0, int switch_patience=5, double pr_first_neigh=1.0);

  vector<Evaluator*> evaluators;
  vector<double> best_evals;

  //Search Function
  Evaluator* SF_evaluator;

  //Objective Function
  Evaluator* OF_evaluator;


};

}

#endif /* ILS_H_ */
