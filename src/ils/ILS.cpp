/*
 * ApertureILS.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "ILS.h"

namespace imrt {

double ILS::iteratedLocalSearch (Plan& current_plan, int max_evaluations,
                              vector<NeighborhoodType>& ls_neighborhoods, int perturbation_size, ofstream& output_stream,
			                  int& used_evaluations, std::clock_t time_begin, bool verbose, 
                              double min_delta_eval, double alpha, int switch_patience, double pr_first_neigh) {

    cout << "Starting iterated local search " << endl;

    double best_eval = SF_evaluator->get_evaluation();
    Plan* best_plan = new Plan(current_plan);

    while (true) {
        // Apply local search
        double eval = FILocalSearch(current_plan, max_evaluations,
                ls_neighborhoods, output_stream, used_evaluations, time_begin, verbose, min_delta_eval, alpha, switch_patience, pr_first_neigh);

        cout << SF_evaluator->eval(current_plan) << endl;
        if (eval < best_eval) {
            best_eval = eval;
            delete best_plan;
            cout << "updated best plan" << endl;
            best_plan = new Plan(current_plan);
            cout << SF_evaluator->eval(*best_plan) << endl;
        }

        current_plan=*best_plan; 
        pr_first_neigh=1.0;

        // Termination criterion
        if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
        if (perturbation_size == 0)  break;

        
        //Perturbation
        cout << "perturbation" << endl;
        perturbation(current_plan, PerturbationType::p_mixed, perturbation_size, verbose);

    }

    return best_eval;
};


double ILS::FILocalSearch (Plan& current_plan, int max_evaluations,
                              vector<NeighborhoodType>& ls_neighborhoods, ofstream& output_stream,
			                  int& used_evaluations, std::clock_t time_begin, bool verbose, 
                              double& min_delta_eval, double alpha, int switch_patience,  double pr_first_neigh) {
    cout << "Starting first improvement local search " << endl;
    vector <NeighborMove> neighborhood;
    double current_eval = SF_evaluator->eval(current_plan);
    double best_OF = OF_evaluator->incremental_eval(); int OF_no_improv=0;

    int count=0;

    // Select initial neighborhood
      
    int id_neigh=0;
    int last_improved_neigh=0;
    while (true) {
        
        NeighborhoodType current_neighborhood = ls_neighborhoods[id_neigh];
       
        //Get the moves in the neighborhood (this is possible because the moves are not that many!)
        neighborhood = getNeighborhood(current_plan, current_neighborhood);
    
        //the first move improving eval more than min_delta_eval is selected
        bool improvement = false;
        list<pair<int, double>> changes;
        NeighborMove move;
        for(int n = 0; n< (double) neighborhood.size()*pr_first_neigh; n++){
            move = neighborhood[n];

            //get the changes required by the move
            changes = get_changes_in_fm(current_plan,move); //changes in fluence_map
            if(changes.empty()) continue;

            //get delta_eval
            double delta_eval = SF_evaluator->get_delta_eval (changes, current_plan.get_station(move.station_id)->getAngle());

            //Evaluation updates
            used_evaluations++;
            min_delta_eval*=alpha; 

            if( -delta_eval >= min_delta_eval ) {
                improvement = true;
                break;
            }

            // Termination criterion
            if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
            
        }

        // Check if there is an improvement
        if (improvement) {
            applyMove(current_plan, move); //NO llama a incremental_eval
            current_eval = SF_evaluator->incremental_eval (changes,  current_plan.get_station(move.station_id)->getAngle());
            last_improved_neigh = id_neigh;

            //update evaluation of other evaluators
            for(int i=0; i<evaluators.size(); i++)
                if (evaluators[i]!=SF_evaluator) evaluators[i]->incremental_eval();
            

            if(SF_evaluator!=OF_evaluator){
                double OF = OF_evaluator->get_evaluation();
                if(best_OF == 0.0 || OF < best_OF) {
                    OF_no_improv=0;
                    best_OF=OF;
                    //guardar solucion best_OF
                    if(best_planOF) *best_planOF=current_plan;
                    else best_planOF=new Plan(current_plan);
                }else OF_no_improv++;  

                if (switch_patience>0 && OF_no_improv>=switch_patience){
                    cout << "switch evaluator" << endl;
                    SF_evaluator=OF_evaluator;
                     //comenzar con solucion best_OF
                     if(best_planOF){
                         current_plan=*best_planOF;
                         SF_evaluator->eval(current_plan);
                         OF_evaluator->incremental_eval();
                     }
                    current_eval = SF_evaluator->get_evaluation();
                }
            }
            

            cout << "Scores: " << used_evaluations << " " ;
            for (auto score : dynamic_cast<EvaluatorGS*>(OF_evaluator)->scores )
                cout << score.value << " ";
            for(int i=0; i<evaluators.size(); i++){
                double ev=evaluators[i]->get_evaluation();
                cout << ev << ";";
                if(best_evals[i]==0 || ev<best_evals[i])
                    best_evals[i]=ev;
            }
            cout << current_neighborhood << ";" << endl;
        }

     
        //save results
        if (improvement && output_stream.is_open()) {
            clock_t time_end = clock();
            double used_time = double(time_end - time_begin) / CLOCKS_PER_SEC;
            output_stream <<  used_evaluations <<";";
            output_stream <<  used_time <<";";
            for (auto score : dynamic_cast<EvaluatorGS*>(OF_evaluator)->scores )
                output_stream << score.value << ";";
            for(auto ev:evaluators)
                output_stream << ev->get_evaluation() << ";";
            output_stream << current_neighborhood << ";";

            output_stream << min_delta_eval <<";";
            output_stream << (SF_evaluator==OF_evaluator) << endl;
        }


        
        // Termination criterion
        if (max_evaluations!=0 && used_evaluations>=max_evaluations) break;
        
        //change neighbourhood
        if(!improvement){
            pr_first_neigh=1.0;
            id_neigh = (id_neigh+1) % ls_neighborhoods.size(); 
        }

        // cycled over all the neighbourhoods without improvement
        if (!improvement && id_neigh == last_improved_neigh) break; 

        
    }

    return(current_eval);
  };



}