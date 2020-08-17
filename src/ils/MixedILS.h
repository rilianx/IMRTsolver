/*
 * MixedILS.h
 *
 *  Created on: 18 mar. 2020
 *      Author: ignacio
 */

#include "ApertureILS.h"
#include "IntensityILS2.h"

#ifndef MIXEDILS_H_
#define MIXEDILS_H_

namespace imrt {

class MixedILS {

IntensityILS2 ibo;
ApertureILS dao;

public:
    double total_evals;
    MixedILS(int bsize, int vsize, double prob_intensity,
                int step_intensity) : dao(bsize, vsize, prob_intensity, step_intensity) { }

    double iteratedLocalSearch (Plan& P, int max_time, int max_evaluations,
                                LSType ls_type_IBO, LSType ls_type_DAO,
				bool continuous,
                                NeighborhoodType ls_neighborhood_IBO,
                                NeighborhoodType ls_neighborhood_DAO,
                                LSTargetType ls_target_type_IBO,
                                LSTargetType ls_target_type_DAO,
                                PerturbationType perturbation_type,
                                int perturbation_size, int tabu_size, string convergence_file,
                                int evaluations=0, std::clock_t begin_time = clock()) {

        double F =P.getEvaluation();
        double bestF = F;
        total_evals = 0;

        Plan* best_plan=new Plan(P);

        //fast convergence
        //NeighborhoodType neighborhood = NeighborhoodType::intensity;
        //F = ibo.iteratedLocalSearch(P, max_time, 30, ls_type_IBO,
        //  neighborhood, ls_target_type_IBO, perturbation_type,
        //  0, 0 /*tabu size*/, convergence_file, total_evals, begin_time);
        //total_evals=ibo.used_evaluations;

        while(true){
          F = ibo.iteratedLocalSearch(P, max_time, max_evaluations, ls_type_IBO, continuous,
            ls_neighborhood_IBO, ls_target_type_IBO, perturbation_type,
            0, tabu_size, convergence_file, total_evals, begin_time);
          total_evals=ibo.used_evaluations;
          //plan for perturbations

          //cout << "P:" << P.eval() << endl;
          Plan P_pert(P);
          //cout << "Ppert:" << P_pert.eval() << endl;

          for(int i=0;i<5;i++)
            P.printIntensity(i);

          P.generateApertures();
          for(auto s:P.get_stations()) s->generateIntensityMatrix();

          F = dao.iteratedLocalSearch(P, max_time, max_evaluations, ls_type_DAO, continuous,
            ls_neighborhood_DAO, ls_target_type_DAO, perturbation_type,
      			  0, tabu_size, convergence_file, total_evals, begin_time);

          total_evals=dao.used_evaluations;

          //TODO: save best plan
          if(F<bestF){
            bestF=F;
            delete best_plan;
            best_plan = new Plan(P);
          }


          // Termination criterion
          std::clock_t time_end = clock();
          double used_time = double(time_end - begin_time) / CLOCKS_PER_SEC;
          if (max_time!=0 && used_time >= max_time) break;
          if (max_evaluations!=0 && total_evals>=max_evaluations) break;
          if (perturbation_size == 0)  break;

          //Perturbation
          cout << "pert" << endl;
          ibo.perturbation(P_pert, perturbation_type, perturbation_size);
          cout << "Ppert:" << P_pert.eval() << endl;
          P.newCopy(P_pert);
          cout << "P:" << P.eval() << endl;
        }
        P.newCopy(*best_plan);
        delete best_plan;
        return bestF;

    }

};

}

#endif /* MIXEDILS_H_ */
