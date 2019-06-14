/*
 * ApertureILS.h
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "ILS.h"

#ifndef APERTUREILS_H_
#define APERTUREILS_H_

namespace imrt {


class ApertureILS : public ILS {
public:
  
  /* bsize: number of beamlets to be used in the target beamlet heuristic
     vsize: number of voxels to be considered when selecting the targeted beamlets
     search_intensity: perform local search over intensity
     search_aperture: perform local search over aperture
     prob_intensity: probability to perform local search over intensity
     step_intensity: step size for intensity
     initial_temperature: initial temperature for acceptance criterion
     alpha: alpha value for acceptance criterion
     do_perturbate: boolean variable that indicates if perturbation must be performed
     acceptance: type of acceptnace criterion to be used
     ls_type: type of local search to be initialized
  */
  ApertureILS(int bsize, int vsize, bool search_intensity, 
              bool search_aperture, double prob_intensity, 
              int step_intensity, double initial_temperature, 
              double alpha, bool do_perturbate, 
              int perturbation_size, int acceptance, int ls_type);

  ApertureILS(const ApertureILS & ils);
  
  virtual ~ApertureILS() {};  

  pair <bool, pair<Station*, int> > getLSBeamlet(Plan& P);
  
  bool isBeamletModifiable(int beamlet, Station* station, bool open_flag) ;
  
  double improvementIntensity(int beamlet, Station& station, bool open_beamlet, 
                                   double c_eval, Plan &P, bool best_improvement); 
  
  double openBeamlet(int beamlet, int aperture, Station& station, double c_eval, Plan& P);
  
  double closeBeamlet(int beamlet, int side, int aperture, Station& station, double c_eval,  Plan& P);
  
  double improvementAperture(int beamlet, Station& station, bool open_beamlet, 
                             double c_eval, Plan& P, bool best_improvement); 
  
  bool acceptanceCriterion(double new_eval, double prev_eval);
  
  double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P);

  double localSearch(int type, int target, Plan& P);
  
  double iLocalSearch(Plan& P, double max_time, bool verbose=true);
  double aLocalSearch(Plan& P, double max_time, bool verbose=true);


  //double aiLocalSearch (Plan& P, double max_time, bool verbose=true);
  double simpleLocalSearch(Plan& P, bool verbose=true);
  
  void updateTemperature();
  
  double perturbation (Plan& P);
  
  bool perturbate(int no_improvement, int iteration);
  
 
  vector <NeighborMove> getNeighborhood(Plan& current_plan, 
                                        NeighborhoodType ls_neighborhood, 
                                        int ls_target);
  vector < NeighborMove > getShuffledIntensityNeighbors(Plan &P);
  vector < NeighborMove > getShuffledApertureNeighbors(Plan &P);
  vector < NeighborMove > getOrderedApertureNeighbors(Plan &P);
  vector < NeighborMove > getShuffledNeighbors(Plan &P);
  double applyMove (Plan &P, NeighborMove move);

  int getStepIntensity ();
  
  static const int FIRST_IMPROVEMENT=0;
  static const int BEST_IMPROVEMENT=1;
private:
  bool search_intensity;
  bool search_aperture;
  double prob_intensity;
  int step_intensity;
  int ls_type;
  
  double temperature;
  double initial_temperature;
  double max_temperature;
  double alpha;
  int perturbation_size;  
  bool do_perturbate;

  list<pair<Station*, int>> tabu;
};

}

#endif /* APERTUREILS_H_ */
