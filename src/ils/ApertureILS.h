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
  
  ApertureILS(int bsize, int vsize, bool search_intensity, bool search_aperture, double prob_intensity, 
              int step_intensity, double initial_temperature, double alpha, bool do_perturbate, 
              int perturbation_size, int acceptance, int ls_type);

  virtual ~ApertureILS() {};  

  pair <bool, pair<Station*, int> > getLSBeamlet(Plan& P);
  
  bool isBeamletModifiable(int beamlet, Station* station, bool open_flag) ;
  
  double improvementIntensity(int beamlet, Station& station, bool open_beamlet, 
                                   double c_eval, Plan &P, bool best_improvement); 
  
  double doOpen(int beamlet, int aperture, Station& station, double c_eval, Plan& P);
  
  double doClose(int beamlet, int aperture, Station& station, double c_eval,  Plan& P);
  
  double improvementAperture(int beamlet, Station& station, bool open_beamlet, 
                             double c_eval, Plan& P, bool best_improvement); 
  
  bool acceptanceCriterion(double new_eval, double prev_eval);
  
  double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P);
  
  void updateTemperature();
  
  double perturbation(Plan& P);
  
  bool perturbate(int no_improvement, int iteration);
  
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
