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
  
  ApertureILS(int bsize, int vsize, bool search_intensity, bool search_aperture, double prob_intensity, double initial_temperature, double alpha, int acceptance);
  
  pair <bool, pair<Station*, int> > getLSBeamlet(Plan& P);
  
  bool isBeamletModifiable(int beamlet, Station* station, bool open_flag) ;
  
  double firstImprovementIntensity(int beamlet, Station& station, bool open_beamlet, 
                                   double c_eval, Plan &P); 
  
  double doOpen(int beamlet, int aperture, Station& station, double c_eval, Plan& P);
  
  double doClose(int beamlet, int aperture, Station& station, double c_eval,  Plan& P);
  
  double firstImprovementAperture(int beamlet, Station& station, bool open_beamlet, 
                             double c_eval, Plan& P); 
  
  bool acceptanceCriterion(double new_eval, double prev_eval);
  
  double localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P);
  
  void updateTemperature();
  
private:
  bool search_intensity;
  bool search_aperture;
  double prob_intensity;
  
  double temperature;
  double initial_temperature;
  double max_temperature;
  double alpha;
  
};

}

#endif /* APERTUREILS_H_ */
