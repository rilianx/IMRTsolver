/*
 * ApertureILS.cpp
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
  
  void printHeader();
  
private:
  bool search_intensity;
  bool search_aperture;
  double prob_intensity;
};

}

#endif /* APERTUREILS_H_ */
