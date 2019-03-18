/*
 * ApertureILS.h
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "ACO.h"

#ifndef ACS_H_
#define ACS_H_

namespace imrt {

class ACS : public ACO {
public:
  
  ACS(vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
      vector<Volume>& volumes, int max_apertures, int max_intensity, int initial_intensity,
      int step_intensity, int _n_ants, double _initial_pheromone, double _alpha, double _beta);

  virtual ~ACS() {};
  
  void generateTours();
  void calculateProbability();
  void printProbability();

  //double ailocalsearch (Plan &P);
  
  static const int FIRST_IMPROVEMENT=0;
  static const int BEST_IMPROVEMENT=1;
private:
 
  double q0;
};

}

#endif /* APERTUREILS_H_ */
