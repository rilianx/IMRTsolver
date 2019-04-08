//
//  Ant.h
//  
//
//  Created by Leslie on 07/04/2019.
//

#ifndef Ant_h
#define Ant_h

#include <stdio.h>
#include "Plan.h"

using namespace std;

namespace imrt {
  class Ant {
  private:
    Plan * p;
    vector < map <int, vector<int>> > ref_plan;
    
  public:
    Collimator & collimator;
    int n_apertures;
    int n_stations;
    
    Ant (vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
         vector<Volume>& volumes, int _max_apertures, int max_intensity, int initial_intensity,
         int step_intensity);
    
    void printIntensity (int i);
    void printReferencePlan();
  };
  
}

#endif /* Ant_h */
