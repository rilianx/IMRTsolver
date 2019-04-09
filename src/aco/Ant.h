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
#include <float.h>
#include "ApertureILS.h"

using namespace std;

namespace imrt {
  class Ant {
  private:
    Plan * p;
    vector < map <int, vector<int>> > ref_plan;
    vector< vector<int>> ref_intensity;
    
  public:
    Collimator & collimator;
    int n_apertures;
    int n_stations;
    int step_intensity;
    int max_intensity;
    
    Ant (vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
         vector<Volume>& volumes, int _max_apertures, int max_intensity, int initial_intensity,
         int step_intensity);
    Ant (const Ant& a);
    ~Ant ();
    
    void copy(Ant& a);
    
    void printIntensity (int i);
    void printIntensities ();
    void printReferencePlan();
    void generateReference();
    int getReferenceAt(int s, int ap, int x);
    int getIntensityReferenceAt(int s, int ap);
    
    void generateTour(vector<map <int, Matrix>> & probability, vector<Matrix> &iprobability);
    
    double getEvaluation();
    Plan * getPlan();
  };
  
}

#endif /* Ant_h */
