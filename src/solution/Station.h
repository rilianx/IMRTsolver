
/*
 * Station.h
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */

#include <map>
#include <vector>
#include <list>
#include <iterator>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "Collimator.h"
#include "tools/Matrix.h"
#include "Volume.h"

#ifndef STATION_H_
#define STATION_H_

using namespace std;
using namespace maths;

using namespace std;
using namespace maths;


namespace imrt {

/**
 * A Station consists of an angle and a set of apertures with intensities
 * The apertures+intensities are mapped to an intensity matrix (I)
 */

class Station {
private:

  Collimator& collimator;

  int angle;

  // Dose deposition matrices for each volume
  map<int, const Matrix*> D;

  // Maximum number of apertures
  int max_apertures;

  /** Apertures (representation 1):
   * Each aperture is represented by a vector of pairs A[i] = (x_ini, x_fin)
   * and an intensity
   *
   */

   // Range open (x_ini, x_fin) of row "r" for aperture d: A[d][r](x_ini, x_fin)
   vector<vector<pair<int,int> > > A;

   // intensity of an aperture i
   vector<double> intensity;

   //  Apertures (representation 2):
   // Intensity for each beam of the collimator
   Matrix I;


   void change_intensity(int i, int j, double intensity, list< pair< int, double > >* diff=NULL );

   void clearIntensity();



public:
  Station(Collimator& _collimator, vector<Volume>& volumes, int _angle, int _aperture);

  // Function to be used to get the index in the location
  // in the matrix I of the rows of matrix D
  pair<int,int> getPos(int beam) const;

  // Get intensity of beam
  int getIntensity(int beam) const{
    pair<int,int> p = getPos(beam);
    return I(p.first,p.second);
  }


  //revert the last change
  //should be followed by a call to the incremental_eval procedure
  void revert(list< pair< int, double > >& diff){
    for(auto let:diff){
      pair <int,int> a = beam2pos[let.first];
      change_intensity(a.first, a.second, I(a.first, a.second) - let.second);
    }
  }


  // Function to generate the intensity matrix from the
  // defined apertures and the intensity vector
  void generateIntensity();

  void setUniformIntensity(double i);

  void printIntensity(bool vector_form=false);

  void printApertures();



  const Matrix& getDepositionMatrix(int o) const;

  int getAngle(){ return angle;}

  int getNbBeamlets() const{
    return collimator.getNangleBeamlets(angle);
  }

  // Increase the intensity of a set of beams
  // (possible movement of a local search algorithm)
  // return a list with the changed beamlets an their changes to be used by the incremental evaluation
  list< pair< int, double > > increaseIntensity(int beam, double intensity, int ratio=0);
  list< pair< int, double > > increaseIntensity_repair(int beam, double intensity, int ratio=0);


  mutable map <pair<int,int>, int > pos2beam;
  mutable map <int, pair<int,int> > beam2pos;

   //maps from intensity to nb of open beamlets
   map< int, int > int2nb;
};
}

#endif /* STATION_H_ */
