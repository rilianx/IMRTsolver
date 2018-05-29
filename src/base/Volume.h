/*
 * Volume.h
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

#include "../tools/Matrix.h"
#include "Collimator.h"

#ifndef VOLUME_H_
#define VOLUME_H_

using namespace std;
using namespace maths;


namespace imrt {

  /**
   * Information of a volume (organ):
   * - Dose deposition matrices for each angle
   * - Number of voxels
   */
  class Volume {
    private:
      /* Information of the general beamlet configuration */
	  Collimator & collimator;

      /* Dose deposition matrices for each angle
       * D[a](k,b): Dose delivered to voxel k of the organ by the beamlet b and angle a */
      map<int, Matrix> D;

      /* Mapping of beamlets into rows
       * R[a][r]<<x,y>,<i1,i2>>: angle a, row r from beamlet i1 to i2 */
      //map<double,vector<pair<double,double> > > R;
      //map<double, vector<int> > r_active_size;

      /* number of voxels of the volume */
      int nb_voxels;


      void set_data(string file);

    public:
      Volume(Collimator& collimator, string deposition_file);

      void print_deposition();
      void print_coordinates();
      int getNbVoxels() {return nb_voxels;}
      //void print_r_active_size();

      Matrix& getDepositionMatrix(int angle);

  };
}

#endif /* VOLUME_H_ */
