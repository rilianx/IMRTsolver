
/*
* Colimador.h
*
*  Created on: 22 may. 2018
*      Author: leslie
*/

#include <map>
#include <set>
#include <vector>
#include <list>
#include <iterator>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include "../tools/Matrix.h"

#ifndef COLIMADOR_H_
#define COLIMADOR_H_

using namespace std;
using namespace maths;

namespace imrt {

class Collimator {
private:

  // Overall coordinates
  //  beam_coord[x][y] coordinates as in all files, x rows and y cols
  //  note: in the format of the files
  map<double, vector<double> > beam_coord;

  // Per angle coordinates
  map<int, vector<pair<double,double> > >  angle_coord;

  // Ordered coordinates
  vector<double> xcoord;
  vector<double> ycoord;

  int nb_beamlets;
  int xdim, ydim;

  // Insert a new X coordinate in the xcoord vector
  void insertXorder(double x);
  // Insert a new Y coordinate in the ycoord vector
  void insertYorder(double y);

  // Number of angles available
  //int nb_angles;

  // Range (i,j) of active beamlets of angle "a" row "r":
  //  angle_row_beam_active[a][r](i,j)
  //  (-1,-1) indicates a full row is not active
  map<int, vector<pair<int,int> > > angle_row_beam_active;

  // Active beamlets per angle
  void setActiveBeamlets(map<int,  vector <pair<double,double> > >& coord);

  // Number of beamlets defined for an angle "a": nb_angle_beamlets[a]
  map<int, int> nb_angle_beamlets;

  list<int> angles;

  int n_angles;

  vector< pair<int, string> > coord_files;

  map < int, pair <int,int> > reference;
  vector<int> levels_to_intensity;
  map<double, int> intensity_to_levels;

  int ref_size;

public:
  Collimator(){};
  Collimator(string coord_filename, set<int> angles);
  Collimator(vector< pair<int, string> >& coord_files);
  Collimator(const Collimator& c);
  Collimator& operator=(Collimator& c);

  void initializeCoordinates(vector < pair<int, string> >& coord_files);
  // Overall coordinates pair (x,y)
  void printCoordinates();
  // Values of the overall coordinates axis
  void printAxisValues();

  // Size of the overall matrix
  int getXdim();
  int getYdim();

  //int getNangles();

  // Number of beamlets active in an angle
  int getNangleBeamlets (int angle) const;

  // Indicates if a beam located in position (x,y) is active
  bool isActiveBeamAngle(int x, int y, int angle);

  pair<int, int> getActiveRange(int x, int angle);

  pair<int,int> indexToPos(int index, int angle);

  void printActiveBeam();

  const list<int>& getAngles() const;

  int getAngle(int i) const;

  int getNbAngles() const;

  static string delimiter;

  void generateReference ();

  void printReference ();

  int getReferenceSize();

  pair <int, int> getReference (int r) ;

  int searchReferenceIndex (pair<int, int>);

  void generateIntensityLevels(int step_intensity, int max_intensity);

  void printIntensityLevels ();

  int getIntensityLevelSize();

  int getIntensityLevel(double v);

  double getLevelIntensity(int l);
};
}

#endif /* COLIMADOR_H_ */
