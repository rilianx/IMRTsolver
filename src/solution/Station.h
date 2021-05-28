
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

enum StationSetup {
  open_all_min = 1,
  open_all_max = 2,
  closed_all_min = 3,
  closed_all_max = 4,
  rand_all_rand = 5,
  rand_int = 6,
  manual_all_manual = 7,
  open_min_min = 8, 
  open_min_k = 9
};


/**
 * A Station consists of an angle and a set of apertures with intensities
 * The apertures+intensities are mapped to an intensity matrix (I)
 */

class Station {
private:
  int angle;

  // Dose deposition matrices for each volume
  map<int, const Matrix*> D;

  // Maximum number of apertures
  int max_apertures;
  int max_intensity;
  int min_intensity;
  int initial_intensity;
  int step_intensity;
  int n_volumes;
  /** Apertures (representation 1):
   * Each aperture is represented by a vector of pairs A[i] = (x_ini, x_fin)
   * and an intensity
   *
   */

   // Range open (x_ini, x_fin) of row "r" for aperture d: A[d][r](x_ini, x_fin)
   vector<vector<pair<int,int> > > A;


   void clearIntensityMatrix();

   //These auxiliar variables to track last changes
   pair<pair<int,int>, pair<int,int>> last_mem;
   list<pair<int,double>> last_diff;
   vector<double> last_intensity;

public:

  Collimator& collimator;

  //  Apertures (representation 2):
  // Intensity for each beam of the collimator
  Matrix I;

   // Constructs a new Station
  Station(Collimator& _collimator, vector<Volume>& volumes, int _angle,
          int max_apertures, int max_intensity=28, int initial_intensity=0,
          int step_intensity=2, StationSetup setup=StationSetup::open_all_min, istringstream* fluence_map=NULL);

  Station(const Station &s);

  Station& operator=(const Station & s);

  virtual ~Station(){ };

  void initializeStation(StationSetup type);

  void setApertureShape (int a, int row, int start, int end);

  void setApertureShape (int a, int row, pair<int, int> p);

  pair<int, int> getApertureShape(int a, int row);

  void generate_random_intensities();


  void change_intensity(int i, int j, double intensity, list< pair< int, double > >* diff=NULL );
  list< pair< int, double > > get_changes_intensity_move(double intensity, double delta) const;

 // intensity of an aperture i
  vector<double> intensity;

  // Function to be used to get the index in the location
  // in the matrix I of the rows of matrix D
  pair<int,int> getPos(int beam) const;

  int getBeamIndex(pair<int,int> coord) const;

  // Get intensity of beam
  int getIntensity(int beam) const{
    pair<int,int> p = getPos(beam);
    return I(p.first,p.second);
  };

  int getMaxIntensityRow(int i) const{
    int max=-1;
    for(int j=0;j<I.nb_cols();j++){
      if(I(i,j)>max)
        max=I(i,j);
    }
    return max;
  };

  // Returns the total intesity of the apertures
  int get_sum_alpha(string strategy) const{
    int intens=0;
    if(strategy=="dao_ls"){
      for(auto i:intensity) intens+=i;
      return (intens);
    }else if(strategy=="ibo_ls"){
	    for (int i=0; i<collimator.getXdim();i++)
		    for (int j=0; j<collimator.getYdim(); j++)
          if(I(i,j)>intens) intens=I(i,j);
      return (intens);
    }
    return(intens);
  };

  int get_nb_apertures(string strategy){
    if(strategy=="dao_ls")
       return (intensity.size());
    else if(strategy=="ibo_ls")
       return (int2nb.size());
    return(-1);
  };

  //revert the last change
  //should be followed by a call to the incremental_eval procedure
  void revert(list< pair< int, double > >& diff){
    for(auto let:diff){
      pair <int,int> a = beam2pos[let.first];
      change_intensity(a.first, a.second, I(a.first, a.second) - let.second);
    }
  };

  // Function to generate the intensity matrix from the
  // defined apertures and the intensity vector
  void generateIntensityMatrix();

  // Generate apertures from an ADEQUATE intensity matrix
  void generateApertures();

  void setUniformIntensity(double i);

  void printIntensity(bool vector_form=false);

  void writeIntensity(ifstream& myfile);

  void printApertures();
  
  string toStringApertures();
  
  string toStringIntensities();

  void printAperture(int aperture);

  const Matrix& getDepositionMatrix(int o) const;

  int getAngle(){
    return angle;
  };

  int getNbApertures() {
    return(max_apertures);
  };

  int getNbBeamlets() const{
    return collimator.getNangleBeamlets(angle);
  };

  list<int> open_beamlets(int aperture);

  list<int> closed_beamlets(int aperture);

  // Increase the intensity of a set of beams
  // (possible movement of a local search algorithm)
  // return a list with the changed beamlets an their changes to be used by the incremental evaluation
  list< pair< int, double > > increaseIntensity(int beam, double intensity, int ratio=0);
  list< pair< int, double > > increaseIntensity_repair(int beam, double intensity, int ratio=0);

  //return the next intensity level of the cell (if it is feasible)
  double next_intensity(int i, int j) const;

  //return the previous intensity level of the cell (if it is feasible)
  double prev_intensity(int i, int j) const;

  //increase (or decrease) the intesities equal to intensity to intensity+1
  //return a list with the changes to be used by the incremental evaluation
  list< pair< int, double > > change_intensity(double intensity, double delta);

  //Generate a new intensity by increasing (or decreasing) the most internal (resp. external) beamlets
  //of the just below (resp. just above) intensity.
  list< pair< int, double > > generate_intensity(double intensity, bool increase);

  //undo the changes in diff
  void diff_undo(list< pair< int, double > >& diff);

  void reduce_apertures(list< pair< int, double > >& diff);


  mutable map <pair<int,int>, int > pos2beam;
  mutable map <int, pair<int,int> > beam2pos;

  //maps from intensity to nb of open beamlets
  map< int, int > int2nb;

  //Aperture info and modification functions
  bool isOpenBeamlet (int beam, int aperture);
  bool isActiveBeamlet(int beam);
  vector<int> getClosed(int beam);
  vector<int> getOpen(int beam);
  bool anyClosed(int beam);
  bool anyOpen(int beam);

  /* Function that opens a beamlet from the left, if lside is true, or
     from the right size otherwise. Return true if the closing was performed.*/
  list<pair<int,double>> openBeamlet(int beam, int aperture);

  /* Function that closes a beamlet from the left, if lside is true, or
       from the right size otherwise. Return true if the closing was performed.*/
  list<pair<int,double>> closeBeamlet(int beam, int aperture, bool lside);
  list<pair<int,double>> closeBeamlet(pair<int,int> coord, int aperture, bool lside);
  list<pair<int,double>> modifyIntensityAperture(int aperture, double size);
  list<pair< int,double>> getModifyIntensityApertureDiff (int aperture, double size);
  void updateIntensity(list<pair<int,double>> diff);
  bool canIncreaseIntensity(int beam);
  bool canReduceIntensity(int beam);
  void setApertureIntensity(int aperture, double value);
  double getApertureIntensity(int aperture);
  int getMaxIntensity() const;

  list<pair <int,double>> closeRow(int row, int aperture, bool lside);
  list<pair <int,double>> getCloseRowDiff(int row, int aperture, bool lside);
  list <pair<int,double>> openRow(int row, int aperture, bool lside);
  list <pair<int,double>> getOpenRowDiff(int row, int aperture, bool lside);

  list <pair<int,double>> undoLast ();

  void clearHistory();

  /*static const int OPEN_MIN_SETUP = 0;
  static const int OPEN_MAX_SETUP = 1;
  static const int CLOSED_MIN_SETUP = 2;
  static const int CLOSED_MAX_SETUP = 3;
  static const int RAND_RAND_SETUP = 4;
  static const int RAND_INTENSITIES = 5; //only for ILS
  static const int MANUAL_SETUP=6;*/

};
}

#endif /* STATION_H_ */
