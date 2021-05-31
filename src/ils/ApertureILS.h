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

  /* bsize: number of beamlets to be used in the target beamlet heuristic
     vsize: number of voxels to be considered when selecting the targeted beamlets
     search_intensity: perform local search over intensity
     search_aperture: perform local search over aperture
     prob_intensity: probability to perform local search over intensity
     step_intensity: step size for intensity
     initial_temperature: initial temperature for acceptance criterion
     alpha: alpha value for acceptance criterion
     do_perturbate: boolean variable that indicates if perturbation must be performed
     acceptance: type of acceptnace criterion to be used
     ls_type: type of local search to be initialized
  */
  ApertureILS(int bsize, int vsize, double _prob_intensity,
              int _step_intensity);

  ApertureILS(const ApertureILS & ils);

  virtual ~ApertureILS() {};

  //pair <bool, pair<Station*, int> > getLSBeamlet(Plan& P);

  bool isBeamletModifiable(int beamlet, Station* station, bool open_flag) ;

  //double openBeamlet(int beamlet, int aperture, Station& station, double c_eval, Plan& P);

  //double closeBeamlet(int beamlet, int side, int aperture, Station& station, double c_eval,  Plan& P);

  //double iLocalSearch(Plan& P, double max_time, bool verbose=true);
  //double aLocalSearch(Plan& P, double max_time, bool verbose=true);

  //double simpleLocalSearch(Plan& P, bool verbose=true);

  //double perturbation (Plan& P);

  //bool perturbate(int no_improvement, int iteration);


  vector <NeighborMove> getNeighborhood(Plan& current_plan,
                                        NeighborhoodType ls_neighborhood,
                                        LSTarget ls_target);
  vector < NeighborMove > getShuffledIntensityNeighbors(Plan &P);
  vector < NeighborMove > getShuffledApertureNeighbors(Plan &P);
    //vector < NeighborMove > getOrderedApertureNeighbors(Plan &P);
  vector < NeighborMove > getShuffledNeighbors(Plan &P);
  pair<vector < NeighborMove >, vector<NeighborMove>> getFriendsIntensityNeighbors(Plan &P, NeighborMove target);
  pair <vector < NeighborMove >, vector < NeighborMove >> getFriendsApertureNeighbors (Plan &P, NeighborMove target);
  vector < NeighborMove > getFriendsNeighbors(Plan &P, NeighborMove target);
  vector < NeighborMove > getOrderedApertureNeighbors(Plan &P);
  //double get_delta_eval(Plan &P, NeighborMove move, list<pair<int, double> >& diff);
  list<pair<int, double> > get_changes_in_fm(Plan &current_plan, NeighborMove move);
  double applyMove (Plan &P, NeighborMove move);
  bool checkMove (Plan &P, NeighborMove move) const;
  string planToString(Plan &P);

  int getStepIntensity ();

private:
  bool search_intensity;
  bool search_aperture;
  double prob_intensity;
  int step_intensity;

  int perturbation_size;
  bool do_perturbate;

};

}

#endif /* APERTUREILS_H_ */
