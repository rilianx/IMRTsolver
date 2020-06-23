/*
 * IntensityILS2.h
 *
 *  Created on: 19 jun. 2019
 *      Author: ignacio
 */


#include "ILS.h"

#ifndef INTENSITYILS2_H_
#define INTENSITYILS2_H_

namespace imrt {


class IntensityILS2 : public ILS {
public:

  static int vsize;
  static double min_improvement;

  IntensityILS2() : ILS() { };

  IntensityILS2(const IntensityILS2 & ils) : ILS(ils) {
  };

  virtual ~IntensityILS2() {};

  vector <NeighborMove> getNeighborhood(Plan& current_plan,
                                        NeighborhoodType ls_neighborhood,
                                        LSTarget ls_target);
  vector < NeighborMove > getShuffledIntensityNeighbors(Plan &P);
  vector < NeighborMove > getShuffledApertureNeighbors(Plan &P);
  vector < NeighborMove > getShuffledApertureNeighbors_target(Plan &P, int nv=5);
  vector < NeighborMove > getOrderedApertureNeighbors(Plan &P);
  vector < NeighborMove > getShuffledNeighbors(Plan &P);


  //return interleaved sets of neighbors: i1, a1, i2, a2,... ik, ak
  vector < NeighborMove > getShuffledNeighbors(Plan &P, int k, bool intensity);

  //recorrer el vector y deja los "friends" al comienzo
  void prioritizeFriends(vector < NeighborMove >& neighborList, LSTarget ls_target, Plan& current_plan);

  double applyMove (Plan &P, NeighborMove move, bool p);

  double applyMove (Plan &P, NeighborMove move){
	  return applyMove(P,move,false);
  }

  double applyMoveP (Plan &P, NeighborMove move){
	  return applyMove(P,move,true);
  }

  string planToString(Plan &P);
private:
  static int myrandom (int i) { return std::rand()%i;}

  vector < NeighborMove > ShuffledApertureNeighbors;
  int last_move;


};

}

#endif /* APERTUREILS_H_ */
