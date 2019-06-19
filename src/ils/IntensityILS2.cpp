/*
 * IntensityILS2.cpp
 *
 *  Created on: 19 jun. 2019
 *      Author: ignacio
 */

#include "IntensityILS2.h"

namespace imrt {

/*
struct NeighborMove {
  // Type
  // 1: beamlet
  // 2: intensity
  int type;
  int station_id;
  int aperture_id; // intensity that will be changed
  // Action
  //+1: increase
  //-1: decrease
  int action;
  int beamlet_id;
};*/


vector < NeighborMove > IntensityILS2::getShuffledIntensityNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  vector< NeighborMove > moves;

  int k=0;
  for(auto s:stations){

	  for(pair <int, int> intensity : s->int2nb){
		  if (intensity.first < s->getMaxIntensity()){
			  NeighborMove move = {2,k,intensity.first,+1,0};
			  moves.push_back(move);
		  }

		  NeighborMove move = {2,k,intensity.first,-1,0};
		  moves.push_back(move);
	  }
	  k++;
  }

  std::random_shuffle ( moves.begin(), moves.end(), myrandom);
  return(moves);
};

vector < NeighborMove > IntensityILS2::getShuffledBeamletNeighbors(Plan &P){

   list<Station*> stations = P.get_stations();
   vector < NeighborMove > moves;

   int k=0;
   for(auto s:stations){
 	   for (int i=0; i<s->I.nb_rows();i++)
			 for (int j=0; j<s->I.nb_cols(); j++)
				 if(s->I(i,j)!=-1){
					 double intensity=s->intensityUp(i,j);
					 if(intensity != s->I(i,j)){
						 NeighborMove move = {1,k,0,+1,s->pos2beam[make_pair(i,j)]};
						 moves.push_back(move);
					 }

					 intensity=s->intensityDown(i,j);
					 if(intensity != s->I(i,j)){
						 NeighborMove move = {1,k,0,-1,s->pos2beam[make_pair(i,j)]};
						 moves.push_back(move);
					 }
				 }
 	  k++;
   }

   std::random_shuffle ( moves.begin(), moves.end(), myrandom);
   return moves;
 }



vector < NeighborMove > IntensityILS2::getShuffledNeighbors(Plan &P) {

  vector< NeighborMove > a_list, final_list;

  final_list = getShuffledIntensityNeighbors(P);
  a_list = getShuffledBeamletNeighbors(P);
  final_list.insert(final_list.end(), a_list.begin(), a_list.end());

  std::random_shuffle(final_list.begin(), final_list.end());

  return(final_list);
};

vector < NeighborMove> IntensityILS2::getNeighborhood(Plan& current_plan,
                                       NeighborhoodType ls_neighborhood,
                                       LSTarget ls_target){
  vector < NeighborMove> neighborList;

  if (ls_neighborhood == intensity) {
    neighborList = getShuffledIntensityNeighbors(current_plan);
  } else if (ls_neighborhood == aperture) {
    neighborList = getShuffledBeamletNeighbors(current_plan);
  } else {
    //mixed
    neighborList = getShuffledNeighbors(current_plan);
  }
  return(neighborList);
}



double IntensityILS2::applyMove (Plan & current_plan, NeighborMove move) {

  double current_eval = current_plan.getEvaluation();

  list<pair<int, double> > diff;
  int type            = move.type;
  Station * s         = current_plan.get_station(move.station_id);
  int beamlet         = move.beamlet_id;
  int action          = move.action;
  double intens	  = move.aperture_id;

  int i = s->beam2pos[beamlet].first;
  int j = s->beam2pos[beamlet].second;

  if (type == 1) {
    // single beam
    double intensity;

    if (action < 0) // decrease intensity
      intensity=s->intensityDown(i,j);
     else // Increase intensity
      intensity=s->intensityUp(i,j);


  	 double delta_eval = current_plan.get_delta_eval (*s, i, j, intensity-s->I(i,j));
     if(delta_eval < 0.001){
	   s->change_intensity(i, j, intensity);
	   current_eval = current_plan.incremental_eval(*s, i, j, intensity-s->I(i,j));
     }

  } else {
    // change intensity
    if (action < 0)
      diff = s->change_intensity(intens, -1.0);
    else
      diff = s->change_intensity(intens, +1.0);

    double delta_eval = current_plan.get_delta_eval (*s, diff);
    if(delta_eval < 0.001)
    	current_eval = current_plan.incremental_eval(*s, diff);
    else
    	s->diff_undo(diff);

  }

  return(current_eval);

};



}



