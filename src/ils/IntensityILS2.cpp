/*
 * IntensityILS2.cpp
 *
 *  Created on: 19 jun. 2019
 *      Author: ignacio
 */

#include "IntensityILS2.h"
#include "Evaluator.h"



namespace imrt {


  double IntensityILS2::vsize=0.01;
  double IntensityILS2::min_improvement=0.05;

vector < NeighborMove > IntensityILS2::getShuffledIntensityNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  vector< NeighborMove > moves;

  int k=0;
  for(auto s:stations){

	  NeighborMove move = {2,k,0,+1,0};
		moves.push_back(move);

	  for(int i=0; i<s->getNbApertures(); i++){
			NeighborMove move = {2,k,i+1,+1,0};
		  moves.push_back(move);

		  move = {2,k,i+1,-1,0};
		  moves.push_back(move);
	  }

	  k++;
  }

  std::random_shuffle ( moves.begin(), moves.end(), myrandom);
  return(moves);
};

vector < NeighborMove > IntensityILS2::getShuffledApertureNeighbors_target(Plan &P, double vsize){
  vector < NeighborMove > best_moves;
  vector < NeighborMove > shuffled_moves;
  multimap < double, pair<int, int>, MagnitudeCompare2> beamlets =
    P.get_evaluator().sorted_beamlets(P, vsize);
  bool first=true;
  for(auto b : beamlets){
    Station *s = P.get_station(b.second.first);
    pair<int,int> ij = s->beam2pos[b.second.second];
    int i=ij.first, j=ij.second;
    double intensity=s->I(i,j);
    if(b.first > 0) {
      intensity=s->next_intensity(i,j);
        NeighborMove move = {1,b.second.first,0,+1,b.second.second};
        if(b.first>=abs(min_improvement)) best_moves.push_back(move);
        else shuffled_moves.push_back(move);

        move = {1,b.second.first,0,-1,b.second.second};
        shuffled_moves.push_back(move);
    }else if(b.first <= 0) {
      intensity=s->prev_intensity(i,j);
        NeighborMove move = {1,b.second.first,0,-1,b.second.second};
        if(b.first<=-abs(min_improvement)) best_moves.push_back(move);
        else shuffled_moves.push_back(move);

        move = {1,b.second.first,0,+1,b.second.second};
        shuffled_moves.push_back(move);
    }
  }

  if(min_improvement < 0.0)
    std::random_shuffle ( best_moves.begin(), best_moves.end(), myrandom);
  std::random_shuffle ( shuffled_moves.begin(), shuffled_moves.end(), myrandom);
  best_moves.insert(best_moves.end(), shuffled_moves.begin(), shuffled_moves.end());

  return best_moves;
}

vector < NeighborMove > IntensityILS2::getShuffledApertureNeighbors(Plan &P){

   list<Station*> stations = P.get_stations();
   vector < NeighborMove > moves;

   int k=0;
   for(auto s:stations){
 	   for (int i=0; i<s->I.nb_rows();i++)
			 for (int j=0; j<s->I.nb_cols(); j++)
				 if(s->I(i,j)!=-1){
					 double intensity=s->next_intensity(i,j);
					 if(intensity != s->I(i,j)){
						 NeighborMove move = {1,k,0,+1,s->pos2beam[make_pair(i,j)]};
						 moves.push_back(move);
					 }

					 intensity=s->prev_intensity(i,j);
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
  a_list = getShuffledApertureNeighbors(P);
  final_list.insert(final_list.end(), a_list.begin(), a_list.end());

  std::random_shuffle(final_list.begin(), final_list.end());


  return(final_list);
};

vector < NeighborMove > IntensityILS2::getShuffledNeighbors(Plan &P, int k, bool intensity) {

  vector< NeighborMove > a_list, i_list, final_list;

  i_list = getShuffledIntensityNeighbors(P);
  a_list = getShuffledApertureNeighbors(P);

  int size_i= (double)i_list.size()/k + 0.5;
  int size_a= (double)a_list.size()/k + 0.5;

  int i=0, a=0;
  while(!i_list.empty() || !a_list.empty()){

	  if(intensity){
		  if(!i_list.empty()) {
			  final_list.push_back( i_list.back() );
		  	  i_list.pop_back();
		  }
		  i++;
	  }else{
		  if(!a_list.empty()) {
			  final_list.push_back( a_list.back() );
		  	  a_list.pop_back();
		  }
		  a++;
	  }

	  if(i==size_i || a==size_a){
		  i=0; a=0;
		  intensity = !intensity;
	  }

  }

  return(final_list);
};

vector < NeighborMove> IntensityILS2::getNeighborhood(Plan& current_plan,
                                       NeighborhoodType ls_neighborhood,
                                       LSTarget ls_target){
  vector < NeighborMove> neighborList;

  if (ls_neighborhood == intensity) {
    neighborList = getShuffledIntensityNeighbors(current_plan);
  } else if (ls_neighborhood == aperture) {
    if(ls_target.target_type!=LSTargetType::target_friends)
       neighborList = getShuffledApertureNeighbors(current_plan);
    else
       neighborList = getShuffledApertureNeighbors_target(current_plan,vsize);
  }else if (ls_neighborhood == mixed) {
    //mixed
    neighborList = getShuffledNeighbors(current_plan);
  }else if (ls_neighborhood == smixed_i) {
	    neighborList = getShuffledNeighbors(current_plan, 2, true);
  }else if (ls_neighborhood == smixed_a) {
	    neighborList = getShuffledNeighbors(current_plan, 2, false);
  }

  if(ls_target.target_type == target_friends)
     prioritizeFriends(neighborList, ls_target, current_plan);

  return neighborList;
}

void IntensityILS2::prioritizeFriends(vector < NeighborMove >& neighborList, LSTarget ls_target, Plan& current_plan){

  int k=0;
  for(int i=0; i<neighborList.size();i++){
    NeighborMove n = neighborList[i];

    if(n.type==1){ //aperture
       //if target has an intensity move the aperture moves
       //in the same station AND APERTURE are  priorized in the neighbourhood
       if(ls_target.target_move.type==1 &&
         ls_target.target_move.station_id == n.station_id ){
    	   Station* s=current_plan.get_station(n.station_id);
    	   pair<int, int> p1=s->beam2pos[n.beamlet_id];
    	   if (ls_target.target_move.aperture_id+ls_target.target_move.action == int(s->I(p1.first,p1.second)+0.5)){
             neighborList[i]=neighborList[k];
             neighborList[k]=n;
             k++;
    	   }
       }

       //if target has an aperture move, friend beamlets
       //are priorized in the neighbourhood
       if(ls_target.target_move.type==2 &&
         ls_target.target_move.station_id == n.station_id){
            Station* s=current_plan.get_station(n.station_id);
            pair<int, int> p1=s->beam2pos[ls_target.target_move.beamlet_id];
            pair<int, int> p2=s->beam2pos[n.beamlet_id];
            if(max (abs(p1.first-p2.first), abs(p1.second-p2.second)) <=1  ){
              neighborList[i]=neighborList[k];
              neighborList[k]=n;
              k++;
            }
       }
     }

    }
}

/*
struct NeighborMove {
  // Type
  // 1: aperture
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


list<pair<int, double> > IntensityILS2::get_changes_in_fm(Plan &current_plan, NeighborMove move){
  int type            = move.type;
  Station* s          = current_plan.get_station(move.station_id);
  int beamlet         = move.beamlet_id;
  int action          = move.action;
  double intens	      = move.aperture_id;

  int i = s->beam2pos[beamlet].first;
  int j = s->beam2pos[beamlet].second;

  list<pair<int, double> > changes;

  if (type == 1) {
    // single beam
    double intensity;

    if (action < 0) // decrease intensity
      intensity=s->prev_intensity(i,j);
     else // Increase intensity
      intensity=s->next_intensity(i,j);
     if(intensity == s->I(i,j)) return changes; //the move is not valid
     changes.push_back(make_pair(s->pos2beam[make_pair(i,j)], intensity-s->I(i,j)));

  } else {
    // change intensity
    if((int)intens > s->int2nb.size() || (intens == 0.0 && s->int2nb.size()==s->getNbApertures()) )
        return changes;
    
    if(intens != 0.0){
       map< int, int >::iterator iter=s->int2nb.begin();
       std::advance( iter, ((int) intens-1) );
       intens = iter->first;
     }

    if (action < 0)  changes = s->get_changes_intensity_move(intens, -1.0);
    else changes = s->get_changes_intensity_move(intens, +1.0);

  }

  return changes;   

}

double IntensityILS2::applyMove (Plan & current_plan, NeighborMove move, bool p) {

  int type            = move.type;
  Station* s          = current_plan.get_station(move.station_id);
  int beamlet         = move.beamlet_id;
  int action          = move.action;
  double intens	  = move.aperture_id;
  

  int i = s->beam2pos[beamlet].first;
  int j = s->beam2pos[beamlet].second;

  double delta_eval= -1.0;

  if (type == 1) {
    // single beam
    double intensity;

    if (action < 0) // decrease intensity
      intensity=s->prev_intensity(i,j);
     else // Increase intensity
      intensity=s->next_intensity(i,j);

     s->change_intensity(i, j, intensity);

  } else {
    // change intensity
    if(intens != 0.0){
       map< int, int >::iterator iter=s->int2nb.begin();
       std::advance( iter, ((int) intens-1) );
       intens = iter->first;
     }

    if (action < 0)  s->change_intensity(intens, -1.0);
    else s->change_intensity(intens, +1.0);

  }


  return(0.0);
}

string IntensityILS2::planToString(Plan &P) {
  return(P.toStringIntensities());
};

}
