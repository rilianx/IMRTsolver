/*
 * ApertureILS.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "ApertureILS.h"

namespace imrt {

/* bsize: number of beamlets to be used in the target beamlet heuristic
   vsize: number of voxels to be considered when selecting the targeted beamlets
   prob_intensity: probability to perform local search over intensity
   step_intensity: step size for intensity
   do_perturbate: boolean variable that indicates if perturbation must be performed
   acceptance: type of acceptnace criterion to be used
*/
ApertureILS::ApertureILS(int bsize, int vsize, double _prob_intensity,
                         int _step_intensity):
                         ILS(bsize, vsize, 0),
                         prob_intensity(_prob_intensity),
                         step_intensity(_step_intensity) {
};

ApertureILS::ApertureILS(const ApertureILS & ils): ILS(ils) {
  prob_intensity=ils.prob_intensity;
  step_intensity=ils.step_intensity;
};

bool ApertureILS::isBeamletModifiable (int beamlet, Station* station, bool open_flag) {

  if (!station->isActiveBeamlet(beamlet)) return(false);
  if (!open_flag) {
    if (station->anyOpen(beamlet) && search_aperture)
      return(true);
    if (station->canReduceIntensity(beamlet) && search_intensity)
      return(true);
  } else {
    if (station->anyClosed(beamlet) && search_aperture)
      return(true);
    if (station->canIncreaseIntensity(beamlet) && search_intensity)
      return(true);
  }
  return(false);
};

// side=1 left
// side=2 right
// other will check both
double ApertureILS::closeBeamlet(int beamlet, int side, int aperture,
                                 Station& station, double c_eval, Plan& P) {
  double aux_eval, l_eval, r_eval;
  list <pair<int, double> > diff;

  if (P.get_delta_eval(station, beamlet, -station.getApertureIntensity(aperture)) > c_eval){
    P.clearLast();
    //station.clearHistory();
    return(c_eval);
  }

  if (side == 1) {
    diff = station.closeBeamlet(beamlet, aperture, true);
    if (diff.size() > 0) {
      c_eval = P.incremental_eval(station, diff);
    } else {
      P.clearLast();
    }
  } else if (side==2) {
    diff = station.closeBeamlet(beamlet, aperture, false);
    if (diff.size() > 0) {
      c_eval = P.incremental_eval(station, diff);
    } else {
      P.clearLast();
    }
  } else {
    l_eval = c_eval;
    diff = station.closeBeamlet(beamlet, aperture, true);
    if (diff.size() > 0) {
      l_eval = P.incremental_eval(station, diff);
      diff = station.undoLast();
      aux_eval = P.incremental_eval(station, diff);
    } else {
      P.clearLast();
    }
    r_eval = c_eval;
    diff = station.closeBeamlet(beamlet, aperture, false);
    if (diff.size() > 0) {
      r_eval = P.incremental_eval(station, diff);
      diff = station.undoLast();
      aux_eval = P.incremental_eval(station, diff);
    } else {
      P.clearLast();
    }

    if (r_eval > l_eval){
      diff = station.closeBeamlet(beamlet, aperture, true);
      c_eval = P.incremental_eval(station, diff);
    } else {
      diff = station.closeBeamlet(beamlet, aperture, false);
      c_eval = P.incremental_eval(station, diff);
    }
  }
  return(c_eval);
};


double ApertureILS::openBeamlet(int beamlet, int aperture, Station& station, double c_eval, Plan& P) {
  double aux_eval=0;
  list<pair<int, double> > diff;

  if (P.get_delta_eval(station, beamlet, station.getApertureIntensity(aperture)) > c_eval){
    P.clearLast();
    //station.clearHistory();
    return(c_eval);
  }

  diff = station.openBeamlet(beamlet, aperture);
  if(diff.size() == 0) {
    P.clearLast();
    return(c_eval);
  }
  c_eval = P.incremental_eval(station, diff);

  return(c_eval);
}

bool ApertureILS::perturbate(int no_improvement, int iteration) {
  if (!do_perturbate) return(false);
  //if (no_improvement >= ((double) iteration)*0.3) {
  //  return(true);
  //}
  if (no_improvement==100) return(true);
    return(false);
};

vector < NeighborMove > ApertureILS::getShuffledIntensityNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  vector< NeighborMove > a_list;
  list<Station*>::iterator s;
  s = stations.begin();

  for (int i = 0; i < stations.size(); i++) {
    for (int j = 0; j < (*s)->getNbApertures();j++){
      //One pair -j (-aperture) for reducing intensity
      //One pair for increasing intensity (+j)
      //Comment the ifs to generate all neighbors
      //if ((*s)->getApertureIntensity(j) > 0)
        a_list.push_back({1,i,j,-1,0});

      //if ((*s)->getApertureIntensity(j) <= (*s)->getMaxIntensity())
        a_list.push_back({1,i,j,1,0});
    }
    std::advance(s,1);
  }
  //cout << "Size neighborhood " << a_list.size()<< endl;
  std::random_shuffle(a_list.begin(), a_list.end());
  return(a_list);
};

pair<vector < NeighborMove >, vector<NeighborMove>> ApertureILS::getFriendsIntensityNeighbors(Plan &P, NeighborMove target){
  list<Station*> stations = P.get_stations();
  vector< NeighborMove > a_list;
  vector< NeighborMove > b_list;
  pair<vector <NeighborMove>, vector<NeighborMove>> all_list;
  bool target_flag = false;
  list<Station*>::iterator s;
  s = stations.begin();

  for (int i = 0; i < stations.size(); i++) {
    for (int j = 0; j < (*s)->getNbApertures();j++) {
      if (target.type == 1) {
        if (target.station_id == i && target.aperture_id == j)
            target_flag = true;
        else
            target_flag = false;
      } else {
        if (target.station_id ==i && target.aperture_id ==j)
            target_flag = true;
        else
            target_flag = false;
      }
      //One pair -j (-aperture) for reducing intensity
      //One pair for increasing intensity (+j)
      if ((*s)->getApertureIntensity(j) > 0){
        if (target_flag)
          a_list.push_back({1,i,j,-1,0});
        else
          b_list.push_back({1,i,j,-1,0});
      }

      if ((*s)->getApertureIntensity(j) <= (*s)->getMaxIntensity()) {
        if (target_flag)
          a_list.push_back({1,i,j,1,0});
        else
          b_list.push_back({1,i,j,1,0});
      }
    }
    std::advance(s,1);
  }
  //cout << "Size neighborhood " << a_list.size()<< endl;
  std::random_shuffle(a_list.begin(), a_list.end());
  std::random_shuffle(b_list.begin(), b_list.end());
  all_list = make_pair(a_list, b_list);
  return(all_list);
};

/*vector < NeighborMove > ApertureILS::getShuffledApertureNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  vector<NeighborMove> a_list;
  list<Station*>::iterator st;
  int beamlet;
  st = stations.begin();

  pair <int,int> pattern;
  pair <int,int> active;

  for (int s = 0; s < stations.size(); s++) {
    for (int a = 0; a < (*st)->getNbApertures(); a++){
      for (int k = 0; k< (*st)->collimator.getXdim() ; k++){
        //One pair -k (-row) for closing aperture
        //One pair (k) (row) for opening aperture
        pattern = (*st)->getApertureShape(a, k);
        active = (*st)->collimator.getActiveRange(k, (*st)->getAngle());

        // this row is not active for this station
        if (active.first == -1) continue;

        if (pattern.first == -1) {
          // all beamlets are closed in this row so we test
          // opening all beamlets
          for (int i=active.first; i<=active.second; i++){
            beamlet = (*st)->getBeamIndex(make_pair(k,i));
            a_list.push_back({2,s,a,1,beamlet});
          }
        } else {
          // open slid, try closing or opening
          beamlet = (*st)->getBeamIndex(make_pair(k,pattern.first));
          if((*st)->isOpenBeamlet(beamlet, a))
            a_list.push_back({2,s,a,-1,beamlet});
          else
            a_list.push_back({2,s,a,1,beamlet});

          if (pattern.first != pattern.second) {
            beamlet = (*st)->getBeamIndex(make_pair(k,pattern.second));
            if((*st)->isOpenBeamlet(beamlet, a))
              a_list.push_back({2,s,a,-2,beamlet});
            else
              a_list.push_back({2,s,a,1,beamlet});
          }
        }
      }
    }
    std::advance(st,1);
  }
  std::random_shuffle(a_list.begin(), a_list.end());
  return(a_list);
};
*/

vector < NeighborMove > ApertureILS::getShuffledApertureNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  vector<NeighborMove> a_list;
  list<Station*>::iterator st;
  st = stations.begin();

  pair <int,int> active;

  for (int s = 0; s < stations.size(); s++) {
    for (int a = 0; a < (*st)->getNbApertures(); a++){
      for (int k = 0; k< (*st)->collimator.getXdim() ; k++){
        active = (*st)->collimator.getActiveRange(k, (*st)->getAngle());

        // this row is not active for this station
        if (active.first == -1) continue;
        a_list.push_back({2,s,a,1,k});
        a_list.push_back({2,s,a,2,k});
        a_list.push_back({2,s,a,-1,k});
        a_list.push_back({2,s,a,-2,k});
      }
    }
    std::advance(st,1);
  }
  std::random_shuffle(a_list.begin(), a_list.end());
  return(a_list);
};


/* Order neighbors in two lists: a_list, b_list. The first
   list is for the preferred neighbors and the rest is in b_list
   - If previous move was an aperture move: schedule in a_list
   moves with the same action (opening or closing) on the same
   station+aperture and only neighbors beamlets to the previous
   move beamlet are considered.
   - If previous move was an intensity move: schedule in a_list
   moves related to the same station+aperture pair.
*/
pair <vector < NeighborMove >, vector < NeighborMove >>
      ApertureILS::getFriendsApertureNeighbors (Plan &P, NeighborMove target){
  list<Station*> stations = P.get_stations();
  vector<NeighborMove> a_list;
  vector<NeighborMove> b_list;
  pair <vector<NeighborMove>, vector<NeighborMove>> all_list;

  list<Station*>::iterator st;
  int beamlet;
  st = stations.begin();

  pair <int,int> pattern;
  pair <int,int> active;

  bool target_flag=false, staap_flag=false, beam_flag=false;

  pair <int,int> target_beamlet;
  if (target.type == 2)
    target_beamlet = P.get_station(target.station_id)->getPos(target.beamlet_id);

  //cout << "Target move: " << target.type << ":" << target.station_id << "," << target.aperture_id << "," << target.beamlet_id <<","  << target.action <<endl;

  for (int s = 0; s < stations.size(); s++) {
    for (int a = 0; a < (*st)->getNbApertures(); a++){
      // Signal a_list neighbors when intensity move is the target
      if (s==target.station_id && a==target.aperture_id)
        staap_flag =true;
      else
        staap_flag =false;

      for (int k = 0; k< (*st)->collimator.getXdim() ; k++){
        //One pair -k (-row) for closing aperture
        //One pair (k) (row) for opening aperture
        pattern = (*st)->getApertureShape(a, k);
        active = (*st)->collimator.getActiveRange(k, (*st)->getAngle());

        // this row is not active for this station
        if (active.first == -1) continue;

        if (pattern.first == -1) {
          // all beamlets are closed in this row so we test
          // opening all beamlets
          for (int i=active.first; i<=active.second; i++){
            beamlet = (*st)->getBeamIndex(make_pair(k,i));
            beam_flag = target_flag = true;
            if (target.type ==2) {
              if (staap_flag && abs(k-target_beamlet.first) <=1
                  && abs(i-target_beamlet.second) <=1) {
                beam_flag = true;
                if (target.action > 0) target_flag=true;
                else target_flag = false;
              } else {
                beam_flag = false;
              }
            }

            if (staap_flag && beam_flag && target_flag)
              a_list.push_back({2,s,a,1,beamlet});
            else
              b_list.push_back({2,s,a,1,beamlet});
          }
        } else {
          // open slid, try closing or opening
          beamlet = (*st)->getBeamIndex(make_pair(k,pattern.first));

          beam_flag = target_flag = true;
          if (target.type ==2) {
            if (staap_flag && abs(k-target_beamlet.first) <=1
                && abs(pattern.first-target_beamlet.second) <=1) {
              beam_flag = true;
            } else {
              beam_flag = false;
            }
          }

          if ((*st)->isOpenBeamlet(beamlet, a)) {
            if (target.type ==2) {
              if (target.action < 0) target_flag = true;
              else target_flag = false;
            }

            if (staap_flag && beam_flag && target_flag)
              a_list.push_back({2,s,a,-1,beamlet});
            else
              b_list.push_back({2,s,a,-1,beamlet});
          } else {
            if (target.type ==2) {
              if (target.action > 0) target_flag = true;
              else target_flag = false;
            }

            if (staap_flag && beam_flag && target_flag)
              a_list.push_back({2,s,a,1,beamlet});
            else
              b_list.push_back({2,s,a,1,beamlet});
          }

          if (pattern.first != pattern.second) {
            beamlet = (*st)->getBeamIndex(make_pair(k,pattern.second));
            beam_flag = target_flag = true;
            if (target.type ==2) {
              if (staap_flag && abs(k-target_beamlet.first) <=1
                  && abs(pattern.second-target_beamlet.second) <=1) {
                beam_flag = true;
              } else {
                beam_flag = false;
              }
            }

            if((*st)->isOpenBeamlet(beamlet, a)) {
              if (target.type ==2) {
                if (target.action < 0) target_flag = true;
                else target_flag = false;
              }

              if (staap_flag && beam_flag && target_flag)
                a_list.push_back({2,s,a,-2,beamlet});
              else
                b_list.push_back({2,s,a,-2,beamlet});
            } else {
              if (target.type ==2) {
                if (target.action > 0) target_flag = true;
                else target_flag = false;
              }

              if (staap_flag && beam_flag && target_flag)
                a_list.push_back({2,s,a,1,beamlet});
              else
                b_list.push_back({2,s,a,1,beamlet});
            }
          }
        }
      }
    }
    std::advance(st,1);
  }

  std::random_shuffle(a_list.begin(), a_list.end());
  std::random_shuffle(b_list.begin(), b_list.end());
  all_list = make_pair(a_list, b_list);
  return(all_list);
};

vector < NeighborMove > ApertureILS::getShuffledNeighbors(Plan &P) {

  vector< NeighborMove > a_list, final_list;

  final_list = getShuffledIntensityNeighbors(P);
  a_list = getShuffledApertureNeighbors(P);
  final_list.insert(final_list.end(), a_list.begin(), a_list.end());

  std::random_shuffle(final_list.begin(), final_list.end());

  return(final_list);
};

vector < NeighborMove > ApertureILS::getFriendsNeighbors(Plan &P, NeighborMove target) {

  pair <vector< NeighborMove >, vector<NeighborMove>> a_list, i_list;
  vector <NeighborMove> aux_list;
  vector <NeighborMove> final_list;

  i_list = getFriendsIntensityNeighbors(P, target);
  a_list = getFriendsApertureNeighbors(P, target);
  final_list.insert(final_list.end(), a_list.first.begin(), a_list.first.end());
  final_list.insert(final_list.end(), i_list.first.begin(), i_list.first.end());
  std::random_shuffle(final_list.begin(), final_list.end());

  aux_list.insert(aux_list.end(), a_list.second.begin(), a_list.second.end());
  aux_list.insert(aux_list.end(), i_list.second.begin(), i_list.second.end());
  std::random_shuffle(aux_list.begin(), aux_list.end());

  final_list.insert(final_list.end(), aux_list.begin(), aux_list.end());

  return(final_list);
};

vector < NeighborMove> ApertureILS::getNeighborhood(Plan& current_plan,
                                       NeighborhoodType ls_neighborhood,
                                       LSTarget ls_target){
  vector < NeighborMove> neighborList;
  pair <vector<NeighborMove>, vector<NeighborMove>> aux;

  if (ls_neighborhood == intensity) {
    if (ls_target.target_type == LSTargetType::target_none )
      neighborList = getShuffledIntensityNeighbors(current_plan);
    else {
      aux = getFriendsIntensityNeighbors(current_plan, ls_target.target_move);
      neighborList.insert(neighborList.end(), aux.first.begin(), aux.first.end());
      neighborList.insert(neighborList.end(), aux.second.begin(), aux.second.end());
    }
  } else if (ls_neighborhood == aperture) {
    if (ls_target.target_type == LSTargetType::target_none )
      neighborList = getShuffledApertureNeighbors(current_plan);
    else {
      aux = getFriendsApertureNeighbors(current_plan, ls_target.target_move);
      neighborList.insert(neighborList.end(), aux.first.begin(), aux.first.end());
      neighborList.insert(neighborList.end(), aux.second.begin(), aux.second.end());
    }
  } else {
    //mixed
    if (ls_target.target_type == LSTargetType::target_none )
      neighborList = getShuffledNeighbors(current_plan);
    else {
      neighborList = getFriendsNeighbors(current_plan, ls_target.target_move);

    }
  }
  return(neighborList);
}

/*double ApertureILS::applyMove (Plan & current_plan, NeighborMove move) {
  double current_eval = current_plan.getEvaluation();
  double aux_eval = current_plan.getEvaluation();
  list<pair<int, double> > diff;
  int type            = move.type;
  int station_id      = move.station_id;
  Station * s         = current_plan.get_station(move.station_id);
  int aperture        = move.aperture_id;
  int beamlet         = move.beamlet_id;
  int action          = move.action;


  if (!checkMove(current_plan, move))
    return (-1);

  if (type == 1) {
    // Intensity move
    if (action < 0) {
      // Reduce intensity
      aux_eval = current_plan.modifyIntensityAperture(station_id, aperture, -step_intensity);
    } else {
      // Increase intensity
      aux_eval = current_plan.modifyIntensityAperture(station_id, aperture, step_intensity);
    }
  } else {
    // Aperture move
    if (action < 0) {
      aux_eval = current_plan.closeBeamlet(station_id, aperture, beamlet, abs(action));
    } else {
      aux_eval = current_plan.openBeamlet(station_id, aperture, beamlet);
    }
  }
  current_eval=aux_eval;
  return(current_eval);

};*/

double ApertureILS::get_delta_eval(Plan &current_plan, NeighborMove move, list<pair<int, double> >& diff){
  double delta_eval = 0.0;
  int type            = move.type;
  int station_id      = move.station_id;
  Station * s         = current_plan.get_station(move.station_id);
  int aperture        = move.aperture_id;
  int row             = move.beamlet_id;
  int beamlet         = move.beamlet_id;
  int action          = move.action;

  if (!checkMove(current_plan, move))
    return (-1);

    if (type == 1) {
      // Intensity move
      if (action < 0) {
        // Reduce intensity
        delta_eval = current_plan.modifyIntensityAperture(station_id, aperture, -step_intensity, true, &diff);
      } else {
        // Increase intensity
        delta_eval = current_plan.modifyIntensityAperture(station_id, aperture, step_intensity, true, &diff);
      }
    } else {
      // Aperture move
      if (action < 0) {
        if (action == -1)
          delta_eval = current_plan.closeRow(station_id, aperture, row, true, true, &diff);
        else
          delta_eval = current_plan.closeRow(station_id, aperture, row, false, true, &diff);
      } else {
        if (action == 1)
          delta_eval = current_plan.openRow(station_id, aperture, row, true, true, &diff);
        else
          delta_eval = current_plan.openRow(station_id, aperture, row, false, true, &diff);
      }
    }

  return delta_eval;
}

double ApertureILS::applyMove (Plan & current_plan, NeighborMove move) {
  list<pair<int, double> > diff;
  double delta_eval = get_delta_eval(current_plan, move, diff);
  double current_eval = 0.0;
  Station *s = current_plan.get_station(move.station_id);

  current_eval = current_plan.incremental_eval(*s, diff);


  return(current_eval);
}


bool ApertureILS::checkMove (Plan & current_plan, NeighborMove move) {
  int type            = move.type;
  int station_id      = move.station_id;
  Station * s         = current_plan.get_station(move.station_id);
  int aperture        = move.aperture_id;
  int row             = move.beamlet_id;
  int action          = move.action;
  int beamlet;
  pair<int,int> pattern;
  pair<int,int> active;

  if (type ==1){
    // Intensity move
    if (action < 0) {
      if (s->getApertureIntensity(aperture) <= 0)
        return(false);
    } else {
      if (s->getApertureIntensity(aperture) >= s->getMaxIntensity())
        return(false);
    }
  } else {
    // Aperture move
    active = s->collimator.getActiveRange(row, s->getAngle());
    if (active.first == -1) return(false);

    pattern = s->getApertureShape(aperture, row);
    if (abs(action)==1)
      beamlet = pattern.first;
    else
      beamlet = pattern.second;

    if (action < 0) {
      // Close a beamlet in the row
      if (pattern.first==-1 || pattern.second==-1) return(false);
      if(!s->isActiveBeamlet(beamlet) || !s->isOpenBeamlet(beamlet, aperture))
        return(false);
    } else {
      // Open a beamlet in the row
      if (pattern.first==-1 || pattern.second==-1) return(true);
      if (action ==1){
        // Open to the left
        beamlet--;
        if (active.first > beamlet)
          return(false);
      } else {
        // Open to the right
        beamlet++;
        if (active.second < beamlet) return (false);
      }
      //if (!s->isActiveBeamlet(beamlet) && s->isOpenBeamlet(beamlet, aperture))
        //return(false);
    }
  }
  return(true);

};

string ApertureILS::planToString(Plan &P) {
  return(P.toStringApertures());
};


// This function performs a local search over all the
// aperture intensities in a treatment plan.
double ApertureILS::iLocalSearch (Plan& P,  double max_time, bool verbose) {
    list<Station*> stations = P.get_stations();
    list<Station*>::iterator s;

    std::clock_t time_end, time_begin;
    double used_time;

    double local_best_eval, current_eval, aux_eval;
    list<pair<int, double> > diff;
    vector<NeighborMove> a_list;
    NeighborMove tabu;

    bool improvement = true;
    bool best_improvement=false;
    bool completed = false;
    int i, j, best_n;

    tabu = {0,0,0,0,0};
    best_n = -1;
    j=-1;
    local_best_eval = current_eval = aux_eval = P.getEvaluation();

    if (verbose)
      cout << "Staring intensity local search..." << endl;
    time_begin=clock();

    // Main local search loop
    while (improvement) {
      improvement = false;
      completed = false;
      j++;
      current_eval = local_best_eval;
      a_list = getShuffledIntensityNeighbors(P);

      if (verbose) {
        cout << "  iLS Neighborhood "<< j << " size "<< a_list.size()
             << "    current " << local_best_eval << endl;
      }

      // Check all the neighbors
      for (i = 0; i < a_list.size(); i++) {
         //skip the tabu neighbor (returns the station to previous state)
         if (a_list[i].station_id == tabu.station_id &&
             a_list[i].aperture_id == tabu.aperture_id &&
             a_list[i].action == tabu.action ) {
           if (i == (a_list.size()-1)) completed = true;
           continue;
         }
         //get the station of the movement
         s = stations.begin();
         std::advance(s, a_list[i].station_id);

         if (verbose)
           cout << "  iLS Neighbor " << i << " over station "
           << (*s)->getAngle() << " aperture " <<  a_list[i].aperture_id;

         aux_eval = current_eval;

         //apply step_size intensity change (-(a+1) or +(a+1))
         if (a_list[i].action < 0 ){
           diff = (*s)->getModifyIntensityApertureDiff(a_list[i].aperture_id,
                                                       -step_intensity);
           if (P.get_delta_eval((*(*s)), diff) > current_eval) {
             (*s)->clearHistory();
           } else {
             diff = (*s)->modifyIntensityAperture(a_list[i].aperture_id,
                                                  -step_intensity);
             if (diff.size() > 0) {
               aux_eval = P.incremental_eval(*(*s), diff);
             }
           }
           if (verbose)
             cout << " (-" << step_intensity << ")";
         } else {
           diff = (*s)->getModifyIntensityApertureDiff(a_list[i].aperture_id,
                                                       step_intensity);
           if (P.get_delta_eval((*(*s)), diff) > current_eval) {
             (*s)->clearHistory();
           } else {
             diff = (*s)->modifyIntensityAperture(a_list[i].aperture_id,
                                                  step_intensity);
             if (diff.size() > 0) {
               aux_eval = P.incremental_eval(*(*s), diff);
             }
           }
           if (verbose)
             cout << " (+"<< step_intensity << ")";
         }

         if (verbose )
             cout << " result " << aux_eval <<endl;
         // First improvement
         if ((local_best_eval - aux_eval) > 0.00001){
           local_best_eval = aux_eval;
           best_n = i;
           improvement = true;
           if (verbose)
             cout << "     improvement " << aux_eval << endl ;


           // Add a tabu movement since we found improvement
           tabu.type = a_list[i].type;
           tabu.station_id = a_list[i].station_id;
           tabu.aperture_id = a_list[i].aperture_id;
           tabu.action = a_list[i].action;

           // If first improvement has been chosen break and
           if (!best_improvement) {
             if (i==(a_list.size()-1)) completed = true;
             break;
           }
         }

         // Undo movement to continue the search
         diff = (*s)->undoLast();
         if (diff.size()>0)
           aux_eval = P.incremental_eval(*(*s), diff);

         if (i == (a_list.size()-1)) completed = true;

         time_end = clock();
         used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
         if (max_time!=0 && used_time >= max_time) {
           break;
         }
       }

       //Apply best neighbor
       if (improvement && best_improvement) {
         if (a_list[best_n].action < 0 ){
           diff = (*s)->modifyIntensityAperture(a_list[best_n].aperture_id,
                                                -step_intensity);
           aux_eval = P.incremental_eval(*(*s), diff);
         } else {
           diff = (*s)->modifyIntensityAperture(a_list[best_n].aperture_id,
                                                step_intensity);
           aux_eval = P.incremental_eval(*(*s), diff);
         }
       }

       time_end = clock();
       used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
       if (max_time!=0 && used_time >= max_time) {
         break;
       }
    }

    cout << "  iLS best: " << local_best_eval ;
    if (!completed) cout << ": [nolo] : ";
    else  cout << ": [lo] : ";
    time_end = clock();
    used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
    cout << max_time << " :" << used_time << endl;
    return(local_best_eval);
};

// This function performs a local search over all the
// aperture patterns in a treatment plan.
double ApertureILS::aLocalSearch(Plan& P, double max_time, bool verbose) {
  Station *s;
  std::clock_t time_end, time_begin;
  double used_time;
  double local_best_eval, aux_eval, current_eval;
  bool improvement=true;
  bool best_improvement=false;

  vector < NeighborMove > a_list;
  NeighborMove tabu;
  NeighborMove best_move;
  list<pair<int, double> > diff;
  int i,j, best_n, aperture, beamlet, sign;
  bool completed = false;
  best_n = -1;
  j=-1;
  local_best_eval = aux_eval= P.getEvaluation();
  if (verbose)
    cout << "Staring aperture local search..." << endl;

  time_begin=clock();

  while (improvement) {
    j++;
    a_list = getShuffledApertureNeighbors(P);
    improvement = false;
    completed = false;

    if (verbose)
      cout << "Neighborhood "<<j<< " size "<< a_list.size()
           << "    current " << local_best_eval << endl;

    current_eval = P.getEvaluation();
    for (i = 0; i < a_list.size(); i++) {

      s = P.get_station(a_list[i].station_id);
      aperture = a_list[i].aperture_id;
      beamlet = a_list[i].beamlet_id;
      sign = a_list[i].action;

      if (sign < 0) {
        aux_eval = closeBeamlet(beamlet, abs(sign), aperture, *s, current_eval, P);
      } else {
        aux_eval = openBeamlet(beamlet, aperture, *s, current_eval, P);
      }

      if (verbose) {
          cout << "  aLS Neighbor " << i << " over station "
               << s->getAngle() << " aperture " << aperture
               << " beamlet " << beamlet << " operator " << sign
               << " result " << aux_eval << endl;
      }

      if ((local_best_eval - aux_eval) > 0.00001){
        local_best_eval = aux_eval;
        best_n = i;
        improvement = true;
        if (verbose)
          cout << "     improvement " << aux_eval << endl ;
        if (!best_improvement) {
          if (i==(a_list.size()-1)) completed=true;
          i = a_list.size();
          break;
        }
      }

      diff = s->undoLast();
      if (diff.size() > 0)
        aux_eval = P.incremental_eval(*s, diff);

      if (i==(a_list.size()-1)) completed=true;

      time_end = clock();
      used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
      if (max_time!=0 && used_time >= max_time) {
        if (verbose)
          cout << "  aLS timed out " << endl;
        break;
      }

    }

    //Apply best neighbor
    if (improvement && best_improvement) {
      s = P.get_station(a_list[best_n].station_id);
      aperture = a_list[best_n].aperture_id;
      beamlet = a_list[best_n].beamlet_id;
      sign = a_list[best_n].action;
      if (sign < 0) {
        aux_eval = closeBeamlet(beamlet, abs(sign), aperture, *s, aux_eval, P);
      } else{
        aux_eval = openBeamlet(beamlet, aperture, *s, aux_eval, P);
      }
    }

    time_end = clock();
    used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
    if (max_time!=0 && used_time >= max_time) {
      if (verbose)
        cout << "  aLS timed out " << endl;
      break;
    }
  }

  cout << "  aLS best: " << local_best_eval ;
  if (!completed) cout << ": [nolo] : ";
  else  cout << ": [lo] : ";
  time_end = clock();
  used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
  cout << max_time << " :" << used_time << endl;

  return(local_best_eval);
};

double ApertureILS::simpleLocalSearch(Plan& P, bool verbose) {
  bool improvement=true;
  double local_best=P.getEvaluation(), aux_best;
  while (improvement) {
    improvement =false;
    aux_best = iLocalSearch(P, verbose);
    if (aux_best<local_best) {
      local_best=aux_best;
    }
    aux_best = aLocalSearch(P, verbose);
    if (aux_best<local_best) {
      local_best=aux_best;
      improvement=true;
    }
  }
  cout << "Local search finished with: " << local_best << " evaluation "
       << P.getEvaluation()<< endl;
  return(P.getEvaluation());
};

int ApertureILS::getStepIntensity () {
  return(step_intensity);
};

}
