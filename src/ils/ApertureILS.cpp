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

pair <bool, pair<Station*, int> > ApertureILS::getLSBeamlet(Plan& P){
  auto sb = P.best_beamlets(bsize, vsize);
  Station* s;
  if (sb.size() < 1) return(make_pair(false, make_pair(s,-1))); 

  auto it = sb.begin();
  int addit = rand()%sb.size();
  std::advance(it,addit);

  s = it->second.first; 
  int beamlet = it->second.second;
  bool sign = it->first.second; 
  
  if (!s->isActiveBeamlet(beamlet)){ 
    cout << "ERROR!: Beamlet "<< beamlet << " not active in station" 
         << s->getAngle()  << endl;
    return(make_pair(false, make_pair(s,-1)));
  }

  int i=0;
  while (!isBeamletModifiable(beamlet, s, !sign)) {
    // If no beamlet was found return -1
    i++;
    if (i==sb.size()) {
      return(make_pair(false, make_pair(s,-1)));
    }
    std::advance(it,1);
    if (it==sb.end()) it = sb.begin();

    s = it->second.first; 
    beamlet = it->second.second;
    sign = it->first.second; 
  }
  return(make_pair(sign, make_pair(s,beamlet)));
}

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

/*double ApertureILS::improvementIntensity(int beamlet, Station& station, 
                                         bool open_beamlet, 
                                         double c_eval, Plan & P, 
                                         bool best_improvement) {
  
  double local_best =-1, aux_eval=0;
  list<pair<int, double> > diff;
  int a, local_a, i;
  bool flag = true;
  vector<int> a_list;
  
  a_list = station.getOpen(beamlet);
  
  if (a_list.size()<1) {
    //cout << endl << "Warning: not possible to make intensity change beamlet:" << beamlet <<endl;
    cout << "[NOPT] ";
    return(c_eval);
  }
  
  a = (int)rand() % a_list.size();
  i=0;

  // Check every aperture 
  while(flag) {
    if (open_beamlet) {
      diff = station.modifyIntensityAperture(a_list[a], step_intensity);
      aux_eval = P.incremental_eval(station, diff);
      //cout << endl<<"  Increasing intensity aperture: " << a <<" eval: " << aux_eval << " list: ";
      //for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      //cout << endl;
    } else {
      diff = station.modifyIntensityAperture(a_list[a], -step_intensity);
      aux_eval = P.incremental_eval(station, diff);
      //cout << endl<< "  Reducing intensity aperture: " << a <<" eval: " << aux_eval <<  " list: ";
      //for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      //cout << endl;
    }
    
    // First improvement
    if (local_best < 0 || (local_best- aux_eval)>0.00001){
      local_best=aux_eval;
      local_a=a;
      i=0; // to support best improvement
      if (!best_improvement)
        return(local_best);
    } 
    
    //Undo change
    diff = station.undoLast();
    if (diff.size()>0)
      aux_eval = P.incremental_eval(station, diff);
    
    i++;
    a++;
    if (a==a_list.size()) a=0;
    if (i==a_list.size()) flag=false;
  }
  
  // If no improvement was found, return the best solution found.
  if (!open_beamlet) {
    diff = station.modifyIntensityAperture(a_list[local_a], 1);
    local_best = P.incremental_eval(station, diff);
  } else{
    diff = station.modifyIntensityAperture(a_list[local_a], -1);
    local_best = P.incremental_eval(station, diff);
  }
  return(local_best);
  }*/

double ApertureILS::openBeamlet(int beamlet, int aperture, Station& station, double c_eval, Plan& P) {
  double aux_eval=0;
  list<pair<int, double> > diff;
  
  if (P.get_delta_eval(station, beamlet, station.getApertureIntensity(aperture)) > c_eval){
    P.clearLast();
    //station.clearHistory();
    return(c_eval);
  }

  diff = station.openBeamlet(beamlet, aperture);
  if(diff.size() <1) {
    P.clearLast();
    return(c_eval);
  }
  c_eval = P.incremental_eval(station, diff);

  return(c_eval);
}


/*double ApertureILS::improvementAperture(int beamlet, Station& station, bool open_beamlet, 
                                             double c_eval, Plan& P, bool best_improvement) {
  double local_best=-1, aux_eval=0;
  list<pair<int, double> > diff;
  int a, local_a, i;
  bool flag=true; 
  vector<int> a_list;
  
  if (!open_beamlet) a_list = station.getOpen(beamlet);
  else a_list = station.getClosed(beamlet);
  
  if (a_list.size()<1) {
    cout << "[NOPT] ";
    //cout << endl << "Warning: not possible to make aperture change beamlet:" << beamlet <<endl;
    return(c_eval);
  }
  
  a = (int)rand() % a_list.size();
  i = 0;
  
  // Check every aperture 
  while(flag) {
    if (!open_beamlet) {
      aux_eval = closeBeamlet(beamlet, 0, a_list[a], station, c_eval, P); 
    } else {
      aux_eval = openBeamlet(beamlet, a_list[a], station, c_eval, P);
    } 
    
    // First improvement
    if (local_best < 0 || (local_best-aux_eval) >0.00001){
      local_best = aux_eval;
      local_a =a ;
      i=0; // to support best improvement
      if (!best_improvement)
        return(local_best);
    } 
    
    diff = station.undoLast();
    if (diff.size()>0)
      aux_eval = P.incremental_eval(station, diff);
    i++;
    a++;
    if (a==a_list.size()) a=0;
    if (i>=a_list.size()) flag=false;
  }
  
  // If no improvement was found, return the best solution found.
  if (!open_beamlet) {
    local_best = closeBeamlet(beamlet, 0, a_list[local_a], station, c_eval, P);
  } else{
    local_best = openBeamlet(beamlet, a_list[local_a], station, c_eval, P);
  }
  return(local_best);
  
  }*/




/* Local search procedure used in the HM2019*/
/*double ApertureILS::localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P) {
  double aux_eval, local_eval=P.getEvaluation();
  bool search_a;
  bool best_improvement=ls_type;
  
  int beamlet = target_beam.second.second;
  Station* s = target_beam.second.first;
  bool sign = target_beam.first;

  if (search_aperture && !search_intensity) {
    cout << ", ls aperture ";
    aux_eval = improvementAperture(beamlet, *s, !sign, local_eval, P, best_improvement);
    cout << ", found: " << aux_eval;
  } else if (!search_aperture && search_intensity) {
    cout << ", ls intensity ";
    aux_eval = improvementIntensity(beamlet, *s, !sign, local_eval, P, best_improvement);
    cout << ", found: " << aux_eval;
  } else if (search_aperture && search_intensity) {
    if (!sign)
      search_a = s->anyClosed(beamlet);
    else 
      search_a = s->anyOpen(beamlet);

    if (((double) rand() / (RAND_MAX)) <= prob_intensity || !search_a) {
      cout << ", ls intensity";
      aux_eval = improvementIntensity(beamlet, *s,  !sign, local_eval, P, best_improvement);
      cout << ", found: " << aux_eval;
    } else {
      cout << ", ls aperture";
      aux_eval = improvementAperture(beamlet, *s, !sign, local_eval, P, best_improvement);
      cout << ", found: " << aux_eval;
    }
  } else {
    cout << "ERROR: Search option not recognized!!!" << endl;
    
  }
  if (local_eval>aux_eval) { 
    tabu.clear();
  }else{
    tabu.push_back(make_pair(s, beamlet));
  }
  return(aux_eval);
}
*/
double ApertureILS::perturbation(Plan& P) {
  list<Station*> stations = P.get_stations();
  int beamlet, aperture;
  list<pair<int,double>> diff;
  double aux_eval=P.getEvaluation();
  
  cout << "##  Perturbation: " ;
  if (perturbation_size==0) cout << "none";   
  //else tabu.clear();
  for (int i=0; i<perturbation_size; i++) {
    list<Station*>::iterator s=stations.begin();
    std::advance(s,rand()%stations.size());
    if (((double) rand() / (RAND_MAX)) > 0.3) {
      //Modify intensity
      aperture = (rand()% (*s)->getNbApertures());
      if (((double) rand() / (RAND_MAX)) > 0.5){ 
        diff = (*s)->modifyIntensityAperture(aperture, -step_intensity);
      }else{ 
        diff = (*s)->modifyIntensityAperture(aperture, step_intensity);
      }
      aux_eval = P.incremental_eval(*(*s), diff);
      cout << "(" << (*s)->getAngle() << ","<< aperture<<") ";
    } else {
      //Modify aperture
      do { 
        beamlet = (rand()% (*s)->getNbBeamlets());
      } while (!(*s)->isActiveBeamlet(beamlet));
      aperture = (rand()% (*s)->getNbApertures());
      if ((*s)->isOpenBeamlet(beamlet, aperture)){
        if (((double) rand() / (RAND_MAX)) > 0.5) 
          diff = (*s)->closeBeamlet(beamlet, aperture, false);
        else
          diff = (*s)->closeBeamlet(beamlet, aperture, true);
      } else {
        diff = (*s)->openBeamlet(beamlet, aperture);
      }
      aux_eval = P.incremental_eval(*(*s), diff);
      cout << "(" << (*s)->getAngle() << "," << aperture<< "," << beamlet << ") ";
    }
    (*s)->clearHistory();
  }
  
  //aux_eval = P.eval();
  cout << ", new eval: "<< aux_eval<< endl;
  return(aux_eval);
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
      if ((*s)->getApertureIntensity(j) > 0)
        a_list.push_back({1,i,j,-1,0});
      
      if ((*s)->getApertureIntensity(j) <= (*s)->getMaxIntensity())
        a_list.push_back({1,i,j,1,0});
    }
    std::advance(s,1);
  }
  //cout << "Size neighborhood " << a_list.size()<< endl;
  std::random_shuffle(a_list.begin(), a_list.end());
  return(a_list);
};

vector < NeighborMove > ApertureILS::getShuffledApertureNeighbors(Plan &P){
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

vector < NeighborMove > ApertureILS::getOrderedApertureNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  
  pair<bool, pair<Station*, int>> target_beam = getBestLSBeamlet(P);
  Station * target_station = target_beam.second.first;
  int s_target;
  
  vector< NeighborMove > a_list, final_list;
  list<Station*>::iterator st;
  int beamlet;
  st = stations.begin();
  
  pair <int,int> pattern;
  pair <int,int> active;
  
  for (int s = 0; s < stations.size(); s++) {
    if ((*st)->getAngle() == target_station->getAngle()){
      s_target=s;
      continue;
    }
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
  
  // Add target station neighbors
  for (int a = 0; a < target_station->getNbApertures(); a++) {
    for (int k = 0; k< target_station->collimator.getXdim() ; k++) {
      active = target_station->collimator.getActiveRange(k, target_station->getAngle());
      if (active.first == -1) continue;
      pattern = target_station->getApertureShape(a, k);
      
      beamlet = target_station->getBeamIndex(make_pair(k,pattern.first));
      if (target_station->isOpenBeamlet(beamlet, a))
        final_list.push_back ({2, s_target, a, -1, beamlet});
        //final_list.push_back (make_pair(make_pair(s_target,a) , make_pair(-1,beamlet)));
      else 
        final_list.push_back ({2, s_target, a, 1, beamlet});
        //final_list.push_back (make_pair(make_pair(s_target,a) , make_pair( 1,beamlet)));
      
      if (pattern.first != pattern.second) {
        beamlet = target_station->getBeamIndex(make_pair(k,pattern.second));
        if (target_station->isOpenBeamlet(beamlet, a))
          final_list.push_back ({2, s_target, a, -2, beamlet});
         //final_list.push_back (make_pair(make_pair(s_target,a) , make_pair(-2,beamlet)));
        else 
          final_list.push_back ({2, s_target, a, 1, beamlet});
         // final_list.push_back(make_pair(make_pair(s_target,a) , make_pair( 1,beamlet)));
      }
    }
  }
  std::random_shuffle(final_list.begin(), final_list.end());
  for (int i=0; i<a_list.size(); i++) {
    final_list.push_back(a_list[i]);
  }
  
  return(final_list);
};

vector < NeighborMove > ApertureILS::getShuffledNeighbors(Plan &P) {

  vector< NeighborMove > a_list, final_list;

  final_list = getShuffledIntensityNeighbors(P);
  a_list = getShuffledApertureNeighbors(P);
  final_list.insert(final_list.end(), a_list.begin(), a_list.end());

  std::random_shuffle(final_list.begin(), final_list.end());

  return(final_list);
};

vector < NeighborMove> ApertureILS::getNeighborhood(Plan& current_plan, 
                                       NeighborhoodType ls_neighborhood, 
                                       LSTarget ls_target){
  vector < NeighborMove> neighborList;

  if (ls_neighborhood == intensity) {
    neighborList = getShuffledIntensityNeighbors(current_plan);
  } else if (ls_neighborhood == aperture) {
    neighborList = getShuffledApertureNeighbors(current_plan);
  } else {
    //mixed
    neighborList = getShuffledNeighbors(current_plan);
  }
  return(neighborList);
}

/*double ApertureILS::applyMove (Plan & current_plan, NeighborMove move) {
  double current_eval = current_plan.getEvaluation();
  double aux_eval = current_plan.getEvaluation();
  list<pair<int, double> > diff;
  int type            = move.type;
  Station * s         = current_plan.get_station(move.station_id);
  int aperture        = move.aperture_id;
  int beamlet         = move.beamlet_id;
  int action          = move.action;

  cout << "  move type " << move.type << ", s:" << 
          move.station_id << ", a:" << aperture << 
          ", b:" << beamlet << ", e:"<< action << endl;
  if (type == 1) {
    // Intensity move
    if (action < 0) {
      // Reduce intensity
      diff = s->modifyIntensityAperture(aperture, 
                                        -step_intensity);
      if (diff.size() > 0) {
        aux_eval = current_plan.incremental_eval(*s, diff);
      } else {
        current_plan.clearLast();
      }
    } else {
      // Increase intensity
      diff = s->modifyIntensityAperture(aperture, 
                                        step_intensity);
      if (diff.size() > 0) {
        aux_eval = current_plan.incremental_eval(*s, diff);
      } else {
        current_plan.clearLast();
      }
    }
  } else {
    // Aperture move
    if (action < 0) {
      aux_eval = closeBeamlet(beamlet, abs(action), aperture, *s, current_eval, current_plan); 
    } else {
      aux_eval = openBeamlet(beamlet, aperture, *s, current_eval, current_plan);
    }
  }
  current_eval=aux_eval;
  return(current_eval);
            
};*/

double ApertureILS::applyMove (Plan & current_plan, NeighborMove move) {
  double current_eval = current_plan.getEvaluation();
  double aux_eval = current_plan.getEvaluation();
  list<pair<int, double> > diff;
  int type            = move.type;
  int station_id      = move.station_id;
  Station * s         = current_plan.get_station(move.station_id);
  int aperture        = move.aperture_id;
  int beamlet         = move.beamlet_id;
  int action          = move.action;

  /*cout << "  move type " << move.type << ", s:" << 
          move.station_id << ", a:" << aperture << 
          ", b:" << beamlet << ", e:"<< action << endl;*/
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



