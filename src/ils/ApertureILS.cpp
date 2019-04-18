/*
 * ApertureILS.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "ApertureILS.h"

namespace imrt {


ApertureILS::ApertureILS(int bsize, int vsize, bool search_intensity, bool search_aperture, 
                         double prob_intensity, int step_intensity , double initial_temperature,  
                         double alpha, bool do_perturbate, int perturbation_size, int acceptance=0, 
                         int ls_type=ApertureILS::FIRST_IMPROVEMENT): 
                         ILS(bsize, vsize, acceptance), search_intensity(search_intensity), 
                         search_aperture(search_aperture), prob_intensity(prob_intensity), 
                         step_intensity(step_intensity) , initial_temperature(initial_temperature), 
                         alpha(alpha), do_perturbate(do_perturbate), perturbation_size(perturbation_size), 
                         ls_type(ls_type){
//  cout << "per:" << perturbation_size << endl;
  temperature=initial_temperature;
};

ApertureILS::ApertureILS(const ApertureILS & ils): ILS(ils) {
  search_intensity=ils.search_intensity;
  search_aperture=ils.search_aperture;
  prob_intensity=ils.prob_intensity;
  step_intensity=ils.step_intensity;
  initial_temperature=ils.initial_temperature;
  alpha=ils.alpha;
  do_perturbate=ils.do_perturbate;
  perturbation_size=ils.perturbation_size;
  ls_type=ils.ls_type;
  temperature=ils.temperature;
};

bool ApertureILS::isBeamletModifiable(int beamlet, Station* station, bool open_flag) {
  for (auto it=tabu.begin(); it!=tabu.end();it++){
    if (it->first->getAngle() == station->getAngle() && it->second == beamlet) return(false);
  }
  
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
}

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
    cout<< "ERROR!: Beamlet "<< beamlet << " not active in station" << s->getAngle()  << endl;
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

double ApertureILS::closeBeamlet(int beamlet, int aperture, Station& station, double c_eval, Plan& P) {
  double aux_eval, t_eval, f_eval;
  list<pair<int, double> > diff;
  
  diff = station.closeBeamlet(beamlet, aperture, true);
  if (diff.size()>=1) {
    c_eval = t_eval = P.incremental_eval(station, diff);
    //cout << "  Closing left eval: " << c_eval << " list: ";
    //for (list<pair<int,double>>::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
    //cout << endl;
    diff = station.undoLast();
    aux_eval = P.incremental_eval(station, diff);
  }
  
  diff = station.closeBeamlet(beamlet, aperture, false);
  if (diff.size()>=1) {
    f_eval = P.incremental_eval(station, diff);
    //cout << "  Closing right eval: " << f_eval << " list: "; 
    //for (list<pair<int,double>>::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
    //cout << endl;
    
    if ( f_eval > t_eval ) {
      diff = station.undoLast();   
      aux_eval = P.incremental_eval(station, diff);
      diff = station.closeBeamlet(beamlet, aperture, true);
      c_eval   = P.incremental_eval(station, diff);
    }else{
      c_eval=f_eval;
    }
  }
  return(c_eval);
};

double ApertureILS::improvementIntensity(int beamlet, Station& station, bool open_beamlet, 
                                         double c_eval, Plan & P, bool best_improvement) {
  
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
}

double ApertureILS::openBeamlet(int beamlet, int aperture, Station& station, double c_eval, Plan& P) {
  double aux_eval=0;
  list<pair<int, double> > diff;
  //cout << "open" << beamlet << "a" << aperture<< endl;  
  diff = station.openBeamlet(beamlet, aperture);
  if(diff.size() <1) return(c_eval);
   //cout <<"passed"<< endl;  
  aux_eval = P.incremental_eval(station, diff);
  //cout << "  Opening eval: " << aux_eval << " size: " << aux_diff.size() << " list: ";
  //for (list<pair<int,double>>::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
  //cout << endl;
  return(aux_eval);
}


double ApertureILS::improvementAperture(int beamlet, Station& station, bool open_beamlet, 
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
      aux_eval = closeBeamlet(beamlet, a_list[a], station, c_eval, P); 
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
    local_best = closeBeamlet(beamlet, a_list[local_a], station, c_eval, P);
  } else{
    local_best = openBeamlet(beamlet, a_list[local_a], station, c_eval, P);
  }
  return(local_best);
  
}

bool ApertureILS::acceptanceCriterion(double new_eval, double prev_eval) {
  if (acceptance==0){
    //Only best solutions
    if (new_eval < prev_eval) return(true);
     return(false);
  } else if (acceptance==1) {
    //SA criterion
    double p = exp((double)-(new_eval-prev_eval)/temperature);
    double r = ((double)rand() / (RAND_MAX));
    if (r <= p){
      cout << "  Accept worst. temperature: "<< temperature<<endl;
      return(true);
    }
    return(false);
  }
}

void ApertureILS::updateTemperature(){
  temperature=alpha*temperature;
}

/* Local search procedure used in the HM2019*/
double ApertureILS::localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P) {
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

double ApertureILS::perturbation(Plan& P) {
  list<Station*> stations = P.get_stations();
  int beamlet, aperture;
  list<pair<int,double>> diff;
  double aux_eval=P.getEvaluation();
  
  cout << "##  Perturbation: " ;
  if (perturbation_size==0) cout << "none";   
  else tabu.clear();
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

vector < pair<int, int> > ApertureILS::getShuffledIntensityNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  vector<pair<int, int>> a_list;
  list<Station*>::iterator s;
  int i,j;
  s = stations.begin();
  for (i = 0; i < stations.size(); i++) {
    for (j=0; j< (*s)->getNbApertures();j++){
      //One pair -j (-aperture) for reducing intensity
      //One pair for increasing intensity (+j)
      a_list.push_back(make_pair(i,-(j+1)));
      a_list.push_back(make_pair(i,(j+1)));
    }
    std::advance(s,1);
  }
  std::random_shuffle(a_list.begin(), a_list.end());
  return(a_list);
};

vector < pair<int, int> > ApertureILS::getShuffledIntensityNeighbors(Plan &P, int station){
  list<Station*> stations = P.get_stations();
  vector<pair<int, int>> a_list;
  list<Station*>::iterator s;
  int j;
  s = stations.begin();
  for (j=0; j< (*s)->getNbApertures();j++){
    //One pair -j (-aperture) for reducing intensity
    //One pair for increasing intensity (+j)
    a_list.push_back(make_pair(station,-(j+1)));
    a_list.push_back(make_pair(station,(j+1)));
  }
  std::random_shuffle(a_list.begin(), a_list.end());
  return(a_list);
};

vector < pair<pair<int, int>, pair<int, int>> > ApertureILS::getShuffledApertureNeighbors(Plan &P){
  list<Station*> stations = P.get_stations();
  vector<pair< pair<int,int> , pair<int, int> >> a_list;
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
        if (active.first == -1) continue;
        beamlet = (*st)->getBeamIndex(make_pair(k,pattern.first));
        if((*st)->isOpenBeamlet(beamlet, a))
           a_list.push_back(make_pair(make_pair(s,a) , make_pair(-1,beamlet)));
        else 
            a_list.push_back(make_pair(make_pair(s,a) , make_pair( 1,beamlet)));

        beamlet = (*st)->getBeamIndex(make_pair(k,pattern.second));
        if((*st)->isOpenBeamlet(beamlet, a))
           a_list.push_back(make_pair(make_pair(s,a) , make_pair(-1,beamlet)));
        else 
           a_list.push_back(make_pair(make_pair(s,a) , make_pair( 1,beamlet)));
      }
    }
    std::advance(st,1);
  }
  std::random_shuffle(a_list.begin(), a_list.end());
  return(a_list);
};

// This function performs a local search over all the
// aperture intensities in a treatment plan.
double ApertureILS::iLocalSearch (Plan& P,  double max_time, bool verbose) {
    list<Station*> stations = P.get_stations();
    list<Station*>::iterator s;
    std::clock_t time_end, time_begin;
    double used_time; 
    
    double local_best_eval , aux_eval;
    list<pair<int, double> > diff;
    vector<pair<int, int>> a_list;
    pair <int,int> tabu;
    int i,j, best_n;
    bool improvement = true, best_improvement=ls_type;
    bool completed = false;
    
    //check all neighbors
    tabu = make_pair(-1,-1);
    best_n = -1;
    j=-1;
    local_best_eval = aux_eval = P.getEvaluation();
    if (verbose)
      cout << "Staring intensity local search..." << endl;
    time_begin=clock();
    
    while (improvement) {
      j++;
      a_list = getShuffledIntensityNeighbors(P);
      improvement = false;
      completed = false;
      if (verbose)
        cout << "Neighborhood size "<< a_list.size() << "    current " << local_best_eval << endl;
      for (i = 0; i < a_list.size(); i++) {
         //skip the tabu neighbor (returns the station to previous state)
         if (a_list[i].first == tabu.first && a_list[i].second == tabu.second) continue;
         //get the station of the movement
         s = stations.begin();
         std::advance(s,a_list[i].first);
         //apply step_size intensity change (-(a+1) or +(a+1))
         if (a_list[i].second < 0 ){
           diff = (*s)->modifyIntensityAperture(abs(a_list[i].second)-1, -step_intensity);
           aux_eval = P.incremental_eval(*(*s), diff);
         } else {
           diff = (*s)->modifyIntensityAperture(a_list[i].second-1, step_intensity);
           aux_eval = P.incremental_eval(*(*s), diff);
         }
         if (verbose)
           cout << "  iLS Neighbor " << j << " over station " << (*s)->getAngle() << " aperture " <<  abs(a_list[i].second)-1  << endl;
         // First improvement
         if ((local_best_eval - aux_eval) > 0.00001){
           local_best_eval = aux_eval;
           best_n = i;
           improvement = true;
           if (verbose)
             cout << "     improvement " << aux_eval << endl ;
           
           
           // Add a tabu movement since we found improvement
           tabu.first = a_list[i].first;
           tabu.second = -1*a_list[i].second;
           
           // If first improvement has been chosen break and 
           if (!best_improvement) { 
             if (i==(a_list.size()-1)) completed=true;
             break;
           }
         }
         
         // Undo movement to continue the search
         diff = (*s)->undoLast();
         if (diff.size()>0)
           aux_eval = P.incremental_eval(*(*s), diff);
         
         if (i== (a_list.size()-1)) completed = true;
         
         time_end = clock();
         used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
         if (max_time!=0 && used_time >= max_time) {
           break;
         }
       } 
      
       //Apply best neighbor
       if (improvement && best_improvement) {
         if (a_list[best_n].second < 0 ){
           diff = (*s)->modifyIntensityAperture(abs(a_list[best_n].second)-1, -step_intensity);
           aux_eval = P.incremental_eval(*(*s), diff);
         } else {
           diff = (*s)->modifyIntensityAperture(a_list[best_n].second-1, step_intensity);
           aux_eval = P.incremental_eval(*(*s), diff);
         }
       }
       
       time_end = clock();
       used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
       if (max_time!=0 && used_time >= max_time) {
         break;
       }
    }
    
    cout << "   ils best: " << local_best_eval ;
    if (!completed) cout << ": [nolo]"<< endl;
    else  cout << ": [lo]"<< endl;
    return(local_best_eval);
};

// This function performs a local search over all the
// aperture patterns in a treatment plan.
double ApertureILS::aLocalSearch(Plan& P, double max_time, bool verbose){
  Station *s;
  std::clock_t time_end, time_begin;
  double used_time;
  double local_best_eval, aux_eval;
  bool improvement=true, best_improvement=ls_type;
  vector < pair< pair<int, int>, pair<int, int>> > a_list;
  pair< pair<int,int>, pair<int,int>> best_move;
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
      cout << "Neighborhood size "<< a_list.size() << "    current " << local_best_eval << endl;
    for (i = 0; i < a_list.size(); i++) {
      s = P.get_station(a_list[i].first.first);
      aperture = a_list[i].first.second;
      beamlet = a_list[i].second.second;
      sign = a_list[i].second.first;
      
      if (verbose)
          cout << "  aLS Neighbor " << j << " over station " << s->getAngle() << " aperture " << aperture << " beamlet " << beamlet << endl;
      
      if (sign==-1) {
        aux_eval = closeBeamlet(beamlet, aperture, *s, aux_eval, P); 
      } else {
        aux_eval = openBeamlet(beamlet, aperture, *s, aux_eval, P);
      } 
      
      if ((local_best_eval - aux_eval) > 0.00001){
        local_best_eval = aux_eval;
        best_n = i;
        improvement = true;
        if (verbose)
          cout << "     improvement " << aux_eval << endl ;
        if (!best_improvement) { 
          if (i==(a_list.size()-1)) completed=true;
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
        break;
      }
    }
    
    //Apply best neighbor
    if (improvement && best_improvement) {
      s = P.get_station(a_list[best_n].first.first);
      aperture = a_list[best_n].first.second;
      beamlet = a_list[best_n].second.second;
      sign = a_list[best_n].second.first;
      if (sign == -1) {
        aux_eval = closeBeamlet(beamlet, aperture, *s, aux_eval, P);
      } else{
        aux_eval = openBeamlet(beamlet, aperture, *s, aux_eval, P);
      }
    }
    
    time_end = clock();
    used_time = double(time_end- time_begin) / CLOCKS_PER_SEC;
    if (max_time!=0 && used_time >= max_time) {
      break;
    }
  }
  
  cout << "   als best: " << local_best_eval ;
  if (!completed) cout << ": [nolo]"<< endl;
  else  cout << ": [lo]"<< endl;
  
  return(local_best_eval);
};

// This function performs a local search over all the
// aperture intensities of a station in a treatment plan.
double ApertureILS::iSLocalSearch (Plan& P, bool verbose) {
  pair<bool, pair<Station*, int>> target_beam = getBestLSBeamlet(P);
  list<Station*> stations = P.get_stations();
  list<Station*>::iterator s;
  
  double local_best_eval =P.getEvaluation(), aux_eval=0;
  list<pair<int, double> > diff;
  vector<pair<int, int>> a_list;
  pair <int,int> tabu;
  int i,j, best_n;
  bool improvement = true, best_improvement=ls_type;
  
  // Get station number
  int station=-1;
  s = stations.begin();
  for (i = 0; i < stations.size(); i++) {
    if ((*s)->getAngle() == target_beam.second.first->getAngle()){
      station = i;
      break;
    }
    std::advance(s,1);
  }
  
  //check all neighbors
  tabu = make_pair(-1,-1);
  best_n = -1;
  j=-1;
  while (improvement) {
    j++;
    if (verbose)
      cout << "  iLS Iteration " << j << " over station " << station << endl;
    a_list = getShuffledIntensityNeighbors(P, station);
    improvement = false;
    for (i = 0; i < a_list.size(); i++) {
      //skip the tabu neighbor (returns the station to previous state)
      if (a_list[i].first == tabu.first && a_list[i].second == tabu.second) continue;
      //get the station
      s = stations.begin();
      std::advance(s,a_list[i].first);
      //apply change
      if (a_list[i].second < 0 ){
        diff = (*s)->modifyIntensityAperture(abs(a_list[i].second)-1, -1.0);
        aux_eval = P.incremental_eval(*(*s), diff);
      } else {
        diff = (*s)->modifyIntensityAperture(a_list[i].second-1, 1);
        aux_eval = P.incremental_eval(*(*s), diff);
      }
      //cout << "    Neighbor " << i << " pair <" << a_list[i].first << "," << a_list[i].second << "> " << "eval " << aux_eval;
      // First improvement
      if ((local_best_eval - aux_eval) > 0.00001){
        local_best_eval = aux_eval;
        best_n = i;
        improvement = true;
        if (verbose)
          cout << "     improvement " << aux_eval << endl ;
        if (!best_improvement) { 
          tabu.first = a_list[i].first;
          tabu.second = -1*a_list[i].second;
          break;
        }
      }
      
      diff = (*s)->undoLast();
      if (diff.size()>0)
        aux_eval = P.incremental_eval(*(*s), diff);
    } 
    
    //Apply best neighbor
    if (improvement && best_improvement) {
      tabu.first = a_list[best_n].first;
      tabu.second = -1*a_list[best_n].second;
      if (a_list[best_n].second < 0 ){
        diff = (*s)->modifyIntensityAperture(abs(a_list[best_n].second)-1, -1);
        aux_eval = P.incremental_eval(*(*s), diff);
      } else {
        diff = (*s)->modifyIntensityAperture(a_list[best_n].second-1, 1);
        aux_eval = P.incremental_eval(*(*s), diff);
      }
    }
  }
  return(local_best_eval);
};

// This function performs a local search over all the
// aperture patterns of a station in a treatment plan.
double ApertureILS::aSLocalSearch(Plan& P, bool verbose){
  pair<bool, pair<Station*, int>> target_beam = getBestLSBeamlet(P);
  //If no beamlet was found there is no reason to continue
  if (target_beam.second.second < 0) return(P.getEvaluation());
  
  int beamlet = target_beam.second.second;
  Station* s = target_beam.second.first;
  bool sign = target_beam.first;
  double local_best_eval =P.getEvaluation(), aux_eval=P.getEvaluation();
  bool improvement=true, best_improvement=ls_type;
  vector<int> a_list;
  list<pair<int, double> > diff;
  int i,j, best_n;
  
  best_n = -1;
  j=-1;
  while (improvement) {
    j++;
    if (verbose)
      cout << "  aLS Iteration " << j << " over station " << s->getAngle() << endl;
    //Get neighbor list randomized
    if (sign) a_list = s->getOpen(beamlet);
    else a_list = s->getClosed(beamlet);
    std::random_shuffle(a_list.begin(), a_list.end());    
    
    improvement = false;
    for (i = 0; i < a_list.size(); i++) {
      if (sign) {
        aux_eval = closeBeamlet(beamlet, a_list[i], *s, aux_eval, P); 
      } else {
        aux_eval = openBeamlet(beamlet, a_list[i], *s, aux_eval, P);
      } 
      
      if ((local_best_eval - aux_eval) > 0.00001){
        local_best_eval = aux_eval;
        best_n = i;
        improvement = true;
        if (verbose)
          cout << "     improvement " << aux_eval << endl ;
        if (!best_improvement) { 
          break;
        }
      }
      
      diff = s->undoLast();
      if (diff.size()>0)
        aux_eval = P.incremental_eval(*s, diff);
    }
    
    //Apply best neighbor
    if (improvement && best_improvement) {
      if (sign) {
        aux_eval = closeBeamlet(beamlet, a_list[best_n], *s, aux_eval, P);
      } else{
        aux_eval = openBeamlet(beamlet, a_list[best_n], *s, aux_eval, P);
      }
    }
    
  }
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
  cout << "Local search finished with: " << local_best << " effectively "<< P.getEvaluation()<< endl;
  return(P.getEvaluation());
};

int ApertureILS::getStepIntensity () {
  return(step_intensity);
};

}



