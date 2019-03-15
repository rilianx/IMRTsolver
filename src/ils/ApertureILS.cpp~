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
                         double alpha, bool do_perturbate, int perturbation_size, int acceptance=0, int ls_type=ApertureILS::FIRST_IMPROVEMENT): 
                         ILS(bsize, vsize, acceptance), search_intensity(search_intensity), 
                         search_aperture(search_aperture), prob_intensity(prob_intensity), 
                         step_intensity(step_intensity) , initial_temperature(initial_temperature), 
                         alpha(alpha), do_perturbate(do_perturbate), perturbation_size(perturbation_size), ls_type(ls_type){
//  cout << "per:" << perturbation_size << endl;
  temperature=initial_temperature;
}

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

double ApertureILS::doClose(int beamlet, int aperture, Station& station, double c_eval, Plan& P) {
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
}

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

double ApertureILS::doOpen(int beamlet, int aperture, Station& station, double c_eval, Plan& P) {
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
      aux_eval = doClose(beamlet, a_list[a], station, c_eval, P); 
    } else {
      aux_eval = doOpen(beamlet, a_list[a], station, c_eval, P);
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
    local_best = doClose(beamlet, a_list[local_a], station, c_eval, P);
  } else{
    local_best = doOpen(beamlet, a_list[local_a], station, c_eval, P);
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

}
