/*
 * ApertureILS.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "ApertureILS.h"

namespace imrt {


ApertureILS::ApertureILS(int bsize, int vsize, bool search_intensity, bool search_aperture, double prob_intensity, double initial_temperature, double alpha, int acceptance=0): 
ILS(bsize, vsize, acceptance), search_intensity(search_intensity), search_aperture(search_aperture), prob_intensity(prob_intensity), initial_temperature(initial_temperature), alpha(alpha){
  temperature=initial_temperature;
}

bool ApertureILS::isBeamletModifiable(int beamlet, Station* station, bool open_flag) {
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
  auto it = sb.begin();
  std::advance(it,rand()%sb.size());
  
  Station * s = it->second.first; 
  int beamlet = it->second.second;
  bool sign = it->first.second; 
  
  int i=0;
  while (!isBeamletModifiable(beamlet, s, !sign)) {
    // If no beamlet was found return -1
    i++;
    if (i==sb.size()) return(make_pair(false, make_pair(s,-1)));
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

double ApertureILS::firstImprovementIntensity(int beamlet, Station& station, bool open_beamlet, 
                                              double c_eval, Plan & P) {
  
  double local_best = c_eval, aux_eval=0;
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
      diff = station.modifyIntensityAperture(a_list[a], 1);
      aux_eval = P.incremental_eval(station, diff);
      //cout << endl<<"  Increasing intensity aperture: " << a <<" eval: " << aux_eval << " list: ";
      //for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      //cout << endl;
    } else {
      diff = station.modifyIntensityAperture(a_list[a], -1);
      aux_eval = P.incremental_eval(station, diff);
      //cout << endl<< "  Reducing intensity aperture: " << a <<" eval: " << aux_eval <<  " list: ";
      //for (list<pair<int,double> >::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
      //cout << endl;
    }
    
    // First improvement
    if (local_best> aux_eval){
      local_best=aux_eval;
      local_a=a;
      i=0; // to support best improvement
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
  
  diff = station.openBeamlet(beamlet, aperture);
  if(diff.size() <1) return(c_eval);
  
  aux_eval = P.incremental_eval(station, diff);
  //cout << "  Opening eval: " << aux_eval << " size: " << aux_diff.size() << " list: ";
  //for (list<pair<int,double>>::iterator it=aux_diff.begin();it!=aux_diff.end();it++) cout << it->first << " ";
  //cout << endl;
  return(aux_eval);
}

double ApertureILS::firstImprovementAperture(int beamlet, Station& station, bool open_beamlet, 
                                             double c_eval, Plan& P) {
  double current=c_eval, aux_eval=0, local_best=-1;
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
    if (local_best > aux_eval){
      local_best = aux_eval;
      local_a =a ;
      i=0; // to support best improvement
      return(local_best);
    } 
    
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
  cout << "UPDATE";
  temperature=alpha*temperature;
}

double ApertureILS::localSearch(pair<bool, pair<Station*, int>> target_beam, Plan& P) {
  double aux_eval, local_eval=P.getEvaluation();
  bool search_a;
  
  int beamlet = target_beam.second.second;
  Station* s = target_beam.second.first;
  bool sign = target_beam.first;

  if (search_aperture && !search_intensity) {
    cout << ", ls aperture ";
    aux_eval = firstImprovementAperture(beamlet, *s, !sign, local_eval, P);
    cout << ", found: " << aux_eval;
  } else if (!search_aperture && search_intensity) {
    cout << ", ls intensity ";
    aux_eval = firstImprovementIntensity(beamlet, *s, !sign, local_eval, P);
    cout << ", found: " << aux_eval;
  } else if (search_aperture && search_intensity) {
    if (!sign)
      search_a = s->anyClosed(beamlet);
    else 
      search_a = s->anyOpen(beamlet);

    if (((double) rand() / (RAND_MAX)) <= prob_intensity || !search_a) {
      cout << ", ls intensity";
      aux_eval = firstImprovementIntensity(beamlet, *s,  !sign, local_eval, P);
      cout << ", found: " << aux_eval;
    } else {
      cout << ", ls aperture";
      aux_eval = firstImprovementAperture(beamlet, *s, !sign, local_eval, P);
      cout << ", found: " << aux_eval;
    }
  } else {
    cout << "ERROR: Search option not recognized!!!" << endl;
    
  }
  
  return(aux_eval);
}


}