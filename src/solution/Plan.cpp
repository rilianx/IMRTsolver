/*
 * Plan.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "Plan.h"

namespace imrt {

  Plan::Plan(EvaluationFunction &ev) : ev(ev) {};

  Plan::Plan(EvaluationFunction &ev, vector<double> w, vector<double> Zmin, vector<double> Zmax) : ev(ev), w(w), Zmin(Zmin), Zmax(Zmax) {};
  
  void Plan::add_station(Station& s){
    stations.push_back(&s);
  }
  
  double Plan::eval(vector<double>& w, vector<double>& Zmin, vector<double>& Zmax) {
    double eval=ev.eval(*this,w,Zmin,Zmax);
    ev.generate_voxel_dose_functions ();
    evaluation_fx=eval;
    return eval;
  };
  
  double Plan::eval() {
    double eval=ev.eval(*this,w,Zmin,Zmax);
    ev.generate_voxel_dose_functions ();
    evaluation_fx=eval;
    return eval;
  };
  
  double Plan::incremental_eval (Station& station, list< pair< int, double > >& diff) {
    evaluation_fx=ev.incremental_eval(station, w, Zmin, Zmax, diff);
    last_changed=&station;
    return(evaluation_fx);
  };
  
  double Plan::getEvaluation(){
    return(evaluation_fx);
  }
  
  void Plan::undoLast() {
    list< pair< int, double > > diff = last_changed->undoLast();
    if (diff.size()>0)
      incremental_eval (*last_changed, diff);
  }
  
  const list<Station*>& Plan::get_stations() const{
    return stations;
  }
  
  void Plan::write_open_beamlets(){
    ofstream myfile;
    myfile.open ("openbeamlets.txt");
    
    myfile << "Angles\t";
    for(auto s:stations)
      myfile << s->getAngle() << "\t";
    myfile << endl;
    int k=0;
    for(auto s:stations){
      set<int> open_beamlets;
      myfile << endl << "Station Angle\t" << s->getAngle() << endl;
      for(int i=0; i<s->getNbApertures(); i++){
        myfile << "Aperture\t" << i << endl;
        myfile << "Intensity\t" << s->intensity[i] << endl;
        
        myfile << "OpenBeamlets\t" ;
        bool first=true;
        for(auto beam:s->open_beamlets(i)){
          if(!first) myfile << "\t";
          first=false;
          myfile << k+beam;
        }
        myfile << endl;
      }
      k+=s->getNbBeamlets();
    }
    myfile.close();
    
  }
  
  set < pair< pair<double,bool>, pair<Station*, int> >,
        std::greater < pair< pair<double,bool>, pair<Station*, int> > > >
    Plan::best_beamlets(int n, int nv, int mode) {
    
    return(ev.best_beamlets(*this, n, nv));
  }

}