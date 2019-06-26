/*
 * Plan.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "Plan.h"

namespace imrt {

  Plan::Plan (EvaluationFunction &ev) : ev(ev), last_changed(NULL) {};

  Plan::Plan (EvaluationFunction &ev, vector<double> w, vector<double> Zmin, vector<double> Zmax):
             ev(ev), w(w), Zmin(Zmin), Zmax(Zmax) {
    last_changed=NULL;
  }

  Plan::Plan(vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& collimator,
             vector<Volume>& volumes, int max_apertures, int max_intensity, int initial_intensity,
             int step_intensity, int setup, char* file) :
             ev(EvaluationFunction::getInstance(volumes, collimator)), w(w), Zmin(Zmin), Zmax(Zmax) {

    cout << "##Initilizing plan."<< endl;

    fstream* myfile=NULL;
    if(file) myfile=new fstream(file, std::ios_base::in);

    for (int i=0;i<collimator.getNbAngles();i++) {
      Station* station = new Station(collimator, volumes, collimator.getAngle(i), max_apertures,
                                     max_intensity, initial_intensity, step_intensity,
                                     setup, myfile);
      add_station(*station);
      angle2station[collimator.getAngle(i)] = station;
    }
    n_stations= stations.size();
    if(myfile){
    	myfile->close();
        delete myfile;
    }
    cout << "##  Created " << stations.size() << " stations."<< endl;
    eval();
    cout << "##  Initial evaluation: " << evaluation_fx << "."<< endl;
    last_changed=NULL;

    //ev.create_beam2voxel_list(*this);
  };

  Plan::Plan(const Plan &p): ev(p.ev), w(p.w), Zmin(p.Zmin), Zmax(p.Zmax), last_changed(NULL) {
    //EvaluationFunction aux_ev(p.ev);
    //ev=aux_ev;
    for (list<Station*>::const_iterator it=p.stations.begin();it!=p.stations.end();it++) {
      Station* aux = new Station(**it);
      if (p.last_changed && p.last_changed->getAngle()==aux->getAngle()) last_changed=aux;
      add_station(*aux);
      angle2station[aux->getAngle()]=aux;
      //real_stations.push_back(*aux);
    }
    evaluation_fx=p.evaluation_fx;
  }

  void Plan::newCopy(Plan& p) {
    last_changed=NULL;
    //ev= EvaluationFunction::getInstance();

    w=p.w;
    Zmin=p.Zmin;
    Zmax=p.Zmax;
    stations.clear();
    //real_stations.clear();
    for (list<Station*>::const_iterator it=p.stations.begin();it!=p.stations.end();it++) {
      Station* aux = new Station (**it);
      if (p.last_changed!=NULL && p.last_changed->getAngle()==aux->getAngle()) last_changed=aux;
      stations.push_back(aux);
      angle2station[aux->getAngle()]=aux;
      //real_stations.push_back(*aux);
    }
    evaluation_fx=p.evaluation_fx;
  }

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
    evaluation_fx = ev.incremental_eval(station, w, Zmin, Zmax, diff);
    last_changed = &station;
    last_diff = diff;
    return(evaluation_fx);
  };

  double Plan::getEvaluation(){
    return(evaluation_fx);
  }

  double Plan::openBeamlet(int station, int aperture, int beamlet) {
    Station * s = get_station(station);
    double aux_eval=0;
    list<pair<int, double> > diff;

    if (get_delta_eval(*s, beamlet, s->getApertureIntensity(aperture)) > evaluation_fx){
			clearLast();
      return(evaluation_fx);
    }

    diff = s->openBeamlet(beamlet, aperture);
    if(diff.size() <1) {
      clearLast();
      return(evaluation_fx);
    }

    return(incremental_eval(*s, diff));
  };

  double Plan::closeBeamlet(int station, int aperture, int beamlet, int side) {
    Station * s = get_station(station);
    list<pair<int, double> > diff;
    double l_eval, r_eval;

    if (get_delta_eval(*s, beamlet, -(s)->getApertureIntensity(aperture)) > evaluation_fx){
      clearLast();
      return(evaluation_fx);
    }

    if (side == 1) {
      diff = s->closeBeamlet(beamlet, aperture, true);
      if (diff.size() > 0) {
        incremental_eval(*s, diff);
      } else {
        clearLast();
      }
    } else if (side==2) {
      diff = s->closeBeamlet(beamlet, aperture, false);
      if (diff.size() > 0) {
        incremental_eval(*s, diff);
      } else {
        clearLast();
      }
    } else {
      cout << "ERROR!! you must specify a  valid side for closing a beamlet!" << endl;
      cout << "side " << side << endl;
      getchar();
    }
    return(evaluation_fx);
  };

  double Plan::modifyIntensityAperture (int station, int aperture, int delta) {
    Station * s = get_station(station);
    list<pair<int, double> > diff;

    diff = s->modifyIntensityAperture(aperture, delta);
    if (diff.size() > 0) {
      incremental_eval(*s, diff);
    } else {
      clearLast();
    }
    return(evaluation_fx);
  }

  void Plan::undoLast() {
    if (last_changed==NULL) return;

    list< pair< int, double > > diff = last_changed->undoLast();
    if (diff.size()>0)
      incremental_eval (*last_changed, diff);

    last_changed->clearHistory();
    last_changed = NULL;
    last_diff.clear();
  };

  void Plan::clearLast() {
		if (last_changed==NULL) return;
    last_changed->clearHistory();
    last_changed = NULL;
    last_diff.clear();
  };

  //Lepi's version
  void Plan::undoLast2() {
    last_changed->revert(last_diff);
    ev.undo_last_eval(w,Zmin,Zmax);
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

  void Plan::printSolution () {
    Station *st;
    for (int s=0; s < n_stations; s++) {
      st = get_station(s);
      st->printAperture(s);
    }
  };
  
  string Plan::toStringApertures () {
    Station *st;
    string solution_str="";
    for (int s=0; s < stations.size(); s++) {
      st = get_station(s);
      solution_str = solution_str + st->toStringApertures();
    }
    return(solution_str);
  };
  
  string Plan::toStringIntensities () {
    Station *st;
    string solution_str="";
    for (int s=0; s < stations.size(); s++) {
      st = get_station(s);
      solution_str = solution_str + st->toStringIntensities();
    }
    return(solution_str);
  };

  void Plan::printIntensity(int n) {
    list<Station*>::iterator s= stations.begin();
    advance(s,n);
    (*s)->printIntensity(false);
  }

  void Plan::writeIntensities(string file, int n) {
  	ifstream myfile;
  	myfile.open(file);

    std::string line;
    for(int i=0; i<n; i++) std::getline(myfile, line);

  	int angle1, angle2, angle3, angle4, angle5;
  	double F;
  	myfile >> angle1;
  	myfile >> angle2;
  	myfile >> angle3;
  	myfile >> angle4;
  	myfile >> angle5;
  	myfile >> F;

  	for(auto station:stations){
  		station->writeIntensity(myfile);
  	}

  }

}
