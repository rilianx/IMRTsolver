/*
 * Plan.cpp
 *
 *  Created on: 1 ago. 2018
 *      Author: leslie
 */

#include "Plan.h"


namespace imrt {


  Plan::Plan(Collimator& collimator, vector<Volume>& volumes, int max_apertures,
             int max_intensity, int initial_intensity,
             int step_intensity, StationSetup setup, istringstream* fm_stream, list<int> bac) {

      cout << "##Initilizing plan."<< endl;

      if(bac.size()>0) n_stations=bac.size();
      else n_stations = collimator.getNbAngles();

      for (int i=0;i<n_stations;i++) {
        int angle=collimator.getAngle(i);
        if(bac.size()>0) {angle=bac.front(); bac.pop_front();}

        Station* station = new Station(collimator, volumes, angle,
                                      max_apertures, max_intensity, initial_intensity,
                                      step_intensity, setup, fm_stream);
        add_station(*station);
        angle2station[angle] = station;
      }

  };

  Plan::Plan(const Plan &p)  {

    for (list<Station*>::const_iterator it=p.stations.begin();it!=p.stations.end();it++) {
      Station* aux = new Station(**it);
      //if (p.last_changed && p.last_changed->getAngle()==aux->getAngle()) last_changed=aux;
      add_station(*aux);
      angle2station[aux->getAngle()]=aux;
      //real_stations.push_back(*aux);
    }
    n_stations= stations.size();
  }

  Plan& Plan::operator=(Plan& p){
    for( auto station :stations ) delete station;
    stations.clear();

    for (list<Station*>::const_iterator it=p.stations.begin();it!=p.stations.end();it++) {
      
      Station* aux = new Station (**it);
      add_station(*aux);
      angle2station[aux->getAngle()]=aux;
    }
    n_stations= stations.size();    
    return *this;
  }

  void Plan::newCopy(Plan& p) {
    for( auto station :stations ) delete station;
    stations.clear();
    //real_stations.clear();
    for (list<Station*>::const_iterator it=p.stations.begin();it!=p.stations.end();it++) {
      
      Station* aux = new Station (**it);
      //if (p.last_changed && p.last_changed->getAngle()==aux->getAngle()) 
      //  last_changed=aux;
      add_station(*aux);
      //stations.push_back(aux);
      angle2station[aux->getAngle()]=aux;
      //real_stations.push_back(*aux);
    }
    n_stations= stations.size();
  }

  void Plan::add_station(Station& s){
    stations.push_back(&s);
  }

  double Plan::openBeamlet(int station, int aperture, int beamlet, bool delta_eval,
      list<pair<int, double> >* diff) {
    Station * s = get_station(station);
    double aux_eval=0;
    if(!delta_eval) diff = new list<pair<int, double> >;

    *diff = s->openBeamlet(beamlet, aperture);

    if(!delta_eval){
      if (diff->size() <= 0)
         clearLast();
      delete diff;
      return(0.0);

    }
    return 0.0;
  };

  double Plan::closeBeamlet(int station, int aperture, int beamlet, int side, bool delta_eval,
      list<pair<int, double> >* diff) {


    Station * s = get_station(station);
    if(!delta_eval) diff = new list<pair<int, double> >;
    double l_eval, r_eval;


    if (side == 1) {
      *diff = s->closeBeamlet(beamlet, aperture, true);
    } else if (side==2) {
      *diff = s->closeBeamlet(beamlet, aperture, false);
    } else {
      cout << "ERROR!! you must specify a  valid side for closing a beamlet!" << endl;
      cout << "side " << side << endl;
      getchar();
    }

    if(!delta_eval){
      if (diff->size() <= 0)
        clearLast();
      delete diff;
      return(0.0);
    }
    return 0.0;
    }
    


  

  double Plan::openRow(int station, int aperture, int row, bool side,
		 bool delta_eval, list<pair<int, double> >* diff){
  Station * s = get_station(station);
  double aux_eval=0;
  pair <int, int> pattern = s->getApertureShape(aperture, row);
  pair <int, int> active  = s->collimator.getActiveRange(row, s->getAngle());
  int beamlet;
  double min_new_eval = DBL_MAX;

  if (active.first == -1) return (0.0);

  if (pattern.first == -1) {
    beamlet = active.first + rand()%(active.second-active.first+1);
    //All row are closed find best opening beamlet
    /*for (int i=active.first; i<=active.second; i++) {
      aux_eval= get_delta_eval(*s, i, s->getApertureIntensity(aperture));
      if (aux_eval < min_new_eval){
        beamlet = i;
        min_new_eval = aux_eval;
      };
    }*/
  } else {
    if (side) {
      beamlet = pattern.first - 1;
      if (active.first > beamlet) return(0.0);
    } else {
      beamlet = pattern.second + 1;
      if (active.second < beamlet) return(0.0);
    }
  }

  *diff = s->openBeamlet(beamlet, aperture);

  if(!delta_eval){
    if (diff->size() ==0)
      clearLast();
    delete diff;
    return(0.0);

  }//else
    //return get_delta_eval(*s, *diff);
  return 0.0;
};
  
  list<pair<int, double> > Plan::getOpenRowDiff(int station, int aperture, int row, bool side) {
    list<pair<int, double> > diff;
    Station * s = get_station(station);
    return(s->getOpenRowDiff(row, aperture, side));
  };

  double Plan::closeRow(int station, int aperture, int row, bool side, bool delta_eval, list<pair<int, double> >* diff){
  Station * s = get_station(station);
  pair <int, int> pattern = s->getApertureShape(aperture, row);
  int beamlet;

  if (side)
    beamlet = pattern.first;
  else
    beamlet = pattern.second;

  *diff = s->closeRow(row, aperture, side);
  if(!delta_eval){
    if (diff->size() == 0)
      clearLast();
    delete diff;
    return(0.0);

  }//else
    //return get_delta_eval(*s, *diff);
  return 0.0;
};

  list<pair<int, double> > Plan::getCloseRowDiff(int station, int aperture, int row, bool side) {
    Station * s = get_station(station);
    return(s->getCloseRowDiff(row, aperture, side));
  };
  
  double Plan::modifyIntensityAperture (int station, int aperture, int delta, bool delta_eval, list<pair<int, double> >* diff) {
    Station * s = get_station(station);
    if(!delta_eval) diff = new list<pair<int, double> >;


    *diff = s->modifyIntensityAperture(aperture, delta);
    if(!delta_eval){
      if (diff->size() == 0)
        clearLast();
      
      delete diff;
      return(0.0);
    }//else
      //return get_delta_eval (*s, *diff);
    return 0.0;
  }
  
  list<pair<int, double> > Plan::getModifyIntensityDiff(int station, int aperture, int delta){
    Station * s = get_station(station);
    return(s->getModifyIntensityApertureDiff (aperture, delta));
  }
  
  void Plan::undoLast() {
    if (last_changed==NULL) return;

    list< pair< int, double > > diff = last_changed->undoLast();
    //if (diff.size()>0)
      //incremental_eval (*last_changed, diff);

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

  const list<Station*>& Plan::get_stations() const{
    return stations;
  }

  void Plan::generateApertures(){
	  for(auto s:stations){
		  s->generateApertures();
    }

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
/*
  set < pair< pair<double,bool>, pair<Station*, int> >,
        std::greater < pair< pair<double,bool>, pair<Station*, int> > > >
    Plan::best_beamlets(int n, int nv, int mode) {

    return(ev.best_beamlets(*this, n, nv));
  }
*/
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

  void Plan::printIntensity(int n, bool vector_form) {
    list<Station*>::iterator s= stations.begin();
    advance(s,n);
    (*s)->printIntensity(vector_form);
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
