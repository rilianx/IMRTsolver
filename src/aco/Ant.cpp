//
//  Ant.cpp
//  
//
//  Created by Leslie on 07/04/2019.
//

#include "Ant.h"

namespace imrt {
  Ant::Ant (vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
       vector<Volume>& volumes, int _max_apertures, int max_intensity, int initial_intensity,
            int step_intensity): collimator(_collimator), n_apertures(_max_apertures), n_stations(collimator.getNbAngles()) {
    
    p = new Plan(w, Zmin, Zmax, collimator, volumes, n_apertures, max_intensity, initial_intensity, step_intensity, -1, Station::OPEN_MIN_SETUP);
    
    Station * st;
    for (int s=0; s<n_stations; s++) {
      map <int, vector<int>> aux;
      st = p->get_station(s);
      for (int a=0; a<n_apertures; a++) {
        vector <int> vaux;
        for (int i=0; i<collimator.getYdim(); i++) {
          vaux.push_back(collimator.searchReferenceIndex(st->getApertureShape(a,i)));
        }
        aux[a]=vaux;
      }
      ref_plan.push_back(aux);
    }
    printReferencePlan();
  };
  
  void Ant::printIntensity (int i) {
    p->printIntensity(i);
  };
  
  void Ant::printReferencePlan() {
    cout << "Reference plan" << endl;
    for (int s=0; s<n_stations; s++) {
      cout << "Station" << s << endl;
      for (int a=0; a<n_apertures; a++) {
        cout << "Aperture " << a << endl;
        for (int i=0; i<collimator.getYdim(); i++) {
          cout << ref_plan[s][a][i] << endl;
          
        }
      }
    }
  };
  
}
