//
//  Ant.cpp
//  
//
//  Created by Leslie on 07/04/2019.
//

#include "Ant.h"

namespace imrt {
  Ant::Ant (vector<double> w, vector<double> Zmin, vector<double> Zmax, Collimator& _collimator,
       vector<Volume>& volumes, int _max_apertures, int _max_intensity, int initial_intensity,
            int _step_intensity): collimator(_collimator), n_apertures(_max_apertures), n_stations(collimator.getNbAngles()),
            max_intensity(_max_intensity), step_intensity(_step_intensity) {
    
    p = new Plan(w, Zmin, Zmax, collimator, volumes, n_apertures, max_intensity, initial_intensity, step_intensity, -1, Station::OPEN_MIN_SETUP);
    generateReference(); //Remove this if intensityLS is used
  };
  
  Ant::~Ant (){
    delete p;
  };
  
  Ant::Ant(const Ant& a): collimator(a.collimator) {
    p = new Plan (*(a.p));
    ref_plan = a.ref_plan;
    ref_intensity= a.ref_intensity;
    n_stations = a.n_stations;
    n_apertures = a.n_apertures;
    step_intensity= a.step_intensity;
    max_intensity=a.max_intensity;
  };
  
  void Ant::copy(Ant& a) {
    p->newCopy(*a.p);
    for (int s=0; s<n_stations; s++) {
      for (int ap=0; ap<n_apertures; ap++) {
        for (int i=0; i < collimator.getIntensityLevelSize(); i++) {
          ref_intensity[s][ap] = a.ref_intensity[s][ap];
        }
        for (int i=0; i<collimator.getXdim(); i++) {
          ref_plan[s][ap][i] = a.ref_plan[s][ap][i];
        }
      }
    }
  };
  
  void Ant::printIntensity (int i) {
    p->printIntensity(i);
  };
  
  void Ant::printIntensities () {
    for (int s=0; s<n_stations; s++)
       p->printIntensity(s);
  };

  void Ant::printReferencePlan() {
    cout << "--------------------------------------------"<< endl;
    cout << "Reference plan:" << endl;
    cout << "--------------------------------------------"<< endl;
    for (int s=0; s<n_stations; s++) {
      cout << "Station" << s << endl;
      for (int a=0; a<n_apertures; a++) {
        cout << "Aperture " << a << endl;
        for (int i=0; i<collimator.getXdim(); i++) {
          cout << ref_plan[s][a][i] << endl;
          
        }
      }
    }
  };
  
  int Ant::getReferenceAt(int s, int ap, int x) {
    return(ref_plan[s][ap][x]);
  };
  
  int Ant::getIntensityReferenceAt(int s, int ap) {
    return(ref_intensity[s][ap]);
  };
  
  void Ant::generateReference() {
    Station * st;
    for (int s=0; s<n_stations; s++) {
      map <int, vector<int>> aux;
      vector <int> iaux;
      st = p->get_station(s);
      for (int a=0; a<n_apertures; a++) {
        // Set reference to intensity
        iaux.push_back(collimator.getIntensityLevel(st->getApertureIntensity(a)));
        
        // Set reference to aperture pattern
        vector <int> vaux;
        for (int i=0; i<collimator.getXdim(); i++) {
          vaux.push_back(collimator.searchReferenceIndex(st->getApertureShape(a,i)));
        }
        aux[a]=vaux;
      }
      ref_intensity.push_back(iaux);
      ref_plan.push_back(aux);
    }
    //printReferencePlan();
  };
  
  void Ant::generateTour(vector<map <int, Matrix>> & probability, vector<Matrix>& iprobability) {
    double aux;
    bool flag;
    Station * st;
    for (int s = 0; s < n_stations; s++) {
      st = p->get_station(s);
      for (int ap=0; ap < n_apertures; ap++) {
        // Choose intensity
        aux = rand() / double(RAND_MAX);
        flag=false;
        for (int i=0; i< collimator.getIntensityLevelSize() && !flag; i++){
          if (iprobability[s](ap,i) >= aux ) {
            ref_intensity[s][ap]=i;
            st->setApertureIntensity(ap, collimator.getLevelIntensity(i));
            flag=true;
          }
        }
        
        // Choose aperture pattern 
        for (int j = 0; j < collimator.getXdim(); j++) {
          aux = rand() / double(RAND_MAX);
          flag = false;
          for (int r = 0; r < collimator.getReferenceSize() && !flag; r++) {
            if (probability[s][ap](j,r) >= aux) {
              ref_plan[s][ap][j] = r;
              st->setApertureShape(ap,j,collimator.getReference(r));
              flag = true;
            }
          }
        }
      }
      st->generateIntensity();
    }
    p->eval();
  };
  
  double Ant::getEvaluation() {
    return(p->getEvaluation());
  };
  
  Plan * Ant::getPlan() {
    return(p);
  };
  
}
