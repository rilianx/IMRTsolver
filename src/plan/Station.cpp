/*
 * Station.cpp
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */

#include "Station.h"

namespace imrt {


  Station::Station(Collimator& collimator, vector<Volume>& volumes, int _angle, int max_apertures):
		collimator(collimator), angle(_angle) , max_apertures(max_apertures){

    for (int i=0; i<volumes.size(); i++)
      D[i]=&volumes[i].getDepositionMatrix(angle);

    // Initialize empty matrix of intensity
    I = Matrix (collimator.getXdim(), collimator.getYdim());
    for (int i=0; i<collimator.getXdim(); i++) {
      for (int j=0; j<collimator.getYdim(); j++) {
        if (collimator.isActiveBeamAngle(i,j,angle)) {
          I(i,j)=0;
        } else {
          I(i,j)=-1;
        }
      }
    }

    // Iniatialize apertures (alternative representation)
    for (int i=0; i<max_apertures; i++) {
      vector<pair<int,int> > aux;
      for (int j=0; j<collimator.getXdim(); j++) {
        aux.push_back(collimator.getActiveRange(j,angle));
      }
      A.push_back(aux);
      intensity.push_back(1);
    }


  };

  void Station::clearIntensity() {
    pair<int,int> aux;
    for (int i=0; i < collimator.getXdim(); i++) {
      aux = collimator.getActiveRange(i,angle);
      if (aux.first<0) continue;
      for (int j=aux.first; j<=aux.second; j++) I(i,j)=0;
    }
  };

  void Station::generateIntensity(){
    pair <int,int>aux;
    clearIntensity();
    //printIntensity();
    for (int a=0 ; a<max_apertures; a++) {
      for (int i=0; i < collimator.getXdim(); i++) {
        aux = collimator.getActiveRange(i,angle);
        if (aux.first<0) continue;
        for (int j=aux.first; j<=aux.second; j++) {
          if (j>=A[a][i].first && j<=A[a][i].second)
            I(i,j)+=intensity[a];
        }
      }
    }
  };

  // Function to be used to get the position
  // in the matrix I of a beam column of matrix D
  pair<int,int> Station::getPos(int index) {
    return(collimator.indexToPos(index, angle));
  }

  void Station::printIntensity() {
    cout << "Angle "<< angle <<" intensity matrix:"<< endl;
    for (int i=0; i<collimator.getXdim();i++) {
      for (int j=0; j<collimator.getYdim(); j++) {
        cout << I(i,j)<< ",";
      }
      cout << endl;
    }
  }

  void Station::printApertures() {
    cout << "Angle "<< angle << endl;
    for (int a=0; a<max_apertures;a++) {
      cout << "aperture " << a << endl;
      for (int i=0; i<collimator.getXdim(); i++) {
        cout << A[a][i].first << "-" << A[a][i].second << endl;
      }
    }
  }

  Matrix& Station::getDepositionMatrix(int o){
    return *D[o];
  }

}
