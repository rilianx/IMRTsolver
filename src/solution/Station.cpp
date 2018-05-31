/*
 * Station.cpp
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */

#include "Station.h"

namespace imrt {


  Station::Station(Collimator& collimator, vector<Volume>& volumes, int _angle, int max_apertures):
		collimator(collimator), angle(_angle) , max_apertures(max_apertures), A(max_apertures), intensity(max_apertures){

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
      A[i]=aux;
      intensity[i]=1;
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
  pair<int,int> Station::getPos(int index) const{
    auto pos = beam2pos.find(index);
    if(pos==beam2pos.end()){
      beam2pos[index]=collimator.indexToPos(index, angle);
      pos2beam[beam2pos[index]]=index;
      return beam2pos[index];
    }
    return(pos->second);
  }

  void Station::printIntensity() {
    cout << "Angle "<< angle <<" intensity matrix:"<< endl;
    for (int i=0; i<collimator.getXdim();i++) {
      for (int j=0; j<collimator.getYdim(); j++) {
        printf("%2.0f ", I(i,j));
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

  const Matrix& Station::getDepositionMatrix(int o) const{
    return *(D.find(o)->second);
  }

  list< pair< int, double > > Station::increaseIntensity(int beam, double intensity, int ratio){

	list< pair< int, double > > diff;
    pair<int,int> p = getPos(beam);
    int x=p.first, y=p.second;
    for(int i=max(x-ratio,0); i<min(x+ratio+1,collimator.getXdim()-1); i++){
      for(int j=max(y-ratio,0); j<min(y+ratio+1,collimator.getYdim()-1); j++){
        if(I(i,j)<-intensity) intensity=-I(i,j);
        I(i,j)+=intensity;
        diff.push_back(make_pair(pos2beam[make_pair(i,j)], intensity));
      }
    }
    return diff;
  }



}
