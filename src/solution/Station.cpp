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

  void Station::setUniformIntensity(double intensity){
    for (int a=0 ; a<max_apertures; a++) {
      for (int i=0; i < collimator.getXdim(); i++) {
        pair <int,int> aux = collimator.getActiveRange(i,angle);
        if (aux.first<0) continue;
        for (int j=aux.first; j<=aux.second; j++) {
          if (j>=A[a][i].first && j<=A[a][i].second){
            change_intensity(i, j, intensity);
          }
        }
      }
    }
  }

  void Station::generateIntensity(){
    pair <int,int>aux;
    clearIntensity();
    //printIntensity();
    for (int a=0 ; a<max_apertures; a++) {
      for (int i=0; i < collimator.getXdim(); i++) {
        aux = collimator.getActiveRange(i,angle);
        if (aux.first<0) continue;
        for (int j=aux.first; j<=aux.second; j++) {
          if (j>=A[a][i].first && j<=A[a][i].second){
            change_intensity(i, j, I(i,j)+intensity[a]);
          }
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

  void Station::change_intensity(int i, int j, double intensity, list< pair< int, double > >* diff ){
    if (intensity==I(i,j)) return;
    if(diff) diff->push_back(make_pair(pos2beam[make_pair(i,j)], intensity-I(i,j)));
    if(I(i,j)>0.0){
       if(int2nb[I(i,j)+0.5]==1) int2nb.erase(I(i,j));
       else int2nb[I(i,j)+0.5]--;
     }

     I(i,j)=intensity;
     if(I(i,j)>0.0) int2nb[I(i,j)+0.5]++;
  }

  list< pair< int, double > > Station::increaseIntensity(int beam, double intensity, int ratio){
	list< pair< int, double > > diff;
    pair<int,int> p = getPos(beam);
    int x=p.first, j=p.second;
    for(int i=max(x-ratio,0); i<min(x+ratio+1,collimator.getXdim()); i++){
      //for(int j=max(y-ratio,0); j<min(y+ratio+1,collimator.getYdim()); j++){
    	  if(I(i,j)==-1) continue;
        if(I(i,j)<-intensity) intensity=-I(i,j);
        change_intensity(i, j, I(i,j)+intensity, &diff);
      //}
    }
    return diff;
  }

  list< pair< int, double > > Station::increaseIntensity_repair(int beam, double intensity, int ratio){
    list< pair< int, double > > diff=increaseIntensity(beam, intensity, ratio);

    pair<int,int> p = getPos(beam);
    int x=p.first, y=p.second;

    if(intensity > 0.0) //reparation (phase 1A)

    for(int i=0; i<collimator.getXdim(); i++){
      int max_j=0; double max_value=-1;
      for(int j=0; j<collimator.getYdim(); j++)
         if(I(i,j)>max_value) {max_j=j; max_value=I(i,j);}

      for(int j=1; j<max_j; j++)
        if(I(i,j) < I(i,j-1))
            change_intensity(i, j, I(i,j-1), &diff);

      for(int j=collimator.getYdim()-2; j>max_j; j--)
        if(I(i,j) < I(i,j+1))
             change_intensity(i, j, I(i,j+1), &diff);
    }

    else //reparation (phase 1B)

    for(int i=0; i<collimator.getXdim(); i++){
      int max_j=0; double max_value=-1;

      if(y<collimator.getYdim()/2){
        for(int j=collimator.getYdim()-1; j>=0; j--)
          if(I(i,j)>max_value) {max_j=j; max_value=I(i,j);}
      }else{
        for(int j=0; j<collimator.getYdim(); j++)
          if(I(i,j)>max_value) {max_j=j; max_value=I(i,j);}
      }

      for(int j=max_j+1; j<collimator.getYdim(); j++)
        if(I(i,j) > I(i,j-1))
             change_intensity(i, j, I(i,j-1), &diff);

      for(int j=max_j-1; j>=0; j--)
        if(I(i,j) > I(i,j+1))
              change_intensity(i, j, I(i,j+1), &diff);
    }

    //reparation (phase 2)
    while(int2nb.size()>max_apertures){
      int from,to, min_diff=10000;
      int prev=0; int nprev=0;
      for(auto i:int2nb  ){
          int n= (prev!=0)? min(i.second,nprev) : i.second;
          if((i.first-prev)*n < min_diff){
            min_diff=(i.first-prev)*n;
            if(nprev<i.second && prev!=0) {from=prev; to=i.first;}
            else {to=prev; from=i.first;}
          }
          prev=i.first;
          nprev=i.second;
      }
    //  printIntensity();
      for(int i=0; i<collimator.getXdim(); i++)
        for(int j=0; j<collimator.getYdim(); j++)
          if(int(I(i,j)+0.5)==from)
            change_intensity(i, j, to, &diff);
    //  printIntensity();

      //cout << int2nb.size() << endl;
      //if( int2nb.size() ==8) exit(0);
    }

    return diff;
  }



}
