/*
 * Station.cpp
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */

#include "Station.h"

namespace imrt {


  Station::Station(Collimator& collimator, vector<Volume>& volumes, int _angle, int max_apertures, 
                   int max_intensity, int initial_intensity, int step_intensity, int open_apertures, int setup):
		collimator(collimator), angle(_angle) , max_apertures(max_apertures), A(max_apertures), 
		intensity(max_apertures), max_intensity(max_intensity), initial_intensity(initial_intensity),
                step_intensity(step_intensity) {
    min_intensity=0;
    if(open_apertures==-1) open_apertures=max_apertures;

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
    initializeStation(setup, open_apertures); 
    /*for (int i=0; i<max_apertures; i++) {
      vector<pair<int,int> > aux;
      for (int j=0; j<collimator.getXdim(); j++) {
        if (open_apertures>0)
          aux.push_back(collimator.getActiveRange(j,angle));
        else
          aux.push_back(make_pair(-1,-1));
      }
      A[i]=aux;
      intensity[i]=initial_intensity;
      open_apertures--;
    }*/
    last_mem= make_pair(make_pair(-1,-1), make_pair(-1,-1));
  }

  Station::Station(const Station &s): collimator(s.collimator){
    angle=s.angle;
    max_apertures=s.max_apertures;
    max_intensity=s.max_intensity;
    min_intensity=s.min_intensity;
    initial_intensity=s.initial_intensity;
    D=s.D;
    I=s.I;
    intensity=s.intensity;
    A=s.A;
    last_mem=s.last_mem;
    last_diff=s.last_diff;
    pos2beam=s.pos2beam;
    beam2pos=s.beam2pos;
    int2nb=s.int2nb;
  }

  Station& Station::operator=(const Station & s) {
    collimator=s.collimator;
    angle=s.angle;
    max_apertures=s.max_apertures;
    max_intensity=s.max_intensity;
    min_intensity=s.min_intensity;
    initial_intensity=s.initial_intensity;
    D=s.D;
    I=s.I;
    intensity=s.intensity;
    A=s.A;
    last_mem=s.last_mem;
    last_diff=s.last_diff;
    pos2beam=s.pos2beam;
    beam2pos=s.beam2pos;
    int2nb=s.int2nb;

    return(*this);
  }

  void Station::initializeStation(int type, int open_apertures) {
    // Generating aperture patterns
    if (type==OPEN_MAX_SETUP || type==OPEN_MIN_SETUP) {
      for (int i=0; i<max_apertures; i++) {
        vector<pair<int,int> > aux;
        for (int j=0; j<collimator.getXdim(); j++)
          aux.push_back(collimator.getActiveRange(j,angle));
	      A[i] =aux;
      }
    } else if (type==CLOSED_MAX_SETUP || type==CLOSED_MIN_SETUP) {
      for (int i=0; i<max_apertures; i++){
        vector<pair<int,int> > aux;
        for (int j=0; j<collimator.getXdim(); j++)
          aux.push_back(make_pair(-1,-1));
        A[i] = aux;
      }
    } else if (type==RAND_RAND_SETUP) {
      for (int i=0; i<max_apertures; i++) {
        vector<pair<int,int> > aux;
        for (int j=0; j<collimator.getXdim(); j++) {
          pair<int, int> range = collimator.getActiveRange(j,angle);
          if (range.first<0) { 
            aux.push_back(make_pair(-1,-1));
            continue;
          }
          int index1 = range.first + rand() % (range.second-range.first+1);
          if (index1 == range.second) {
            aux.push_back(make_pair(range.second,range.second));
          } else {
            int index2 = index1 + rand() % (range.second-index1+1);
            aux.push_back(make_pair(index1,index2));
          }
        }
        A[i]=aux;
      }
    } else {
      //Keeping this for backwards compatibility
      for (int i=0; i<max_apertures; i++) {
        vector<pair<int,int> > aux;
        for (int j=0; j<collimator.getXdim(); j++) {
          if (open_apertures>0)
            aux.push_back(collimator.getActiveRange(j,angle));
          else
            aux.push_back(make_pair(-1,-1));
        }
        open_apertures--; 
        A[i]=aux;
      }
    }

    // Generating intensity
    if (type==OPEN_MIN_SETUP || type==CLOSED_MIN_SETUP) {
      for (int i=0; i<max_apertures; i++)
        intensity[i] = min_intensity;
    } else if (type==OPEN_MAX_SETUP || type==CLOSED_MAX_SETUP) {
      for (int i=0; i<max_apertures; i++)
        intensity[i] = max_intensity;
    } else if (type==RAND_RAND_SETUP) {
      int n_levels = ((max_intensity-min_intensity) / step_intensity);
      for (int i=0; i<max_apertures; i++)
        intensity[i] = min_intensity + step_intensity * (rand() %  (n_levels+1));
    } else {
      for (int i=0; i<max_apertures; i++)
        intensity[i] = initial_intensity;
    }
    generateIntensity();
  }

  void Station::clearIntensity() {
    pair<int,int> aux;
    for (int i=0; i < collimator.getXdim(); i++) {
      aux = collimator.getActiveRange(i,angle);
      if (aux.first<0) continue;
      for (int j=aux.first; j<=aux.second; j++) I(i,j)=0;
    }
  }

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
        if (aux.first<0 || A[a][i].first<0) continue;
        //for (int j=aux.first; j<=aux.second; j++) {
        for (int j=A[a][i].first; j<=A[a][i].second; j++) {
          //if (j>=A[a][i].first && j<=A[a][i].second){
            //cout << "a:" <<a << " i:" << i << " j:" <<j << " to add:" <<  intensity[a] << " on:" << I(i,j) <<  endl;
            change_intensity(i, j, I(i,j)+intensity[a]);
          //}
        }
      }
    }
  }

  list<int> Station::open_beamlets(int a){
	  list<int> ob;
	  pair <int,int>aux;

    for (int i=0; i < collimator.getXdim(); i++) {
      aux = collimator.getActiveRange(i,angle);
      if (aux.first<0) continue;
      for (int j=aux.first; j<=aux.second; j++)
        if (j>=A[a][i].first && j<=A[a][i].second)
        	ob.push_back(pos2beam[make_pair(i,j)]);

    }
    return ob;
  }

  list<int> Station::closed_beamlets(int a){
	  list<int> ob;
	  pair <int,int>aux;
    for (int i=0; i < collimator.getXdim(); i++) {
      aux = collimator.getActiveRange(i,angle);
      if (aux.first<0) continue;
      for (int j=aux.first; j<=aux.second; j++)
        if (j<A[a][i].first || j>A[a][i].second)
        	 ob.push_back(pos2beam[make_pair(i,j)]);
    }
    return ob;
  }

  /* Function to be used to get the position
  in the matrix I of a beam column of matrix D*/
  pair<int,int> Station::getPos(int index) const{
    auto pos = beam2pos.find(index);
    if(pos==beam2pos.end()){
      beam2pos[index]=collimator.indexToPos(index, angle);
      pos2beam[beam2pos[index]]=index;
      return beam2pos[index];
    }
    return(pos->second);
  }

  void Station::printIntensity(bool vector_form) {
	  if(!vector_form){
		  cout << "Angle "<< angle <<" intensity matrix:"<< endl;
		  for (int i=0; i<collimator.getXdim();i++) {
			  for (int j=0; j<collimator.getYdim(); j++) {
				  printf("%2.0f ", I(i,j));
			  }
			  cout << endl;
		  }
	  }else{
		  for (int i=0; i<collimator.getXdim();i++) {
			  for (int j=0; j<collimator.getYdim(); j++) {
				  if(I(i,j)!=-1.0) printf("\t%2.1f", I(i,j));
			  }
		  }
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

  void Station::printAperture(int aperture) {
    cout << "Angle: "<< angle << " aperture: " << aperture << endl;
    for (int i=0; i<collimator.getXdim(); i++) {
      cout << A[aperture][i].first << "-" << A[aperture][i].second << endl;
    }
  }

  const Matrix& Station::getDepositionMatrix(int o) const{
    const Matrix* aux=D.find(o)->second;
    return (*aux);
  }

  void Station::change_intensity(int i, int j, double intensity, list< pair< int, double > >* diff ) {
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
        if((int) I(i,j)+ (int) intensity > (int) max_intensity) intensity=max_intensity-I(i,j);
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

  bool Station::isOpenBeamlet(int beam, int aperture) {
  	auto coord = collimator.indexToPos(beam, angle);

    //cout << "Checking beamlet: " << beam << " angle:" << aperture << endl;
    //cout << "Position: " << coord.first << "," << coord.second << endl;
  	if (A[aperture][coord.first].first < 0 || A[aperture][coord.first].first > coord.second || A[aperture][coord.first].second< coord.second)
  	  return(false);
  	//cout << "Not closed" << endl;
  	return(true);
  }

  vector<int> Station::getClosed(int beam) {
    vector<int> closed;
    for (int i =0;i<max_apertures;i++){
      if (!isOpenBeamlet(beam,i)) closed.push_back(i);
    }
    return(closed);
  }

  vector<int> Station::getOpen(int beam){
    vector<int> open;
    for (int i =0;i<max_apertures;i++){
      if (isOpenBeamlet(beam,i)) open.push_back(i);
    }
    return(open);
  }

  bool Station::anyOpen(int beam){
    for (int i =0;i<max_apertures;i++)
      if (isOpenBeamlet(beam,i)) return(true);
    return(false);
  }

  bool Station::anyClosed(int beam){
    for (int i =0;i<max_apertures;i++)
      if (!isOpenBeamlet(beam,i)) return(true);
    return(false);
  }

  bool Station::isActiveBeamlet(int beam) {
    auto coord = collimator.indexToPos(beam, angle);
    if(!collimator.isActiveBeamAngle(coord.first, coord.second,angle)) {
    	//cout << "Not active" << endl;
    	return(false);
    }
    return(true);
  }

  /* Function that closes a beamlet from the left, if lside is true, or
     from the right size otherwise. Return true if the closing was performed.*/
  list<pair <int,double> > Station::closeBeamlet(int beam, int aperture, bool lside) {
    //cout << "Attempt to close beam: " << beam << endl;
    list<pair <int, double> > diff;
    //cout << "active: " << isActiveBeamlet(beam) << "open: " << isOpenBeamlet(beam, aperture) << endl;

    if (isActiveBeamlet(beam) && isOpenBeamlet(beam, aperture)) {
      auto coord = collimator.indexToPos(beam, angle);
      int row= coord.first;
      //cout << "Coordinates: " << coord.first << "," << coord.second << endl;
      if (lside) {
        for (int i=0;i<=coord.second-A[aperture][row].first;i++) {
          diff.push_back(make_pair(beam-(coord.second-A[aperture][row].first)+i, -intensity[aperture]));
        }
        last_mem = make_pair(make_pair(aperture,row), A[aperture][row]);
        if (A[aperture][row].second == coord.second) {
          A[aperture][row].first=-1;
          A[aperture][row].second=-1;
        } else {
          A[aperture][row].first=coord.second+1;
        }
      } else {
        for (int i=0;i<=A[aperture][row].second-coord.second;i++) {
          diff.push_back(make_pair(beam+i, -intensity[aperture]));
	      }
        last_mem = make_pair(make_pair(aperture,row), A[aperture][row]);
	      if (A[aperture][row].first == coord.second) {
	        A[aperture][row].first=-1;
	        A[aperture][row].second=-1;
	      } else {
          A[aperture][row].second=coord.second-1;
        }
      }
      updateIntensity(diff);
    }

    last_diff=diff;
    return(diff);
  }

  /* Function that opens a beamlet from the left, if lside is true, or
     from the right size otherwise. Return true if the closing was performed.*/
  list <pair<int,double> > Station::openBeamlet(int beam, int aperture) {
   // cout << "Attempt to open beam: " << beam << endl;
    list<pair<int, double>> diff;
    //cout << "active: " << isActiveBeamlet(beam) << " open: " << isOpenBeamlet(beam, aperture) << " aperture: " << aperture << endl;
    if (isActiveBeamlet(beam) && !isOpenBeamlet(beam, aperture)) {
      auto coord = collimator.indexToPos(beam, angle);
      int row= coord.first;
      //cout << "Coordinates: " << coord.first << "," << coord.second << endl;
      last_mem = make_pair(make_pair(aperture,row), A[aperture][row]);
      if (A[aperture][row].first < 0) {
        //When the row is completely closed
        diff.push_back(make_pair(beam, intensity[aperture]));
        A[aperture][row].first=coord.second;
        A[aperture][row].second=coord.second;
      } else {
        if (A[aperture][row].first > coord.second) {
          for (int i=0;i<(A[aperture][row].first-coord.second);i++)
            diff.push_back(make_pair(beam+i, intensity[aperture]));
          A[aperture][row].first = coord.second;
        } else {
          for (int i=0;i<(coord.second-A[aperture][row].second);i++)
            diff.push_back(make_pair(beam-(coord.second-A[aperture][row].second)+1+i, intensity[aperture]));
          A[aperture][row].second = coord.second;
        }
      }
		 	updateIntensity(diff);
		 //	printIntensity(true);
		}/* else {
		  cout << "Not active or not closed" << endl;
		} */
		last_diff=diff;
		return(diff);
  }

  list <pair< int,double> > Station::modifyIntensityAperture(int aperture, double size) {
    list < pair <int, double > > diff;
    if (intensity[aperture]+size < 0 || intensity[aperture]+size>max_intensity) {
      if (intensity[aperture]+size < 0  && intensity[aperture]!=0) size = intensity[aperture];
      else if (intensity[aperture]+size>max_intensity && intensity[aperture]!=max_intensity) size = max_intensity-intensity[aperture];
      else return (diff);
    }

    intensity[aperture]=intensity[aperture]+size;

    for (int i=0; i<collimator.getXdim(); i++) {
      if (A[aperture][i].first<0) continue;
      int beamlet = pos2beam[make_pair(i, A[aperture][i].first)];
      for(int j=A[aperture][i].first; j<=A[aperture][i].second; j++){
        diff.push_back(make_pair(beamlet, size));
        beamlet++;
      }
    }
    updateIntensity(diff);
    last_diff=diff;
    return (diff);
  }

  void Station::updateIntensity(list<pair<int,double> > diff) {
    pair<int,int> coord;
    for (list<pair<int,double> >::iterator it=diff.begin();it!=diff.end();it++) {
      coord=collimator.indexToPos(it->first, angle);
     // cout << "THIS: " << coord.first << " " << coord.second << endl;
      I(coord.first,coord.second) = I(coord.first,coord.second) + it->second;
    }
  }

  bool Station::canIncreaseIntensity(int beam){
    for (int i =0;i<max_apertures;i++)
      if (isOpenBeamlet(beam,i) && intensity[i]<max_intensity) return(true);
      return(false);
  }

  bool Station::canReduceIntensity(int beam){
    for (int i =0;i<max_apertures;i++)
      if (isOpenBeamlet(beam,i) && intensity[i]>0) return(true);
    return(false);
  }

  list <pair<int,double> > Station::undoLast () {
    list <pair<int,double> > undo_diff;
    //if (last_mem.first.first<0) {
      //return(undo_diff);
    //}
    pair<int,int> coord;

    if (last_mem.first.first>=0)
      A[last_mem.first.first][last_mem.first.second] = last_mem.second;

    for (list<pair<int,double> >::iterator it=last_diff.begin();it!=last_diff.end();it++) {
      coord=collimator.indexToPos(it->first, angle);
      I(coord.first,coord.second) = I(coord.first,coord.second) - (it->second);
      undo_diff.push_back(make_pair(it->first, -(it->second)));
    }

    last_mem= make_pair(make_pair(-1,-1), make_pair(-1,-1));
    last_diff.clear();

    return(undo_diff);
  }

  void Station::clearHistory() {
    last_mem= make_pair(make_pair(-1,-1), make_pair(-1,-1));
    last_diff.clear();
  }
}
