/*
 * Station.cpp
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */

#include "Station.h"

namespace imrt {


  Station::Station(Collimator& collimator, vector<Volume>& volumes,
                   int _angle, int max_apertures, int max_intensity,
                   int initial_intensity, int step_intensity,
                   StationSetup setup, fstream* myfile): collimator(collimator),
                   angle(_angle) , max_apertures(max_apertures),
                   A(max_apertures), intensity(max_apertures),
                   max_intensity(max_intensity),
                   initial_intensity(initial_intensity),
                   step_intensity(step_intensity) {
    min_intensity=1;

    n_volumes=volumes.size();
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
    initializeStation(setup);

    // LESLIE: Is this compatible with the aperture-based
    // representation?? Don't think so, no?
    // TODO: define intensity vector to adjust the following
    //       initialization to aperture-based mode.
    // TODO: this should be in the initialize Station
    //       function
    // File intensity initialization
    if(myfile) {
    	for (int i=0; i<collimator.getXdim(); i++) {
    	   for (int j=0; j<collimator.getYdim(); j++) {
    		   int a;*myfile>>a;
    		   if(a!=-1) change_intensity(i, j, a);
    	   }
    	}
    }

    last_mem = make_pair(make_pair(-1,-1), make_pair(-1,-1));
    last_intensity = intensity;
  };

  Station::Station(const Station &s): collimator(s.collimator){
    collimator=s.collimator;
    angle=s.angle;
    max_apertures=s.max_apertures;
    max_intensity=s.max_intensity;
    min_intensity=s.min_intensity;
    initial_intensity=s.initial_intensity;
    step_intensity=s.step_intensity;
    n_volumes=s.n_volumes;
    for (int i=0; i<n_volumes; i++){
      const Matrix * aux= s.D.find(i)->second;
      D[i]=aux;
    }
    //I=s.I;
    I= Matrix (s.I);
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
    step_intensity=s.step_intensity;
    n_volumes=s.n_volumes;
    for (int i=0; i<n_volumes; i++) {
      const Matrix * aux= s.D.find(i)->second;
      D[i]=aux;
    }
    I= Matrix (s.I);
    //I= s.I;
    intensity=s.intensity;
    A=s.A;
    last_mem=s.last_mem;
    last_diff=s.last_diff;
    pos2beam=s.pos2beam;
    beam2pos=s.beam2pos;
    int2nb=s.int2nb;

    return(*this);
  }



  void Station::initializeStation(StationSetup type) {
    // Generating aperture patterns
    if (type==StationSetup::open_all_max || type==StationSetup::open_all_min) {
      for (int i=0; i<max_apertures; i++) {
        vector<pair<int,int> > aux;
        for (int j=0; j<collimator.getXdim(); j++){
          aux.push_back(collimator.getActiveRange(j,angle));
        }
	      A[i] =aux;
      }
    } else if (type==StationSetup::closed_all_max ||
               type==StationSetup::closed_all_min) {
      for (int i=0; i<max_apertures; i++){
        vector<pair<int,int> > aux;
        for (int j=0; j<collimator.getXdim(); j++)
          aux.push_back(make_pair(-1,-1));
        A[i] = aux;
      }
    } else if (type==StationSetup::rand_all_rand) {
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
    } else if (type==StationSetup::open_min_min
               || type==StationSetup::open_min_k) {
      // First aperture is open
      vector<pair<int,int> > aux;
      for (int j=0; j<collimator.getXdim(); j++){
        aux.push_back(collimator.getActiveRange(j,angle));
      }
	    A[0] =aux;
      // All the other apertures are closed
      for (int i=1; i<max_apertures; i++){
        vector<pair<int,int> > aux;
        for (int j=0; j<collimator.getXdim(); j++)
          aux.push_back(make_pair(-1,-1));
        A[i] = aux;
      }
    }

    // Generating intensity
    if (type==StationSetup::open_all_min || type==StationSetup::closed_all_min
        || type==StationSetup::open_min_min) {
      for (int i=0; i<max_apertures; i++)
        intensity[i] = min_intensity;
    } else if (type==StationSetup::open_all_max || type==StationSetup::closed_all_max) {
      for (int i=0; i<max_apertures; i++)
        intensity[i] = max_intensity;
    } else if (type==StationSetup::rand_all_rand) {
      vector<int> levels;
      for (int k=0; k<=max_intensity; k=k+step_intensity)
        levels.push_back(k);
      int sel;
      for (int i=0; i<max_apertures; i++) {
        sel = (rand() %  levels.size());
        intensity[i] = levels[sel];
      }
    } else if (type==StationSetup::open_min_k) {
      intensity[0] = initial_intensity;
      for (int i=1; i<max_apertures; i++)
        intensity[i] = min_intensity;
    } else {
      for (int i=0; i<max_apertures; i++)
        intensity[i] = initial_intensity;
    }

    if(type==StationSetup::rand_int) {
      generate_random_intensities();
    } else {
      generateIntensityMatrix();
    }

  }

  //only for ILS
  void Station::generate_random_intensities(){
    clearIntensityMatrix();

    vector<int> values(max_apertures+1);
    values[0] = 0;
    for(int i=1; i<max_apertures+1;i++){
      int in=rand()%(max_intensity+1)/2.0;
      values[i]=in;
    }

    for (int i=0; i < collimator.getXdim(); i++) {
      list<int> in;
      pair <int,int> aux = collimator.getActiveRange(i,angle);
      if(aux.first<0) continue;

      for (int j=aux.first; j<=aux.second; j++)
        in.push_back(values[rand()%(max_apertures+1)]);
      in.sort();

      int left=aux.first,right=aux.second;
      while(!in.empty()){
        if(rand()%2)
          change_intensity(i, left++, in.front());
        else
          change_intensity(i, right--, in.front());
        in.pop_front();
      }
    }
  //  printIntensity();
  }

  void Station::clearIntensityMatrix() {
    pair<int,int> aux;
    for (int i=0; i < collimator.getXdim(); i++) {
      aux = collimator.getActiveRange(i,angle);
      if (aux.first<0) continue;
      for (int j=aux.first; j<=aux.second; j++)
        change_intensity(i, j, 0);
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
  };

  void Station::change_intensity(int i, int j, double intensity, list< pair< int, double > >* diff ) {
    if (intensity == I(i,j)) return;
    if(diff)
      diff->push_back(make_pair(pos2beam[make_pair(i,j)], intensity-I(i,j)));
    if(I(i,j)>0.0){
       if(int2nb[I(i,j)+0.5]==1) int2nb.erase(I(i,j));
       else int2nb[I(i,j)+0.5]--;
     }

     I(i,j)=intensity;
     if(I(i,j)>0.0)
       int2nb[I(i,j)+0.5]++;
  };

  void Station::generateApertures(){
	  int k=0;
	  int intens_old=0;
	  for(auto int_n : int2nb){
		  int intens = int_n.first;
		  for (int i=0; i < collimator.getXdim(); i++){
			  bool flag=true;
			  A[k][i] = make_pair(-1,-1);
			  for (int j=0; j<collimator.getYdim(); j++){
				  if(int(I(i,j)+0.5) == intens){
					if(flag){
						A[k][i].first=j; //left
						flag=false;
					}
					A[k][i].second=j; //right
				  }
			  }
		  }
		  intensity[k]=intens-intens_old;
		  intens_old=int_n.first;
		  k++;
	  }

	  for(;k<max_apertures;k++)
		  for (int i=0; i < collimator.getXdim(); i++)
			  A[k][i] = make_pair(-1,-1);
  }

  void Station::generateIntensityMatrix() {
    pair <int,int>aux;
    clearIntensityMatrix();
    for (int a=0 ; a<max_apertures; a++) {
      for (int i=0; i < collimator.getXdim(); i++) {
        aux = collimator.getActiveRange(i,angle);
        if ( aux.first<0 || A[a][i].first<0) continue;
        for (int j=A[a][i].first; j<=A[a][i].second; j++) {
          change_intensity(i, j, I(i,j) + intensity[a]);
        }
      }
    }
    //printIntensity();
  };

  void Station::printIntensity(bool vector_form) {
	  if(!vector_form){
		  cout << "Angle "<< angle << endl;
	    cout << " aperture intensities: ";
	    cout << intensity[0];
	    for (int i=1;i<max_apertures;i++)
	      cout << ", " << intensity[i];
	    cout << endl <<" intensity matrix:"<< endl;
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

    cout << "nb_intensities:" << int2nb.size() << endl;
  }

  void Station::setApertureShape (int a, int row, int start, int end) {
    A[a][row].first = start;
    A[a][row].second = end;
  };

  void Station::setApertureShape (int a, int row, pair<int,int> p) {
    A[a][row].first = p.first;
    A[a][row].second = p.second;
  };

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
  };

  pair<int, int> Station::getApertureShape(int a, int row) {
    return(A[a][row]);
  };

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

  /* Function to be used to get the index of a coordinate
   in the matrix I */
  int Station::getBeamIndex(pair<int,int> coord) const {
    auto pos = pos2beam.find(coord);
    return(pos->second);
  };


  void Station::writeIntensity(ifstream& myfile) {

	  for (int i=0; i<collimator.getXdim();i++) {
		  for (int j=0; j<collimator.getYdim(); j++) {
			  if(I(i,j)!=-1){
				  double intensity;
				  myfile >> intensity;
          intensity=(int) ((intensity + step_intensity/2.0)/step_intensity)*step_intensity;

				  change_intensity(i, j, (int) (intensity));
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

  string Station::toStringApertures () {
    string s = "s:" + to_string(angle) + ";";
    for (int a=0; a < max_apertures; a++) {
      s = s + "a:" + to_string(a) + ":" + to_string(intensity[a]) + ";";
      for (int i=0; i<collimator.getXdim(); i++) {
        s = s + to_string(A[a][i].first) + ":" + to_string(A[a][i].second) + ";";
      }
    }
    return(s);
  };

  string Station::toStringIntensities () {
    string s = "\"s" + to_string(angle) + "\"" + ":" + "  [";
    for (int i=0; i<collimator.getXdim();i++) {
      s += "[";
      for (int j=0; j<collimator.getYdim(); j++) {
        s = s + to_string(I(i,j));
        if(j<collimator.getYdim()-1) s+= ",";
      }
      s+= "]";
      if(i<collimator.getXdim()-1) s+=",";
    }
    s+="]";
    return(s);
  };

  void Station::printAperture(int aperture) {
    cout << "Station: "<< angle << " aperture: " << aperture << " intensity "<< intensity[aperture]<< endl;
    for (int i=0; i<collimator.getXdim(); i++) {
      cout << A[aperture][i].first << "-" << A[aperture][i].second << endl;
    }
  }

  const Matrix& Station::getDepositionMatrix(int o) const{
    const Matrix* aux=D.find(o)->second;
    return (*aux);
  }


  double Station::intensityUp(int i, int j){
    if(I(i,j)==-1) return I(i,j);
    if(getMaxIntensityRow(i) == int(I(i,j)+0.5) ||
        (j-1>=0 && I(i,j-1)>I(i,j)) || (j+1<I.nb_cols() && I(i,j+1)>I(i,j)) ){
            //next available intensity
            if(I(i,j)==0) return int2nb.begin()->first;
            else{
              map< int, int >::iterator it = int2nb.find(I(i,j)+0.5); it++;
              if(it!=int2nb.end())
                return it->first;
            }
    }

    if(int2nb.size() < max_apertures && getMaxIntensityRow(i) == int(I(i,j)+0.5) && I(i,j)<getMaxIntensity())
    	return I(i,j) + 1.0; // (getMaxIntensity()-I(i,j)+1)/2; // 1.0;



    return I(i,j);
  }

  double Station::intensityDown(int i, int j){
    if(I(i,j)<=0) return I(i,j);
    if(j==0 || j==I.nb_cols()-1 ||
        I(i,j-1)<I(i,j) || I(i,j+1)<I(i,j) ){
            //previous available intensity
            map< int, int >::iterator it = int2nb.find(I(i,j)+0.5);
            if(it!=int2nb.begin()){
              it--; return it->first;
            }else
              return (int2nb.size() < max_apertures)? int2nb.begin()->first /2 : 0.0;
     }


//    if(int2nb.size() < max_apertures && (j==0 || j==I.nb_cols()-1))
  //  	return I(i,j) - 1.0;


    return I(i,j);
  }

  list< pair< int, double > > Station::generate_intensity(double intensity, bool increment){
	  list< pair< int, double > > diff;
	  if(int2nb.size()==max_apertures) return diff;

	  map< int, int >::iterator it = int2nb.upper_bound(intensity+0.5);
	  double aux_int;
	  if(increment){
		  if(it == int2nb.begin()) aux_int=0.0;
		  else{it--; aux_int=it->first;}
	  }else
		  aux_int=it->first;

      //cout << "aux:" << aux_int << endl;
      //cout << "int:" << intensity << endl;
	  //decrement
	  if(!increment)
	  for (int i=0; i<collimator.getXdim();i++){
		  bool increasing = true;
		 /* cout << "decrement input:  " ;
		  for (int j=0; j<collimator.getYdim(); j++)
			  cout << I(i,j) << " " ;
		  cout << endl;*/

		  for (int j=0; j<collimator.getYdim(); j++){
			  if(j>0 && I(i,j)<I(i,j-1)) increasing =false;

			  if(increasing && I(i,j)==aux_int && (j==0 || I(i,j-1)!=aux_int || j==collimator.getYdim()-1)){
				  change_intensity(i, j, intensity, &diff);
			      j++; continue;
			  }


			  if(!increasing && I(i,j)==aux_int && (j==collimator.getYdim()-1 || I(i,j+1)!=aux_int))
			  	  change_intensity(i, j, intensity, &diff);
		  }

		  /*cout << "decrement output: " ;
		  for (int j=0; j<collimator.getYdim(); j++)
			  cout << I(i,j) << " " ;
		  cout << endl;*/
	  }

	  if(increment)
	  for (int i=0; i<collimator.getXdim();i++){
		  bool increasing = true;

		  /*cout << "increment input:  " ;
		  for (int j=0; j<collimator.getYdim(); j++)
			  cout << I(i,j) << " " ;
		  cout << endl;*/

		  for (int j=0; j<collimator.getYdim(); j++){
			  if(j>0 && I(i,j)<I(i,j-1)) increasing =false;

			  if(!increasing && I(i,j)==aux_int && I(i,j-1)>aux_int){
				  change_intensity(i, j, intensity, &diff);
				  continue;
			  }

			  int k=0;
			  while(increasing && I(i,j)==aux_int && j<collimator.getYdim()){
				  k++; j++;
			  }

			  if(k>0){
				  if(j<collimator.getYdim() && I(i,j)>aux_int) //increasing
					  change_intensity(i, j-1, intensity, &diff);

				  if(j==collimator.getYdim() || I(i,j)<aux_int)
					  change_intensity(i, j-(k+1)/2, intensity, &diff);
			  }
		  }

		  /*cout << "increment output: " ;
		  for (int j=0; j<collimator.getYdim(); j++)
			  cout << I(i,j) << " " ;
		  cout << endl;*/

	  }

	  return diff;
  }

  list< pair< int, double > > Station::change_intensity(double intensity, double delta){
	  list< pair< int, double > > diff;

	  for (int i=0; i<collimator.getXdim();i++)
		  for (int j=0; j<collimator.getYdim(); j++)
			  if(I(i,j)==intensity)
				  change_intensity(i, j, I(i,j)+delta, &diff);

	  return diff;
  }

  void Station::diff_undo(list< pair< int, double > >& diff){
	  for(pair< int, double > d:diff){
		  int i = beam2pos[d.first].first;
		  int j = beam2pos[d.first].second;
		  change_intensity(i, j, I(i,j)-d.second);
	  }
	  diff.clear();
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

  void Station::reduce_apertures(list< pair< int, double > >& diff){
	  if(int2nb.size()==0) return;

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
    while(int2nb.size()>max_apertures)
    	reduce_apertures(diff);


    return diff;
  }

  bool Station::isOpenBeamlet(int beam, int aperture) {
  	auto coord = collimator.indexToPos(beam, angle);
        if (!collimator.isActiveBeamAngle(coord.first, coord.second,angle)) return(false);
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

    list<pair <int, double> > diff;
    clearHistory();

    if (isActiveBeamlet(beam) && isOpenBeamlet(beam, aperture)) {
      auto coord = collimator.indexToPos(beam, angle);
      int row = coord.first;

      last_mem = make_pair(make_pair(aperture,row), A[aperture][row]);

      if (lside) {
        // Close from the left side

        for (int i=0; i<=coord.second-A[aperture][row].first; i++) {
          auto coord2 = collimator.indexToPos(beam-(coord.second-A[aperture][row].first)+i, angle);
          diff.push_back(make_pair(beam-(coord.second-A[aperture][row].first)+i, -intensity[aperture]));
          change_intensity(coord2.first, coord2.second, I(coord2.first, coord2.second) - intensity[aperture]);
          //I(coord2.first,coord2.second) = I(coord2.first,coord2.second) - intensity[aperture];
        }

        if (A[aperture][row].second == coord.second) {
          A[aperture][row].first = -1;
          A[aperture][row].second = -1;
        } else {
          A[aperture][row].first = coord.second+1;
        }
      } else {
        // Close from the right side

        for (int i=0;i<=A[aperture][row].second-coord.second;i++) {
          auto coord2 = collimator.indexToPos(beam+i, angle);
          diff.push_back(make_pair(beam+i, -intensity[aperture]));
          change_intensity(coord2.first, coord2.second, I(coord2.first, coord2.second) - intensity[aperture]);
          //I(coord2.first,coord2.second) = I(coord2.first,coord2.second) - intensity[aperture];
	      }

	      if (A[aperture][row].first == coord.second) {
	        A[aperture][row].first=-1;
	        A[aperture][row].second=-1;
	      } else {
          A[aperture][row].second=coord.second-1;
        }
      }
      //updateIntensity(diff);
      last_diff=diff;
    }

    return(diff);
  }

  /* Function that closes a beamlet in coordinate <x,y> from the left, if lside is true, or
   from the right size otherwise. Return true if the closing was performed.*/
  list<pair <int,double> > Station::closeBeamlet(pair<int,int> coord, int aperture, bool lside) {

    list<pair <int, double> > diff;
    int beam = getBeamIndex(coord);

    clearHistory();


    if (isActiveBeamlet(beam) && isOpenBeamlet(beam, aperture)) {
      auto coord = collimator.indexToPos(beam, angle);
      int row= coord.first;

      last_mem = make_pair(make_pair(aperture,row), A[aperture][row]);

      if (lside) {
        // Close from the left side

        for (int i=0; i<=coord.second-A[aperture][row].first; i++) {
          auto coord2 = collimator.indexToPos(beam-(coord.second-A[aperture][row].first)+i, angle);
          diff.push_back(make_pair(beam-(coord.second-A[aperture][row].first)+i, -intensity[aperture]));
          change_intensity(coord2.first, coord2.second, I(coord2.first, coord2.second) - intensity[aperture]);
          //I(coord2.first,coord2.second) = I(coord2.first,coord2.second) - intensity[aperture];
        }

        if (A[aperture][row].second == coord.second) {
          A[aperture][row].first=-1;
          A[aperture][row].second=-1;
        } else {
          A[aperture][row].first=coord.second+1;
        }
      } else {
        // Close from the right side

        for (int i=0;i<=A[aperture][row].second-coord.second;i++) {
          auto coord2 = collimator.indexToPos(beam+i, angle);
          diff.push_back(make_pair(beam+i, -intensity[aperture]));
          change_intensity(coord2.first, coord2.second, I(coord2.first, coord2.second) - intensity[aperture]);
          //I(coord2.first,coord2.second) = I(coord2.first,coord2.second) - intensity[aperture];
        }

        if (A[aperture][row].first == coord.second) {
          A[aperture][row].first=-1;
          A[aperture][row].second=-1;
        } else {
          A[aperture][row].second=coord.second-1;
        }
      }
      //updateIntensity(diff);
      last_diff = diff;
    }
    return(diff);
  }

  /* Function that opens a beamlet from the left, if lside is true, or
     from the right size otherwise. Return true if the closing was performed.*/
  list <pair<int,double> > Station::openBeamlet(int beam, int aperture) {
    list<pair<int, double>> diff;
    clearHistory();

    if (isActiveBeamlet(beam) && !isOpenBeamlet(beam, aperture)) {
      auto coord = collimator.indexToPos(beam, angle);
      int row = coord.first;

      // Record the configuration of the row to be open
      last_mem = make_pair(make_pair(aperture,row), A[aperture][row]);

      if (A[aperture][row].first < 0) {
        //When the row is completely closed
        diff.push_back(make_pair(beam, intensity[aperture]));
        change_intensity(coord.first, coord.second, I(coord.first, coord.second) + intensity[aperture]);
        //I(coord.first,coord.second) = I(coord.first,coord.second) + intensity[aperture];

        A[aperture][row].first=coord.second;
        A[aperture][row].second=coord.second;
      } else {
        // Row has open beamlets
        if (A[aperture][row].first > coord.second) {
          for (int i=0;i<(A[aperture][row].first-coord.second);i++) {
            auto coord2 = collimator.indexToPos(beam+i, angle);
            diff.push_back(make_pair(beam+i, intensity[aperture]));
            change_intensity(coord2.first, coord2.second, I(coord2.first,coord2.second) + intensity[aperture]);
            //I(coord2.first,coord2.second) = I(coord2.first,coord2.second) + intensity[aperture];
          }
          A[aperture][row].first = coord.second;
        } else {
          for (int i=0;i<(coord.second-A[aperture][row].second);i++) {
            auto coord2 = collimator.indexToPos(beam-(coord.second-A[aperture][row].second)+1+i, angle);
            diff.push_back(make_pair(beam-(coord.second-A[aperture][row].second)+1+i, intensity[aperture]));
            change_intensity(coord2.first, coord2.second, I(coord2.first,coord2.second) + intensity[aperture]);
            //I(coord2.first,coord2.second) = I(coord2.first,coord2.second) + intensity[aperture];
          }
          A[aperture][row].second = coord.second;
        }
      }
		 	//updateIntensity(diff);
		  //printIntensity(true);
      last_diff = diff;
		}

		return(diff);
  }

  list <pair< int,double> > Station::getModifyIntensityApertureDiff (int aperture, double size) {
    list < pair <int, double > > diff;
    if ((intensity[aperture]+size) < 0 || (intensity[aperture]+size) > max_intensity) {
      if (intensity[aperture]+size < 0) {
        // Too low intensity
        if (intensity[aperture]>0)
          size = intensity[aperture];
        else
          return(diff);
      } else if ((intensity[aperture]+size) > max_intensity) {
        // Too high intensity
        if (intensity[aperture] < max_intensity)
          size = max_intensity - intensity[aperture];
        else
          return(diff);
      } else {
        return (diff);
      }
    }
    for (int i=0; i<collimator.getXdim(); i++) {
      if (A[aperture][i].first<0 || A[aperture][i].second<0)
        continue;
      int beamlet = pos2beam[make_pair(i, A[aperture][i].first)];
      for(int j=A[aperture][i].first; j<=A[aperture][i].second; j++){
        diff.push_back(make_pair(beamlet, size));
        beamlet++;
      }
    }
    return(diff);

  }

  list <pair< int,double> > Station::modifyIntensityAperture(int aperture, double size) {
    list < pair <int, double > > diff;

    // Since we will apply a change (even if it is not possible)
    // we clear previous changes.
    clearHistory();

    // Check intensity bounds
    if ((intensity[aperture]+size) < 0 || (intensity[aperture]+size) > max_intensity) {
      if (intensity[aperture]+size < 0) {
        // Too low intensity
        if (intensity[aperture]>0)
          size = intensity[aperture];
        else
          return(diff);
      } else if ((intensity[aperture]+size) > max_intensity) {
        // Too high intensity
        if (intensity[aperture] < max_intensity)
          size = max_intensity - intensity[aperture];
        else
          return(diff);
      } else {
        return (diff);
      }
    }

    // Save previous intensity
    for (int i=0; i<max_apertures; i++)
      last_intensity.push_back(intensity[i]);

    // Increase intensity
    intensity[aperture] = intensity[aperture] + size;


    // Create diff and update
    for (int i=0; i<collimator.getXdim(); i++) {
      if (A[aperture][i].first<0 || A[aperture][i].second<0)
        continue;
      int beamlet = pos2beam[make_pair(i, A[aperture][i].first)];
      for(int j=A[aperture][i].first; j<=A[aperture][i].second; j++){
        auto coord = collimator.indexToPos(beamlet, angle);
        change_intensity(coord.first, coord.second, I(coord.first, coord.second) + size);
        diff.push_back(make_pair(beamlet, size));
        beamlet++;
      }
    }

    last_diff = diff;
    return (diff);
  }

  void Station::updateIntensity(list<pair<int,double> > diff) {
    pair<int,int> coord;
    if (diff.size()<1){
      return;
    }
    for (auto it=diff.begin(); it!=diff.end(); it++) {
      coord = collimator.indexToPos(it->first, angle);
      change_intensity(coord.first, coord.second,  I(coord.first,coord.second) + it->second);
      //I(coord.first,coord.second) = I(coord.first,coord.second) + it->second;
      if (I(coord.first,coord.second) < 0) {
        cout <<"Error intensidad negativa en updateIntensity!!!!" << endl;
        cout << "  Cordenadas matriz " << coord.first << ","<<coord.second << endl;
        cout << "  Intensidad actual " << I(coord.first,coord.second)-it->second << endl;
        cout << "  Intensity change "<< it->second << endl;
        //cout << "  Intensidad de la apertura "<< getApertureIntensity(0) << endl;
        getchar();
      }
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

  void Station::setApertureIntensity(int aperture, double value) {
    intensity[aperture]=value;
  };

  double Station::getApertureIntensity(int aperture) {
    return(intensity[aperture]);
  };

  int Station::getMaxIntensity() {
    return(max_intensity);
  };

  list <pair<int,double> > Station::undoLast () {
    list <pair<int,double> > undo_diff;
    pair<int,int> coord;

    if (last_diff.size() < 1) {
      clearHistory();
      return(undo_diff);
    }

    //Last mem is not empty
    if (last_mem.first.first>=0){
      A[last_mem.first.first][last_mem.first.second].first = last_mem.second.first;
      A[last_mem.first.first][last_mem.first.second].second = last_mem.second.second;
    }

    for (list<pair<int,double> >::iterator it=last_diff.begin();it!=last_diff.end();it++) {
      coord=collimator.indexToPos(it->first, angle);
      change_intensity(coord.first, coord.second,  I(coord.first,coord.second) - it->second);
      //I(coord.first,coord.second) = I(coord.first,coord.second) - (it->second);
      undo_diff.push_back(make_pair(it->first, -(it->second)));
    }

    for (int i=0;i<max_apertures && last_intensity.size() > 0;i++) {
      intensity[i] = last_intensity[i];
    }

    clearHistory();

    return(undo_diff);
  }

  void Station::clearHistory() {
    last_mem = make_pair(make_pair(-1,-1), make_pair(-1,-1));
    last_diff.clear();
    last_intensity.clear();
    if (last_intensity.size()>1 && last_diff.size()>1)
     cout << "Error: elementos de undo no se borraron" << endl;
  };
}
