/*
 * Colimador.cpp
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */

#include "Collimator.h"

namespace imrt {
  // There is a coordinate file per angle per angle
Collimator::Collimator(vector < pair<int, string> >& coord_files) {
    ifstream myfile;
    string line, linec;
    stringstream ss;
    double x, y;
    bool flag;

    //nb_angles=coord_files.size();

    for (int i=0 ; i < coord_files.size(); i++) {
      // Open file
      int angle = coord_files[i].first;
      angles.push_back(angle);
      myfile.open(coord_files[i].second);
      if (!myfile.is_open())
        throw runtime_error("error reading file.");

      // Read lines of the file
      while (getline (myfile,line)) {
        flag=false;
        ss.clear(); ss.str(line);
        getline (ss,linec,'\t'); //the first number is discarded (beamlet id)
        getline (ss,linec,'\t');
        x = atof(linec.c_str());
        getline (ss,linec,'\t');
        y = atof(linec.c_str());


        angle_coord[angle].push_back(make_pair(x,y));

        if (beam_coord.find(x) == beam_coord.end()) {
          // New x coordinate
          beam_coord[x].push_back(y);
          insertXorder(x);
          insertYorder(y);
        } else {
          if (find(beam_coord[x].begin(), beam_coord[x].end(), y) != beam_coord[x].end())
            flag=true;
          // New y coordinate: insert in order
          for (int j=0; j < beam_coord[x].size() && !flag; j++) {
            if (beam_coord[x][j] > y) {
              beam_coord[x].insert(beam_coord[x].begin()+j, y);
              flag=true;
            }
          }
          if (!flag) beam_coord[x].push_back(y);
          insertYorder(y);
        }
      }
      myfile.close();
      nb_angle_beamlets[angle]=(angle_coord[angle].size());
    }
    xdim = xcoord.size();
    ydim = ycoord.size();
    nb_beamlets= xdim*ydim;
    //printCoordinates();

    setActiveBeamlets(angle_coord);
    //printActiveBeam();
  }

  // Active beamlets per angle
  void Collimator::setActiveBeamlets(map<int,  vector<pair<double,double> > >& coord){
    map<int, map<string, vector<double> > > ::iterator it;
    vector<double> aux;
    double nmax, nmin;
    int selmax, selmin;
    bool flag;

    for (auto angle : coord){
      int i=angle.first;
    //for (int i=0; i<nb_angles; i++) {
      // Angle cycle
      for (int j=0; j<xcoord.size(); j++) {
        flag=false;
        // Row cycle
        nmax=-9999999; nmin=9999999;
        for (int s=0; s<coord[i].size(); s++) {
          if (angle.second[s].first == xcoord[j]) {
            if (angle.second[s].second<nmin)
              nmin=angle.second[s].second;
            if (angle.second[s].second>nmax)
              nmax=angle.second[s].second;
            flag=true;
          }
        }

        if (flag) {
          for (int s=0; s<ycoord.size(); s++) {
            if (ycoord[s]==nmin) selmin=s;
            else if(ycoord[s]==nmax) selmax=s;
          }
        } else {
          selmin=-1;selmax=-1;
        }
        angle_row_beam_active[i].push_back(make_pair(selmin, selmax));
      }
    }
  }

  // Insert a new X coordinate in the xcoord vector
  void Collimator::insertXorder(double x) {
    bool flag=false;
    if (find(xcoord.begin(), xcoord.end(), x) != xcoord.end())
      flag=true;
    for (int j=0; j< xcoord.size() && !flag; j++) {
      if (xcoord[j] > x) {
        xcoord.insert(xcoord.begin()+j, x);
        flag=true;
      }
    }
    if (!flag) xcoord.push_back(x);
  }

  // Insert a new Y coordinate in the ycoord vector
  void Collimator::insertYorder(double y) {
    bool flag=false;
    if (find(ycoord.begin(), ycoord.end(), y) != ycoord.end())
      flag=true;
    for (int j=0; j< ycoord.size() && !flag; j++) {
      if (ycoord[j] > y) {
        ycoord.insert(ycoord.begin()+j, y);
        flag=true;
      }
    }
    if (!flag) ycoord.push_back(y);
  }

  int Collimator::getXdim() {
    return(xdim);
  }

  int Collimator::getYdim() {
    return(ydim);
  }

  int Collimator::getNangleBeamlets(int angle) {
    return(nb_angle_beamlets[angle]);
  }

  pair<int,int> Collimator::indexToPos(int index, int angle){
    double x= angle_coord[angle][index].first;
    double y= angle_coord[angle][index].second;

    ptrdiff_t posx = find(xcoord.begin(), xcoord.end(), x) - xcoord.begin();
    ptrdiff_t posy = find(ycoord.begin(), ycoord.end(), y) - ycoord.begin();
    return(make_pair(posx, posy));
  }

  void Collimator::printCoordinates() {
    for (map<double, vector<double> >::iterator it = beam_coord.begin() ; it != beam_coord.end() ; it ++ ) {
      for (int i=0; i < it->second.size(); i++) {
        cout << it->first << " : " << it->second[i] << endl;
      }
    }
  }

  void Collimator::printAxisValues() {
    cout << "x:";
    for (int i=0;i<xcoord.size(); i++ )
      cout << xcoord[i] << ",";
    cout << endl;
    cout << "y:";
    for (int i=0;i<ycoord.size(); i++ )
      cout << ycoord[i] << ",";
    cout << endl;
  }

  void Collimator::printActiveBeam() {
    for(auto angle:angle_row_beam_active){
      cout << "Angle: " << angle.first << endl;
      for (int i =0; i < angle.second.size(); i++) {
         cout << angle.second[i].first << " , " << angle.second[i].second<< endl;
      }
    }
  }

  bool Collimator::isActiveBeamAngle(int x, int y, int angle) {
    if (angle_row_beam_active[angle][x].first <= y && angle_row_beam_active[angle][x].second >= y)
      return(true);
    return(false);
  }

  pair<int, int> Collimator::getActiveRange(int x, int angle) {
    return(angle_row_beam_active[angle][x]);

  }

  list<int>& Collimator::getAngles(){
    return angles;
  }

  int Collimator::getAngle(int i){
    std::list<int>::iterator it = angles.begin();
    std::advance(it, i);
    return *it;
  }

}
