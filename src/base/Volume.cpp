/*
 * Volume.cpp
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */

#include "Volume.h"

namespace imrt {

Volume::Volume(Collimator& collimator, string deposition_file) :
		collimator(collimator) {
  set_data(deposition_file);
}

void Volume::set_data(string file) {
  string line, linec;
  ifstream myfile (file);
  stringstream ss;
  double aux1, aux2;
  vector <pair<double,double> >::iterator it;

  //nb_beamlets=-1;
  nb_voxels=-1;

  if (!myfile.is_open())
    throw runtime_error("error reading file.");

  // Get dimensions of the matrix voxels/beamlets matrix
  vector<string> lines;
  // First line not considered
  // TODO: read name of the volume!
  getline (myfile,line);
  // Initial reading
  while (getline (myfile,line) )
    lines.push_back(line);
  myfile.close();
  nb_voxels = lines.size()-1;
  ss.str(lines[0]);

  for(auto angle:collimator.getAngles())
    D[angle]=Matrix(nb_voxels, collimator.getNangleBeamlets(angle));

  for (int i=0; i<nb_voxels; i++) {
    ss.clear(); ss.str(lines[i]);
    getline(ss, line, '\t');
    int a=0, j=0;
    int count=0;
    while (getline(ss, line, '\t') ) {
      if (j >= collimator.getNangleBeamlets( collimator.getAngle(a) )) {
        a++; j=0;
      }
      D[ collimator.getAngle(a) ](i,j)=atof(line.c_str());
      j++;
    }
  }

}

/*//TODO: this should be automatically computed in the future
bool Volume::set_data(string file, vector<string> angle_coord_files){
    string line, linec;
    ifstream myfile (file);
    stringstream ss;
    double aux1, aux2;
    vector <pair<double,double> >::iterator it;

    nb_angles=angle_coord_files.size();
    nb_beamlets=-1;
    nb_voxels=-1;

    if (!myfile.is_open())
      throw runtime_error("error reading file.");

    // Get dimensions of the matrix voxels/beamlets matrix
    vector<string> lines;
    // Initial reading
    while (getline (myfile,line) )
      lines.push_back(line);
    myfile.close();
    nb_voxels = lines.size();
    ss.str(lines[0]);
    while (getline(ss, line, '\t') )
      nb_beamlets++;

    // Reading each angle to obtain number of beamlets
    // and their position in the beam array

    int ccol=0, nbeam=0;
    bool first;
    double prev_row;

    for (int i=0; i<angle_coord_files.size(); i++) {
      nbeam=nbeam+ccol; ccol=0; first=true;
      prev_row=999999.0;

      myfile.open(angle_coord_files[i]);
      if (!myfile.is_open())
        throw runtime_error("error reading file.");

      while (getline (myfile,line)) {
        ss.clear(); ss.str(line);
        getline (ss,linec,'\t');
        aux1 = atof(linec.c_str());
        getline (ss,linec,'\t');
        aux2 = atof(linec.c_str());

        if (first) {
          prev_row = (double) aux1;
          first=false;
        }

        if (prev_row!=aux1) {
          R[i].push_back(make_pair(nbeam,nbeam+ccol-1));
          //r_active_size[i].push_back(ccol);
          nbeam=nbeam+ccol;
          ccol=1;
          prev_row= (double)aux1;
        } else {
          ccol++;
        }
      }
      R[i].push_back(make_pair(nbeam,nbeam+ccol-1));
      //r_active_size[i].push_back(ccol);
      myfile.close();
    }
    //print_r_active_size();
    //print_coordinates();


    //TODO: assuming only one angle!
    D[1] = Matrix(nb_voxels,nb_beamlets);
    // Save matrix
    for (int i=0; i<nb_voxels; i++) {
      ss.clear(); ss.str(lines[i+1]);
      int j=0;
      while (getline(ss, line, '\t') ) {
        if(j!=0)
          D[1](i,j-1) = atof(line.c_str());
        j++;
      }
    }

    //print_deposition();
}*/

Matrix& Volume::getDepositionMatrix(int angle){
  return(D[angle]);
}

void Volume::print_deposition() {
  for(auto d:D){
  //for (int a=0; a< nb_angles; a++) {
    int angle=d.first;
    Matrix* dm=&d.second;
    cout << "ANGLE: " << angle << endl;
    for (int i=0;i<nb_voxels;i++) {
      for (int j=0;j<collimator.getNangleBeamlets(angle);j++)
        cout << (*dm)(i,j) << ",";
      cout << endl;
    }
  }
}

void Volume::print_coordinates() {
/*  for (int i=0; i<nb_angles; i++) {
    cout << "Angle: " << i <<endl;
    for (int j=0; j<R[i].size();j++){
      cout << R[i][j].first <<","<< R[i][j].second << endl;
    }
  }*/
}

/*void Volume::print_r_active_size() {
  for (int i=0; i<nb_angles; i++) {
    cout << "Angle: " << i <<endl;
    for (int j=0; j< r_active_size[i].size();j++){
      cout << j<<":"<<r_active_size[i][j] << endl;
    }
  }
}*/

}
