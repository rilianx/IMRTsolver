/*
 * Volume.cpp
 *
 *  Created on: 22 may. 2018
 *      Author: leslie
 */


#include <regex>

#include <regex>
#include "Volume.h"

namespace imrt {

list<int> get_angles(string str) 
{ 
    list<int> angles;
    std::replace(str.begin(), str.end(), '_', ' ');
    std::replace(str.begin(), str.end(), '-', ' ');
    std::replace(str.begin(), str.end(), '/', ' ');
    stringstream ss;     
  
    /* Storing the whole string into string stream */
    ss << str; 
  
    /* Running loop till the end of the stream */
    string temp; 
    int found; 
    bool flag=false;
    while (!ss.eof()) { 
  
        /* extracting word by word from stream */
        ss >> temp; 
  
        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found) {
            angles.push_back(found); 
            flag=true;
        }else if(flag==true){
          //stringstream(temp) >> organ_name;
          //flag=false;
        }
  
        /* To save from space at the end of string */
        temp = ""; 
    } 
    return angles;
} 


Volume::Volume(Collimator& collimator, string deposition_file, int max_voxels_per_organ) :
		collimator(collimator) {

    string file_map_voxels=deposition_file;
    file_map_voxels = std::regex_replace(file_map_voxels, std::regex("_DDM_"), "_voxelIndex");
    file_map_voxels = std::regex_replace(file_map_voxels, std::regex("-DDM_"), "-voxelIndex_");
    file_map_voxels = std::regex_replace(file_map_voxels, std::regex(".dat"), ".txt");

    if(deposition_file!="")
      set_data(deposition_file, file_map_voxels, max_voxels_per_organ);
}

void Volume::add_data(string file){
  list<int> angles = get_angles(file);


  string file_map_voxels=file;
  file_map_voxels = std::regex_replace(file_map_voxels, std::regex("_DDM_"), "_voxelIndex");
  file_map_voxels = std::regex_replace(file_map_voxels, std::regex("-DDM_"), "-voxelIndex_");
  file_map_voxels = std::regex_replace(file_map_voxels, std::regex(".dat"), ".txt");
  
  set_data(file, file_map_voxels, 0 , angles);
}

void Volume::add_data(string depo_file, string map_voxels_file, int angle){
  list<int> angles;
  angles.push_back(angle);
  set_data(depo_file, map_voxels_file, 0 , angles);
}

void Volume::set_data(string file, string file_map_voxels, int max_voxels_per_organ, list<int> angles) {
  cout << "reading file " << file << "..." <<endl;
  cout << "reading file " << file_map_voxels << "..." <<endl;
  string line, linec;
  ifstream myfile (file);
  ifstream indexfile (file_map_voxels);
  stringstream ss;
  double aux1, aux2;
  vector <pair<double,double> >::iterator it;

  if(angles.size()==0) angles=collimator.getAngles();


  //nb_beamlets=-1;
  nb_voxels=-1;

  if (!myfile.is_open() || !indexfile.is_open())
    throw runtime_error("error reading file.");

  // Get dimensions of the matrix voxels/beamlets matrix
  vector<string> lines;
  // First line not considered
  // TODO: read name of the volume!
  getline (myfile,line);
  // Initial reading
  map<int,int> line2row;
  while (getline (myfile,line) ){
    int lin, id, row;
    indexfile >> lin; indexfile >> id; indexfile >> row; 
    line2row[row]=lin-1;
    //cout << lin << "," << row << endl;
    lines.push_back(line);
  }
  
  myfile.close();
  nb_voxels = lines.size()-1;
  //cout << nb_voxels << endl;

  
  double step=1;
  if(max_voxels_per_organ > 0 && max_voxels_per_organ < nb_voxels){
    step = (double) nb_voxels/ (double) max_voxels_per_organ;
    nb_voxels = max_voxels_per_organ;
  }
  
  ss.str(lines[0]);

  for(auto angle:angles)
    D[angle]=Matrix(nb_voxels, collimator.getNangleBeamlets(angle));
  

  auto ang = angles.begin();
  for (int i=0; i<nb_voxels; i++) {
    
    
    int ii = step* (double) i;
    ss.clear(); ss.str(lines[ii]);
    getline(ss, line, '\t');
    int a=0, j=0;
    ang = angles.begin();
    int count=0; 
    bool flag=false;

    

    while (getline(ss, line, '\t') ) {
      flag=true;
      if (j >= collimator.getNangleBeamlets( *ang )) {
         a++; ang++; j=0; 
      }

      D[ *ang ](line2row[i],j)=atof(line.c_str());

      D[ *ang ](line2row[i],j)=atof(line.c_str());
      j++;
    }
    
    
  }

  //print some cells
  // cout << "Matrix of deposition for angle " << *angles.begin() << endl;
  // for (int i=1000;i<1010;i++) {
  //   for (int j=0;j<20;j++)
  //     cout << D[*angles.begin()](i,j) << ",";
  //   cout << endl;
  // }

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
