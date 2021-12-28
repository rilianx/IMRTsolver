#include <iostream>
#include <iterator>
#include <set>
#include <stack>
#include <unistd.h>
#include <string.h>
#include "Collimator.h"
#include "Volume.h"
#include "IntensityGenerator.h"

using namespace std;
using namespace imrt;

namespace imrt{

vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}


void load_scores(list<Score>& scores, ifstream& indata){
    while ( !indata.eof() ) {
      string type;
      double x, min_value, max_value, weight;
      int organ;
      Score::Type t;
      indata >> type;
      if (type[0]=='#'){
        string aux;
        getline(indata,aux);
        continue;
      }

      indata >>  x >> organ >> min_value >> max_value >> weight;
      //cout <<  type << " " << x  << " "<< organ << " " << min_value <<  " " << max_value << " " <<  weight << endl;
      if (type=="D") t = Score::D; 
      else if (type=="V") t=Score::V;
      else if (type=="Dmean") t= Score::Dmean;

      scores.push_back(Score(t,x,organ,min_value,max_value,weight));
    }
}

void load_scores(list<Score>& scores, string scores_file){
    ifstream indata; // indata is like cin
    indata.open(scores_file); 
    load_scores(scores, indata);
}


set<int> get_angles(string file, int n){
  ifstream _file(file.c_str(), ios::in);
  string line;

  if (! _file)
    cerr << "ERROR: unable to open instance file: " << file << ", stopping! \n";

  //cout << "##Reading volume files." << endl;
  getline(_file, line);
  _file.close();

  set<int> angles;
  stack<string> q;
  char delim=' ';
  std::size_t current, previous = 0;
  current = line.find(delim);
  //cout << line << endl;
  while (current != std::string::npos) {
    angles.insert(atoi(line.substr(previous, current - previous).c_str()));
    previous = current + 1;
    current = line.find(delim, previous);
  }
  angles.insert(atoi(line.substr(previous, current - previous).c_str()));

  return angles;
}


vector<Volume> createVolumes (string organ_filename, Collimator& collimator){
  ifstream organ_file(organ_filename.c_str(), ios::in);
  vector<string> organ_files;
  vector<Volume> volumes;
  string line;

  if (! organ_file)
    cerr << "ERROR: unable to open instance file: " <<
             organ_filename << ", stopping! \n";

  cout << "##Reading volume files." << endl;
  getline(organ_file, line);
  while (organ_file) {
    getline(organ_file, line);
    if (line.empty()) continue;
    cout << "##  " << line << endl;
    //Assuming one data point
    organ_files.push_back(line);
  }
  organ_file.close();
  cout << "##  Read " << organ_files.size() << " files"<< endl;

  for (int i=0; i<organ_files.size(); i++)
    volumes.push_back(Volume(collimator, organ_files[i]));

  return(volumes);
}
}