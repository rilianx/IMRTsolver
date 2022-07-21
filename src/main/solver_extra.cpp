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

double iterated_local_search(Collimator& collimator, vector<Volume>& volumes, int max_apertures,
       int max_intensity, int step_intensity,
       StationSetup setup, list<int> bac, vector<Evaluator*>& evaluators, int sf_eval, int of_eval,
      int maxeval, int& used_evaluations, vector<NeighborhoodType>& neighborhoods, int perturbation_size,
      ofstream& output_stream, std::clock_t begin_time, double min_delta_eval, double alpha,
      int switch_patience, vector<double> pr_neigh, bool _verbose){

    // Create an initial plan
    Plan P (collimator, volumes, max_apertures,
            max_intensity, 0, step_intensity,
            setup, NULL, bac); 

    double best_eval = evaluators[sf_eval]->eval(P);
    cout << "## Initial solution: " << best_eval << endl;
 
    //IBO_LS
    ILS* ils = new IntensityILS2(evaluators, sf_eval, of_eval);
    
    double cost = ils->iteratedLocalSearch(P, maxeval, neighborhoods, perturbation_size, 
                    output_stream,used_evaluations, begin_time, _verbose, min_delta_eval, alpha, 
                    switch_patience, pr_neigh, 1.0);
    
    cout << "## Best solution found: " << ils->best_evals[of_eval] << endl;
    evaluators[0]->eval(P);
    return ils->best_evals[of_eval];
}

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


string get_organ_name(string str) 
{
    list<int> angles;
    std::replace(str.begin(), str.end(), '_', ' ');
    std::replace(str.begin(), str.end(), '-', ' ');
    std::replace(str.begin(), str.end(), '/', ' ');
    stringstream ss;     
    string organ_name;
  
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
            flag=true;
        }else if(flag==true){
          stringstream(temp) >> organ_name;
          
        }
  
        /* To save from space at the end of string */
        temp = ""; 
    }

    return organ_name;
} 

vector<Evaluator*> createEvaluators(Collimator& collimator,  vector<Volume>& volumes, string evaluators_str ){
    vector<double> w={1,1,5};
    vector<double> Zmin={0,0,76};
    vector<double> Zmax={65,65,76};

    FluenceMap* fm = new FluenceMap(volumes, collimator);

    vector<Evaluator*> evaluators;
    vector<std::string> files = split(evaluators_str, ',');

    for (string file: files){
        cout << "reading:" << file << endl;
        ifstream indata; // indata is like cin
        indata.open(file); 

        string function; indata >> function;
        list<Score> scores;
        load_scores(scores, indata);

        EvaluatorGS::Type t;
        cout << function << endl;
        if( function == "gs_relu") t=EvaluatorGS::GS_RELU;
        else if( function == "gs") t=EvaluatorGS::GS;
        else if( function == "gs2") t=EvaluatorGS::GS2;
        
        evaluators.push_back(new EvaluatorGS(*fm,w,Zmin,Zmax, scores, t));
    }
    evaluators.push_back(new EvaluatorF(*fm,w,Zmin,Zmax));

    return evaluators;
}


vector<Volume> createVolumes (string organ_filename, Collimator& collimator){
  ifstream organ_file(organ_filename.c_str(), ios::in);
  vector<string> depo_files;
  string line;

  if (! organ_file)
    cerr << "ERROR: unable to open instance file: " <<
             organ_filename << ", stopping! \n";

  cout << "##Reading volume files." << endl;
  getline(organ_file, line); //line of initial angles
  while (organ_file) {
    getline(organ_file, line);
    if (line.empty()) continue;
    cout << "##  " << line << endl;
    //Assuming one data point
    depo_files.push_back(line);
  }
  organ_file.close();
  cout << "##  Read " << depo_files.size() << " files"<< endl;

  map< string, Volume*> volumes_map;
  vector<Volume> volumes;

  for (int i=0; i<depo_files.size(); i++){
    string organ_name = get_organ_name(depo_files[i]);
    if (volumes_map.find(organ_name)!=volumes_map.end()) volumes_map[organ_name]->add_data(depo_files[i]);
    else {
      volumes.push_back(Volume(collimator, "", 0));
      volumes_map.insert(make_pair(organ_name, &volumes.back() ));
      volumes_map[organ_name]->add_data(depo_files[i]);
    }
  }
    
  return(volumes);
}


/*
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
*/
}