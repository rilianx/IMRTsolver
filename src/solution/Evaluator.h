/*
 * EvaluatorF.h
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <iterator>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <queue>
#include <functional>

#include "Plan.h"
#include "Volume.h"
#include "Matrix.h"
#include "Station.h"

#ifndef Evaluator_H_
#define Evaluator_H_

using namespace std;
using namespace maths;


namespace imrt {

class Plan; 

struct pair_hash
{
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};


class FluenceMap{

public:
	FluenceMap(vector<Volume>& volumes, const Collimator& collimator);
	virtual ~FluenceMap();
	void create_voxel2beamlet_list(vector<Volume>& volumes, const Collimator& collimator);

	// Generate the dose distribution matrices Z for each organ
	void computeFM(const Plan& p);

	//Get deltaFM related to the list of changes
	std::vector<list < pair<int, double>>>& compute_deltaFM(list< pair< int, double > >& changes, double angle) const;
	std::vector<list < pair<int, double>>>& get_deltaFM() const {return deltaFM;}
	void updateFM(list< pair< int, double > >& changes, double angle, std::vector<list < pair<int, double>>>& deltaFM);
	
	const vector< vector<double> >& getFM() {return FM;}

private:
	//dose distribution vectors for each organ
	vector< vector<double> > FM;
	mutable std::vector<list < pair<int, double>>> deltaFM;

	//voxel_dose[o][d] indicates the number of voxels in organ o with a dose between d and d+1
	vector<vector<double> > voxel_dose;

  	vector<Volume>& volumes;

	//number of organs, including the tumor
	int nb_organs;

	//number of voxels for each organ
	vector<int> nb_voxels;

  	//voxels[angle][o,k] --> beamlet list (b)
  	unordered_map < int, unordered_map < pair<int, int>, list <int> , pair_hash > > voxel2beamlet_list;

  	//beamlet [angle][b] --> lista de voxels (o,k)
  	unordered_map < int, unordered_map <int, multimap<double, pair<int,int> > > > beamlet2voxel_list;

};


 struct MagnitudeCompare2
 {
     bool operator()(const double& lhs, const double& rhs) const
     {
         return abs(lhs) > abs(rhs);
     }

     bool operator()(const pair< double, pair<int,int> >& lhs, const pair< double, pair<int,int> >& rhs) const
     {
    	 if(abs(lhs.first) > abs(rhs.first)) return true;
    	 if(abs(lhs.first) == abs(rhs.first)){
    		 //if(pair_hash()(lhs.second)< pair_hash()(rhs.second)) return true;
    		 if(lhs.second.first < rhs.second.first) return true;
    		 if(lhs.second.first == rhs.second.first &&  lhs.second.second < rhs.second.second) return true;
    	 }
         return false;
     }
};

/* Abstract Evaluator
 */


class Evaluator{


public:

	//Constructor of the evaluator.
	Evaluator(FluenceMap& fm_structure, vector<double>& w, vector<double>& Zmin, 
    vector<double>& Zmax) : fm_structure(fm_structure), Zmin(Zmin), Zmax(Zmax), w(w), FM(fm_structure.getFM()) { };


    FluenceMap& fm_structure;
    const vector< vector<double> >& FM;

    vector<double>& w;
    vector<double>& Zmin;
    vector<double>& Zmax;


	virtual ~Evaluator() { }

	// Eval the cost F based on the dose deposition matrix Z
	virtual double eval(const Plan& p)=0;

	virtual double incremental_eval()=0;

	virtual double incremental_eval(list< pair< int, double > >& changes, double angle)=0;

    double get_delta_eval(int angle, int b, double delta_intensity) const{
			  list< pair< int, double > > changes;
			  changes.push_back(make_pair(b,delta_intensity));
			  return  get_delta_eval(changes, angle);
	}

    virtual double get_delta_eval(list< pair< int, double > >& changes, double angle) const=0;

	virtual double get_evaluation()=0;

	const vector< vector<double> >& get_FM(){return FM;}


	void save_sorted_FMs(string filename){
		vector< vector<double> > sortedFM=FM;
		for(int i=0; i<sortedFM.size();i++)
            std::sort(sortedFM[i].begin(), sortedFM[i].end());

		ofstream myfile;
  		myfile.open (filename);
  		myfile << "{ \"fmos\":\n[";
		for(int i=0;i<sortedFM.size();i++){
			myfile << "[";
			for(int j=0;j<sortedFM[i].size();j++){
				myfile << sortedFM[i][j];
				if(j!=sortedFM[i].size()-1) myfile << ",";
			}
			if(i!=sortedFM.size()-1) myfile << "]," << endl;
			else myfile << "]" << endl;
		}
		myfile << "]}\n";
		myfile.close();	
	}

	virtual multimap < double, pair<int, int>, MagnitudeCompare2 > sorted_beamlets(const Plan& p, double vsize)=0;






};

} /* namespace imrt */

#endif /* EVALUATIONFUNCTION_H_ */
