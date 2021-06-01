/*
 * EvaluatorF.h
 *
 *  Created on: 1 june 2021
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
#include "Evaluator.h"

#ifndef EvaluatorGS_H_
#define EvaluatorGS_H_

using namespace std;
using namespace maths;


namespace imrt {

class Plan;


/* The Evaluator GS
 */


class EvaluatorGS : public Evaluator{

public:

    list<double> computeV(vector<double>& sortedFM, const vector<double>& V) const {
        double prevValue, nextValue;
        list<double> result;
        unsigned distrSize = sortedFM.size();


        for (vector<double>::const_iterator test = V.begin(); test != V.end(); test++)
        {
            
            if (*test <= sortedFM.front())
            {
                result.push_back((double) 1 / distrSize); // min percentile returned (not important)
            }
            else if (sortedFM.back() <= *test)
            {
                result.push_back(1); // max percentile returned (not important)
            }
            else
            {
                auto it = lower_bound(sortedFM.begin(), sortedFM.end(), *test);
                prevValue = *(it - 1);
                nextValue = *(it + 1);
                // linear interpolation
                result.push_back(((*test - prevValue) / (nextValue - prevValue) + (it - V.begin())) / distrSize);
            }
        }
        return result;
    }

    list<double> computeD(vector<double>& sortedFM, const vector<double>& D ) const{
        double prevValue, nextValue;
        list<double> results;

        for (int i=0; i<D.size(); i++)
        {
            int idx = (sortedFM.size()-1)*D[i];
            results.push_back(sortedFM[idx]);
        }
        return results;
    }

	//Constructor of the evaluator.
	EvaluatorGS(FluenceMap& fm_structure, vector<double>& w, 
    vector<double>& Zmin, vector<double>& Zmax) : Evaluator(fm_structure,w,Zmin,Zmax), GS(0.0){
		
	};


	virtual ~EvaluatorGS() { }

    double eval(vector< vector<double> >& sortedFM) const{
        for(int i=0; i<sortedFM.size();i++)
            std::sort(sortedFM[i].begin(), sortedFM[i].end());


        double D50_0 = computeD(sortedFM[0],vector<double>({0.5})).front();
        double D50_1 = computeD(sortedFM[1],vector<double>({0.5})).front();
        double D05_2 = computeD(sortedFM[2],vector<double>({0.05})).front();
        double V65_0 = computeV(sortedFM[0],vector<double>({65})).front();
        double V65_1 = computeV(sortedFM[1],vector<double>({65})).front();

        return (D50_0/Zmax[0] + D50_1/Zmax[1] + Zmin[2]/D05_2);
    }

	// Eval the cost F based on the dose deposition matrix Z
	virtual double eval(const Plan& p){
	    fm_structure.computeFM(p);
        vector< vector<double> > FMM = FM;
        GS =  eval(FMM);   
        return GS;   
    }

	virtual double incremental_eval(list< pair< int, double > >& changes, double angle){
        vector< vector<double> > FMM = FM;
        GS=eval(FMM);   
        return GS;   
    }

    virtual double get_delta_eval(list< pair< int, double > >& changes, double angle) const { 
        std::vector<list < pair<int, double>>>& deltaFM = fm_structure.compute_deltaFM(changes,angle);

        vector< vector<double> > FMM = FM;
       
        //update FM
        for(int i=0;i<deltaFM.size();i++){
            for(auto p : deltaFM[i])
                FMM[i][p.first]+=p.second;
        }
        return eval(FMM)-GS;        
    }

	virtual double get_evaluation(){ return GS; }



  //Additional functions

    //Return the beamlets sorted by impact on F taking into account the nv worst voxels.
    //Each returned beamlet is a pair eval,(station, beamlet)
    virtual multimap < double, pair<int, int>, MagnitudeCompare2 > sorted_beamlets(const Plan& p, double vsize){
        return multimap < double, pair<int, int>, MagnitudeCompare2 >();
    }

    double GS;

};

} /* namespace imrt */

#endif 
