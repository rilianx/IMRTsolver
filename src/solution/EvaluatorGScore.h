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
        int cont=0;


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
                int index = it - sortedFM.begin();
                result.push_back((double)(index+1)/distrSize);
                // linear interpolation?
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
		
	}

    EvaluatorGS(Evaluator& ev) : Evaluator(ev.fm_structure,ev.w,ev.Zmin,ev.Zmax), GS(0.0){
		
	}

	virtual ~EvaluatorGS() { }

    mutable double D10_0, D20_0, D25_0, D30_0, D50_0, Dmax_1, D20_1, D30_1, D50_1, Dmax_2, D80_2, D95_2, D99_2;

    void print(){
        cout << (D10_0/75) << "," << (D20_0/70) << "," << (D25_0/65) << "," << (D30_0/60) << "," << (D50_0/50) <<
        "," << (Dmax_1/77.5) << "," << (D20_1/75) << "," << (D30_1/65) << "," << (D50_1/50) <<
        "," << (76/D95_2) << "," << (69.9/D99_2) << "," << (Dmax_2/85.8);// << "," << (max(76/D80_2-1,0.0));
    }

    double eval(vector< vector<double> >& sortedFM, bool verbose=false) const{
        for(int i=0; i<sortedFM.size();i++)
            std::sort(sortedFM[i].begin(), sortedFM[i].end());


        std::list<double> D2 = computeD(sortedFM[2],vector<double>({0.20, 0.05, 0.01, 1.0}));
        auto D2_it = D2.begin();
        D80_2 = *D2_it; D2_it++; //  >76 ** auxiliar
        D95_2 = *D2_it; D2_it++; //  >76
        D99_2 = *D2_it; D2_it++; //  >69.9
        Dmax_2 = *D2_it;         //  <85.8

        std::list<double> D0 = computeD(sortedFM[0],vector<double>({0.90,0.80,0.75,0.70,0.50}));
        auto D0_it = D0.begin();
        D10_0 = *D0_it; D0_it++; //75
        D20_0 = *D0_it; D0_it++; //70
        D25_0 = *D0_it; D0_it++; //65
        D30_0 = *D0_it; D0_it++; //60
        D50_0 = *D0_it;          //50

        std::list<double> D1 = computeD(sortedFM[1],vector<double>({1.0,0.80,0.70,0.5}));
        auto D1_it = D1.begin();
        Dmax_1 = *D1_it; D1_it++; //77.5
        D20_1 = *D1_it; D1_it++; //75
        D30_1 = *D1_it; D1_it++; //65
        D50_1 = *D1_it;          //50

        if(verbose) cout << D10_0 <<"," << D20_0 <<"," << D25_0 <<"," << D30_0 << ","<< D50_0 << "," << 
        Dmax_1 << "," << D20_1 << "," << D30_1 << "," << D50_1 << "," << Dmax_2 << "," << D95_2 << "," << D99_2 << endl;

        return 0.05*(max(D10_0/75,1.0)) + 0.05*(max(D20_0/70,1.0)) + 0.05*(max(D25_0/65,1.0)) + 0.05*(max(D30_0/60,1.0)) + 0.05*(max(D50_0/50,1.0)) + 
        0.0625*(max(Dmax_1/77.5,1.0)) + 0.0625*(max(D20_1/75,1.0)) + 0.0625*(max(D30_1/65,1.0)) + 0.0625*(max(D50_1/50,1.0)) + 
        0.1666*(max(76/D95_2,1.0)) + 0.1666*(max(69.9/D99_2,1.0)) + 0.1666*(max(Dmax_2/85.8,1.0)) + 0.5*(max(76/D80_2-1,0.0)) ;
        
        //0.125*(max(V65_0/0.25-1,0.0)) + 0.125*(max(V40_1/0.35-1,0.0)) + 0.125*(max(V65_1/0.17-1,0.0)) + 0.5*(max(Zmin[2]/D98_2-1,0.0)) + 0.5*(max(Zmin[2]/D80_2-1,0.0)) );
   

        //if(!gs2)
        //    return (0.125*(V40_0/0.50) + 0.125*(V65_0/0.25) + 0.125*(V40_1/0.35) + 0.125*(V65_1/0.17) + 0.5*(Zmin[2]/D98_2) );
        //else
        //    return (0.125*(max(V40_0/0.50-1,0.0)) + 0.125*(max(V65_0/0.25-1,0.0)) + 0.125*(max(V40_1/0.35-1,0.0)) + 0.125*(max(V65_1/0.17-1,0.0)) + 0.5*(max(Zmin[2]/D98_2-1,0.0)) + 0.5*(max(Zmin[2]/D80_2-1,0.0)) );
    }

    virtual double eval(const Plan& p){
        return eval(p, false);
    }
	// Eval the cost F based on the dose deposition matrix Z
	virtual double eval(const Plan& p, bool verbose){
	    fm_structure.computeFM(p);
        vector< vector<double> > FMM = FM;
        GS =  eval(FMM, verbose);   
        return GS;   
    }

    virtual double incremental_eval(){
        vector< vector<double> > FMM = FM;
        GS=eval(FMM);   
        return GS;           
    }

	virtual double incremental_eval(list< pair< int, double > >& changes, double angle){
        std::vector<list < pair<int, double>>>& deltaFM = fm_structure.compute_deltaFM(changes,angle);
        fm_structure.updateFM(changes, angle, deltaFM);

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

    static bool gs2;

};


} /* namespace imrt */

#endif 
