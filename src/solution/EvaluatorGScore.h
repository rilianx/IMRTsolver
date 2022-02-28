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

class Score{
    public: 

    enum Type{D,V,H, Dmean};
    Score(Type t, double x, int organ, double min_value, double max_value, double weight) : t(t), x(x), organ(organ), 
    min_value(min_value), max_value(max_value), weight(weight) { }


    Type t;
    double x;
    int organ;
    
    double weight;
    double max_value=0.0;
    double min_value=0.0;

    mutable double value=0.0;
};


class EvaluatorGS : public Evaluator{

public:

    list<Score> scores;

    double computeDmean(vector<double>& sortedFM) const{
        double mean=0.0;
        for (auto d : sortedFM) mean+=d;
        mean/= (double) sortedFM.size();
        return mean;
    }

    double computeD(vector<double>& sortedFM, double x) const{
        x/=100.0;
        int idx = (sortedFM.size()-1)*(1.0-x);
        return sortedFM[idx];
    }

    double computeV(vector<double>& sortedFM, double x) const {
        unsigned distrSize = sortedFM.size();

        if (x <= sortedFM.front()) return (double) 1 / distrSize; // min percentile returned (not important)
        else if (sortedFM.back() <= x) return 1.0; // max percentile returned (not important)
        else
        {
            auto it = lower_bound(sortedFM.begin(), sortedFM.end(), x);
            int index = it - sortedFM.begin();
            return (double)(index+1)/distrSize;
            // linear interpolation?
        }
    }


    double compute_admissible_gscore() const{
        double s=0.0;
        for(const Score& score:scores){
            if(score.t == Score::D || score.t == Score::Dmean){
                if(score.max_value>0.0){
                    if(type==GS) s+=score.weight *  score.value/score.max_value;
                    else if(type==GS_SQUARED) s+=score.weight *  pow(score.value/score.max_value,2);
                    else if (type==GS_RELU) s+=score.weight * (score.value/score.max_value-1.0);
                }else{
                    if(type==GS) s+=score.weight *  score.min_value/score.value;
                    else if(type==GS_SQUARED) s+=score.weight *  pow(score.min_value/score.value,2);
                }
            }
        }
        return s;
    }

    double compute_gscore(vector<vector<double>>& sortedFM) const{
        double s=0.0;

        //identify if the solution is addmisible
        bool admissible = true;

        for(const Score& score:scores){
            if(score.t == Score::D)
                score.value= computeD(sortedFM[score.organ], score.x);
            else if(score.t == Score::Dmean)
                score.value= computeDmean(sortedFM[score.organ]);
            
            //tumor under-irradiated
            if(score.min_value && score.value < score.min_value)
                admissible&=false;
            
            /*if(score.max_value >0){
                if (score.value > score.max_value)
                    admissible&=false;
            }else if (score.value < score.min_value){
                    admissible&=false;
            }*/
        }


        for(const Score& score:scores){
            if(score.t == Score::D || score.t == Score::Dmean){
                if(score.max_value>0.0){ //organs
                    if(type==GS || type==GS2) s+=score.weight *  score.value/score.max_value;
                    else if(type==GS_SQUARED) s+=score.weight *  pow(score.value/score.max_value,2);
                    //else if (type==GS_RELU && !admissible) s+=score.weight *  max(score.value/score.max_value-1.0,0.0);
                    else if (type==GS_RELU && admissible) s+=score.weight * (score.value/score.max_value-1.0);
                }else{ //tumor
                    if(type==GS || (type==GS2 && !admissible) ) s+=score.weight *  score.min_value/score.value;
                    //else if(type==GS && admissible) s+= score.weight;
                    else if(type==GS_SQUARED) s+=score.weight *  pow(score.min_value/score.value,2);
                    else if (type==GS_RELU) s+=score.weight * max(score.min_value/score.value-1.0,0.0);
                }
            }
        }
        if(!admissible) s+=1.0;

        return s;
    }

	//Constructor of the evaluator.
    enum Type{GS, GS2, GS_RELU, GS_SQUARED};
    Type type;

	EvaluatorGS(FluenceMap& fm_structure, vector<double>& w, 
    vector<double>& Zmin, vector<double>& Zmax, list<Score>& scores, Type t) : Evaluator(fm_structure,w,Zmin,Zmax), 
        gs(0.0), scores(scores), type(t) {

	}

    EvaluatorGS(Evaluator& ev, list<Score>& scores, Type t) :   Evaluator(ev.fm_structure,ev.w,ev.Zmin,ev.Zmax), 
    gs(0.0), scores(scores), type(t){
		
	}

	virtual ~EvaluatorGS() { }

    double eval(vector< vector<double> >& sortedFM, bool verbose=false) const{
        for(int i=0; i<sortedFM.size();i++)
            std::sort(sortedFM[i].begin(), sortedFM[i].end());

        return compute_gscore(sortedFM);
   }

    virtual double eval(const Plan& p){
        return eval(p, false);
    }
	// Eval the cost F based on the dose deposition matrix Z
	virtual double eval(const Plan& p, bool verbose){
	    fm_structure.computeFM(p);
        vector< vector<double> > FMM = FM;
        gs =  eval(FMM, verbose);   
        return gs;   
    }

    virtual double incremental_eval(){
        vector< vector<double> > FMM = FM;
        gs=eval(FMM);   
        return gs;           
    }

	virtual double incremental_eval(list< pair< int, double > >& changes, double angle){
        std::vector<list < pair<int, double>>>& deltaFM = fm_structure.compute_deltaFM(changes,angle);
        fm_structure.updateFM(changes, angle, deltaFM);

        vector< vector<double> > FMM = FM;
        gs=eval(FMM);   
        return gs;   
    }

    virtual double get_delta_eval(list< pair< int, double > >& changes, double angle) const { 
        std::vector<list < pair<int, double>>>& deltaFM = fm_structure.compute_deltaFM(changes,angle);

        vector< vector<double> > FMM = FM;
       
        //update FM
        for(int i=0;i<deltaFM.size();i++){
            for(auto p : deltaFM[i])
                FMM[i][p.first]+=p.second;
        }
        return eval(FMM)-gs;        
    }

	virtual double get_evaluation(){ return gs; }



  //Additional functions

    //Return the beamlets sorted by impact on F taking into account the nv worst voxels.
    //Each returned beamlet is a pair eval,(station, beamlet)
    virtual multimap < double, pair<int, int>, MagnitudeCompare2 > sorted_beamlets(const Plan& p, double vsize){
        return multimap < double, pair<int, int>, MagnitudeCompare2 >();
    }

    double gs;

};


} /* namespace imrt */

#endif 
