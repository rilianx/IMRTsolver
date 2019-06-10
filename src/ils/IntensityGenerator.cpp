/*
 * IntensityGenerator.cpp
 *
 *  Created on: Apr 23, 2019
 *      Author: iaraya
 */
#include <algorithm>

#include "IntensityGenerator.h"
#include "Plan.h"
#include "Station.h"


namespace imrt {

void IntensityGenerator::generate(Plan& P, double alpha, int max_intensity){
	Plan Pprime(P);
	cout <<"OldPlan:"<< Pprime.eval() << endl  ;
	Start:
		Pprime.Zsavepoint();
		multimap < double, pair<Station*, int>, MagnitudeCompare > sorted_beamlets;
		Pprime.get_vc_sorted_beamlets(Pprime,sorted_beamlets); //ordeandos v/cx
		for (auto elem:sorted_beamlets){
			Station* s = elem.second.first;
			int b = elem.second.second;
			pair<double,double> vc=Pprime.get_value_cost(s,b);

			int i=s->beam2pos[b].first; int j=s->beam2pos[b].second;
			s->change_intensity(i, j, s->I(i,j) + std::min(int(alpha/vc.second+0.5),max_intensity));

		}

		cout << alpha << endl;
		cout <<"NewPlan:"<< Pprime.eval() << endl  ;
		cout << endl;
  		for(int i=0;i<5;i++)
  		Pprime.printIntensity(i);
  		cout << endl;

		IntensityRepair(Pprime);

}

void IntensityGenerator::IntensityRepair(Plan& P){
	struct potencial_Gap{
		double GAP;
		double Suplier1; //biggest to the rigth
		double Suplier2; //biggest to the left
	}Gaps[8];
	
	list<Station*>& stations = P.get_stations();
	
	for(Station* s:stations){
		for (int i=0;i<10;i++){
			 s->I.get_row(i);

		}
	}

	for(int i;i<=8;i++){ //go over the struct array
		for(int j;j</*row_size */;j++){ // go over the intensity row
			Gaps[i].GAP=j;
			Gaps[i].Suplier1=j+1;
			Gaps[i].Suplier2=j-1;
			for (int k=j+1;k</*row_size */;k++){ //rigth of intensity row
				if(p[k]>p[Gaps[i].Suplier1]){
					Gaps[i].Suplier1=k;
				}
			}
			for(int k=j-1;k>=0;k--){ //left of intensity row
				if(p[k]>p[Gaps[i].Suplier2]){
					Gaps[i].Suplier2=k;
				}
			}
		
		}
	}
	//search for the worst gap
	int RealGap=Gaps[0].GAP;
	for(int i=1;i<=8;i++){
		if(p[RealGap]>p[Gaps[i].GAP]){
			RealGap=Gaps[i].GAP;
		}
	}
}
IntensityGenerator::IntensityGenerator() {
	// TODO Auto-generated constructor stub

}

IntensityGenerator::~IntensityGenerator() {
	// TODO Auto-generated destructor stub
}

} /* namespace imrt */
