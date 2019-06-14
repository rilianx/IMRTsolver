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
	list<Gap> gaps;
	Gap gap;
	list<Station*>& stations=P.get_stations();
	double* intensities;
	for(Station* s:stations){
		for (int i=0;i<s->I.nb_rows();i++){
			intensities=s->I.get_row(i);
			for(int j=1;j<s->I.nb_cols();j++){ // go over the intensity row
				gap.GAP=j;
				gap.Suplier1=j+1;
				gap.Suplier2=j-1;
				for (int k=j+1;k<s->I.nb_cols();k++){ //rigth of intensity row
					if(intensities[k]>intensities[gap.Suplier1]){
						gap.Suplier1=k;
					}
				}
				for(int k=j-1;k>=0;k--){ //left of intensity row
					if(intensities[k]>intensities[gap.Suplier2]){
						gap.Suplier2=k;
					}
				}
				if(intensities[gap.GAP] < intensities[gap.Suplier1] && intensities[gap.GAP] < intensities[gap.Suplier2]){
					gaps.push_back(gap);	
				}
	
				if(gaps.empty()) break;
				changeworst(intensities, gaps);
			}
		}
	}
}
//search for the worst gap and changeit
void IntensityGenerator::changeworst(double* intensities,list<Gap> gaps){
	Gap RealGap=gaps.front();
	for(Gap gap: gaps){
		if(intensities[RealGap.GAP] > intensities[gap.GAP]){
			RealGap=gap;
		}
	}
	//cambiar intensidades RealGap
	if(intensities[RealGap.Suplier1]<intensities[RealGap.Suplier2]){
		intensities[RealGap.GAP]=intensities[RealGap.GAP] - 1;
		intensities[RealGap.Suplier2]=intensities[RealGap.Suplier2] - 1;
	} else {
		intensities[RealGap.GAP]=intensities[RealGap.GAP] - 1;
		intensities[RealGap.Suplier1]=intensities[RealGap.Suplier1] - 1;
	}

}
IntensityGenerator::IntensityGenerator(){
	// TODO Auto-generated constructor stub
}

IntensityGenerator::~IntensityGenerator() {
	// TODO Auto-generated destructor stub
}

} /* namespace imrt */
