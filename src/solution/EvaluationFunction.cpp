/*
 * EvaluationFunction.cpp
 *
 *  Created on: 8 may. 2018
 *      Author: iaraya
 */

#include "EvaluationFunction.h"
#include <string>
#include <functional>

namespace imrt {

int EvaluationFunction::n_evaluations=0;
EvaluationFunction* EvaluationFunction::instance=NULL;

void EvaluationFunction::create_voxel2beamlet_list(vector<Volume>& volumes, const Collimator& collimator){
	for (int i=0;i<collimator.getNbAngles();i++){
		for(int o=0; o<nb_organs; o++){
			 int angle = collimator.getAngle(i);
			 const Matrix&  D = volumes[o].getDepositionMatrix(angle);
			 for(int k=0; k<nb_voxels[o]; k++){
				 for(int b=0; b<collimator.getNangleBeamlets(angle); b++)
					 if(D(k,b) > 0.0) {
						 voxel2beamlet_list[make_pair(o,k)].push_back(make_pair(angle,b));
						 beamlet2voxel_list[make_pair(angle,b)].insert(make_pair(-abs(D(k,b)), make_pair(o,k)));
					 }
			 }
		}
	}
}

EvaluationFunction::EvaluationFunction(vector<Volume>& volumes, const Collimator& collimator) : prev_F(0.0), F(0.0),
	volumes(volumes), nb_organs(volumes.size()), nb_voxels(volumes.size()), voxel_dose(volumes.size(), vector<double>(150)) {
  n_volumes =volumes.size();
	for(int i=0; i<nb_organs; i++){
		nb_voxels[i]=volumes[i].getNbVoxels();
		Z.insert(Z.end(), vector<double>(nb_voxels[i]));
		D.insert(D.end(), vector<double>(nb_voxels[i]));
	}

	create_voxel2beamlet_list(volumes, collimator);

}

EvaluationFunction::~EvaluationFunction() { }

void EvaluationFunction::generate_linear_system(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax){
	bool flag=false;
	for(int o=0; o<nb_organs; o++){
		double pen=0.0;
		for(int k=0; k<nb_voxels[o]; k++){
			if(flag) cout << " + ";
			flag = true;
			if(k==0) cout << "(";
			if(Zmin[o] == 0 )
				cout << "max(Z_" << o << "_" << k << "-" << Zmax[o] << ", 0)^2" ;
			else
				cout << "max(" << Zmin[o] << "- Z_" << o << "_" << k << ", 0)^2" ;
			if(k==nb_voxels[o]-1) cout << ")";
		}
		cout << "/" << nb_voxels[o];
	}
	cout << endl;

	for(int o=0; o<nb_organs; o++){
		for(int k=0; k<nb_voxels[o]; k++){
			cout << "Z_" << o << "_" << k << "=" ;
			int ap=0; //aperture variable id
			bool flag = false;
			for(auto s : p.get_stations()){
				const Matrix&  D = s->getDepositionMatrix(o);
				for(int a=0; a<s->getNbApertures(); a++){
					double ap_dose_per_intensity=0.0;
					for(auto b: s->open_beamlets(a)) ap_dose_per_intensity+=D(k,b);

					if (ap_dose_per_intensity>0.0){
						if(flag) cout << " + ";
						cout << ap_dose_per_intensity << " * I_" << ap;
						flag = true;

					}
					ap++;
				}
			}
			cout << endl;
		}
	}
}


void EvaluationFunction::generate_Z(const Plan& p){
	const list<Station*>& stations=p.get_stations();

	for(int o=0; o<nb_organs; o++)
	 	 std::fill(Z[o].begin(), Z[o].end(), 0.0);

	//we update the dose distribution matrices Z with the dose delivered by the station
	for(int o=0; o<nb_organs; o++){
		 for(int k=0; k<nb_voxels[o]; k++){
		   double dose=0.0;
			 list<pair <int,int> > beamlets = voxel2beamlet_list[make_pair(o,k)];
			 for(auto bl:beamlets){
				 int angle=bl.first; int b=bl.second;
				 const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
				 dose += Dep(k,b)*p.angle2station.at(angle)->getIntensity(b);
			 }
		   Z[o][k] += dose;
		 }
	}
}


double EvaluationFunction::eval(const Plan& p, vector<double>& w, vector<double>& Zmin, vector<double>& Zmax){
	voxels.clear();
	generate_Z(p);
	F=0.0;
	for(int o=0; o<nb_organs; o++){
		double pen=0.0;
		for(int k=0; k<nb_voxels[o]; k++){
			if(Z[o][k] < Zmin[o] )
				 pen += w[o] * ( pow(Zmin[o]-Z[o][k], 2) );

			if(Z[o][k] > Zmax[o] )
				 pen += w[o] * ( pow(Z[o][k]-Zmax[o], 2) );

			update_sorted_voxels(w, Zmin, Zmax, o, k);
			update_beamlets_impact(o, k);
		}
		//cout << pen/nb_voxels[o] << endl;
		F+=pen/nb_voxels[o];
	}

	//for(auto b:beamlet_impact)
		//cout << b.first.first << "," << b.first.second << ":" << b.second << " " << endl;

	n_evaluations++;
	return F;
}

//Update D (partial derivative of F w.r.t. the voxels)
//And resorts the set of voxels
void EvaluationFunction::update_sorted_voxels(vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax, int o, int k){
		voxels.erase(make_pair(D[o][k],make_pair(o,k)));
    D[o][k] = 0.0;

		if(Zmin[o]>0 && Z[o][k] < Zmin[o]) D[o][k]=-w[o]*(Z[o][k]-Zmin[o])/nb_voxels[o];

		else if(Z[o][k] >= Zmax[o]) D[o][k]=w[o]*(Z[o][k]-Zmax[o])/nb_voxels[o];

		if(D[o][k]!=0.0)
			voxels.insert(make_pair(D[o][k],make_pair(o,k)));
}

void EvaluationFunction::update_beamlets_impact(int o, int k, double prev_Dok){
	list<pair <int,int> >& beamlets = voxel2beamlet_list[make_pair(o,k)];
	for(auto bl:beamlets){
		beamlet_impact[bl] -= prev_Dok;
		int angle = bl.first;
		int b = bl.second;
		const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
		beamlet_impact[bl] += D[o][k] * Dep(k,b);
	}


}

double EvaluationFunction::get_ratio_beamlet(vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax, int angle, int b){
			double ev=0;
			multimap<double, pair<int,int> > voxels = beamlet2voxel_list[make_pair(angle,b)];
			double T=0.0,O=0.0;
			for(auto voxel:voxels){
				int o=voxel.second.first;
				int k=voxel.second.second;
				const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);

				if(Zmin[o]>0){
					 T += max(Dep(k,b) * w[o]*(1.0 - Z[o][k]/Zmin[o])/nb_voxels[o],0.0);
				}else
					 O += Dep(k,b) * w[o]*(1.0 + Z[o][k]/Zmax[o])/nb_voxels[o];

			}
			//cout << T << "," << O << endl;
			return (T/O);
}

double EvaluationFunction::get_impact_beamlet(int angle, int b){
			double ev=0;
			multimap<double, pair<int,int> > voxels = beamlet2voxel_list[make_pair(angle,b)];
			for(auto voxel:voxels){
				int o=voxel.second.first;
				int k=voxel.second.second;
				const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
				ev += D[o][k] * Dep(k,b);
			}
			return ev;
}

multimap < double, pair<Station*, int>, MagnitudeCompare >
EvaluationFunction::get_sorted_beamlets(Plan& p){
	//the sorted set of beamlets
	multimap < double, pair<Station*, int>, MagnitudeCompare > sorted_set;

	double max_ev=0.0;
	for(auto s:p.get_stations()){
		for(int b=0; b<s->getNbBeamlets(); b++){
			double ev=0; int i=0;
			for(auto voxel:voxels){
				int o=voxel.second.first;
				int k=voxel.second.second;
				const Matrix&  Dep = s->getDepositionMatrix(o);
				ev += D[o][k] * Dep(k,b);
			}
			if(ev==0) continue;
			sorted_set. insert(make_pair (ev,  make_pair(s,b)));
		}
	}

	return sorted_set;
}

double EvaluationFunction::get_delta_eval(int angle, int b, double delta_intensity,
	vector<double>& w, vector<double>& Zmin, vector<double>& Zmax,int n_voxels) const{

	 multimap<double, pair<int,int> > voxels =
	                                   beamlet2voxel_list.at(make_pair(angle,b));
	 double delta_F=0.0;

	 int i=0;
	 for(auto voxel:voxels){
		 int o=voxel.second.first, k=voxel.second.second;
		 const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
		 double delta = Dep(k,b)*(delta_intensity);
		 if(delta==0.0) continue;

		 double pen=0.0;
 		 //with the change in the dose of a voxel we can incrementally modify the value of F
 		 if(Z[o][k] < Zmin[o] && Z[o][k] + delta < Zmin[o]) //update the penalty
 			pen += w[o]*delta*(delta+2*(Z[o][k]-Zmin[o]));
 		 else if(Z[o][k] < Zmin[o]) //the penalty disappears
 			pen -=  w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
 		 else if(Z[o][k] + delta < Zmin[o]) //the penalty appears
 			pen +=  w[o] * ( pow(Zmin[o]-(Z[o][k]+delta), 2) );

 		 if(Z[o][k] > Zmax[o] && Z[o][k] + delta > Zmax[o]) //update the penalty
 			 pen += w[o]*delta*(delta+2*(-Zmax[o] + Z[o][k]));
 		 else if(Z[o][k] > Zmax[o]) //the penalty disappears
 			 pen -=  w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
 		 else if(Z[o][k] + delta > Zmax[o]) //the penalty appears
 			 pen +=  w[o] * ( pow(Z[o][k]+delta - Zmax[o], 2) );

 		 delta_F += pen/nb_voxels[o];
		 i++; if(i>n_voxels) break;
	 }
	 return delta_F;
}

double EvaluationFunction::get_delta_eval(list< pair< int, double > >& diff, double angle,
                                          vector<double>& w, vector<double>& Zmin,
                                          vector<double>& Zmax,int n_voxels) const{

  double delta_F=0.0;
  int b;
  double delta_intensity, delta;

  int i=0;
  for(int o=0; o<nb_organs; o++){
		for(int k=0; k<nb_voxels[o]; k++){
  //for(auto voxel:voxels){
    //int o=voxel.second.first;
    //int k=voxel.second.second;
    const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);

    delta = 0;
	for (auto let:diff){
	    int b=let.first;
		if(Dep(k,b)==0.0) continue;
	    delta+= Dep(k,b)*let.second;
	}
    if(delta==0.0) continue;

	double pen=0.0;
	//with the change in the dose of a voxel we can incrementally modify the value of F
	if(Z[o][k] < Zmin[o] && Z[o][k] + delta < Zmin[o]) //update the penalty
		pen += w[o]*delta*(delta+2*(Z[o][k]-Zmin[o]));
	else if(Z[o][k] < Zmin[o]) //the penalty disappears
		pen -=  w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
	else if(Z[o][k] + delta < Zmin[o]) //the penalty appears
		pen +=  w[o] * ( pow(Zmin[o]-(Z[o][k]+delta), 2) );

	if(Z[o][k] > Zmax[o] && Z[o][k] + delta > Zmax[o]) //update the penalty
		pen += w[o]*delta*(delta+2*(-Zmax[o] + Z[o][k]));
	else if(Z[o][k] > Zmax[o]) //the penalty disappears
		pen -=  w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
	else if(Z[o][k] + delta > Zmax[o]) //the penalty appears
		pen +=  w[o] * ( pow(Z[o][k]+delta - Zmax[o], 2) );

    delta_F += pen/nb_voxels[o];
    i++; if(i>n_voxels){break;}
  }
  }
  return delta_F;
}


double EvaluationFunction::incremental_eval(Station& station, vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax, list< pair< int, double > >& diff){

  prev_F=F; Z_diff.clear();
  double delta_F=0.0;

  //for each voxel we compute the change produced by the modified beamlets
  //while at the same time we compute the variation in the function F produced by all these changes

  for(int o=0; o<nb_organs; o++){
		const Matrix&  Dep = station.getDepositionMatrix(o);
		for(int k=0; k<nb_voxels[o]; k++){

		//we compute the change in the delivered dose in voxel k of the organ o
		double delta=0.0;

		//cout << station.changed_lets.size() << endl;
		for (auto let:diff){
		    int b=let.first;
			if(Dep(k,b)==0.0) continue;
		    delta+= Dep(k,b)*let.second;
		}
		if(delta==0.0) continue; //no change in the voxel


		double pen=0.0;
		//with the change in the dose of a voxel we can incrementally modify the value of F
		if(Z[o][k] < Zmin[o] && Z[o][k] + delta < Zmin[o]) //update the penalty
			pen += w[o]*delta*(delta+2*(Z[o][k]-Zmin[o]));
		else if(Z[o][k] < Zmin[o]) //the penalty disappears
			pen -=  w[o] * ( pow(Zmin[o]-Z[o][k], 2) );
		else if(Z[o][k] + delta < Zmin[o]) //the penalty appears
			pen +=  w[o] * ( pow(Zmin[o]-(Z[o][k]+delta), 2) );

		if(Z[o][k] > Zmax[o] && Z[o][k] + delta > Zmax[o]) //update the penalty
			pen += w[o]*delta*(delta+2*(-Zmax[o] + Z[o][k]));
		else if(Z[o][k] > Zmax[o]) //the penalty disappears
			pen -=  w[o] * ( pow(Z[o][k]-Zmax[o], 2) );
		else if(Z[o][k] + delta > Zmax[o]) //the penalty appears
			pen +=  w[o] * ( pow(Z[o][k]+delta - Zmax[o], 2) );

		delta_F += pen/nb_voxels[o];
		Z[o][k] += delta;

		double prev_Dok=D[o][k];
		update_sorted_voxels(w, Zmin, Zmax, o, k);
		//update_beamlets_impact(o, k, prev_Dok);

		//we save the last changes (see undo_last_eval)
		Z_diff.push_back(make_pair(make_pair(o,k),delta));
	}
  }

  F+=delta_F;
  n_evaluations++;
  return F;

}

void EvaluationFunction::undo_last_eval(vector<double>& w,
	vector<double>& Zmin, vector<double>& Zmax){

	for(auto z:Z_diff){
		int o=z.first.first;
		int k=z.first.second;
		Z[o][k]-=z.second;

		double prev_Dok=D[o][k];
		update_sorted_voxels(w, Zmin, Zmax, o, k);
		//update_beamlets_impact(o, k, prev_Dok);
	}
	F=prev_F;
	prev_F=F; Z_diff.clear();
}



multimap < double, pair<int, int>, MagnitudeCompare >
EvaluationFunction::best_beamlets(Plan& p, int nv){
	multimap < double, pair<int, int>, MagnitudeCompare > bestb;

	double max_ev=0.0;
	int k=0;
	for(auto s:p.get_stations()){
		for(int b=0; b<s->getNbBeamlets(); b++){
			double ev=0; int i=0;

      for (auto voxel : voxels){
				  double d = voxel.first;
					int o = voxel.second.first;
					int k= voxel.second.second;

					const Matrix&  Dep = s->getDepositionMatrix(o);
					ev += d * Dep(k,b);
					i++; if(i==nv) break;
			 }

				bestb.insert(make_pair( ev, make_pair(k,b)));

			}
		k++;
	}

	return bestb;
}

void EvaluationFunction::generate_voxel_dose_functions (){

	for(int o=0; o<nb_organs; o++){
		std::set<double, std::greater<double> > dose;
		std::fill(voxel_dose[o].begin(), voxel_dose[o].end(), 0.0);
		for(int k=0; k<nb_voxels[o]; k++){

			dose.insert(Z[o][k]);
			if(Z[o][k]<150)
				voxel_dose[o][(int) Z[o][k]]+=1;
		}

		ofstream myfile;
		myfile.open ("plotter/organ_xls"+std::to_string(o)+".dat");
		for(auto v:dose){
			myfile << v << endl;
		}
		myfile.close();

	}


	for(int o=0; o<nb_organs; o++){
		ofstream myfile;
		myfile.open ("plotter/organ"+std::to_string(o)+".txt");

		double cum=0.0;
		for(int k=149; k>=0; k--){
			cum+= voxel_dose[o][k];
			myfile << k+1 << "," << cum/nb_voxels[o] << endl;
		}

		myfile.close();
	}
}

bool EvaluationFunction::Zupdate(int angle, int b, double delta_intensity, bool return_if_unfeasible,
			vector<double>& Zmax){

	multimap<double, pair<int,int> > voxels = beamlet2voxel_list.at(make_pair(angle,b));
	bool feasible = true;
	 for(auto voxel:voxels){
		 int o=voxel.second.first, k=voxel.second.second;
		 const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
		 double delta = Dep(k,b)*(delta_intensity);
		 if(delta==0.0) continue;

		 if(Z[o][k] + delta > Zmax[o]){
			 feasible=false;
			 if(return_if_unfeasible) return feasible;
		 }

		 Z[o][k] += delta;
		 Z_diff.push_back(make_pair(make_pair(o,k),delta));
	 }

	 return feasible;
}

void EvaluationFunction::Zrollback(){
	for(auto z:Z_diff){
		int o=z.first.first;
		int k=z.first.second;
		Z[o][k]-=z.second;
	}
	Z_diff.clear();
}

pair<double,double> EvaluationFunction::get_value_cost(int angle, int b, vector<double>& Zmin, vector<double>& Zmax){
	multimap<double, pair<int,int> > voxels = beamlet2voxel_list.at(make_pair(angle,b));
	double cost=0.0, value=0.0;
	 for(auto voxel:voxels){
		 int o=voxel.second.first, k=voxel.second.second;
		 const Matrix&  Dep = volumes[o].getDepositionMatrix(angle);
		 if(Zmin[o]==0.0){//organ
			 if(Zmax[o]-Z[o][k] == 0) cost=1e20;
			 else{
				 double c=Dep(k,b)/(Zmax[o]-Z[o][k]);
			 	 if(c>cost) cost=c;
			 }
		 }else
			 if(Z[o][k] < Zmin[o]) value += Dep(k,b);
	 }

	 return make_pair(value,cost);
}

void EvaluationFunction::get_vc_sorted_beamlets(Plan& p, vector<double>& Zmin, vector<double>& Zmax,
		multimap < double, pair<Station*, int>, MagnitudeCompare >& sorted_set){
	//the sorted set of beamlets

	double max_ev=0.0;
	for(auto s:p.get_stations()){
		for(int b=0; b<s->getNbBeamlets(); b++){
			pair<double,double> value_cost=get_value_cost(s->getAngle(), b, Zmin, Zmax);
            double ev=value_cost.first/value_cost.second;
			if(ev==0) continue;
			sorted_set.insert(make_pair (ev,  make_pair(s,b)));
		}
	}

}

} /* namespace imrt */
