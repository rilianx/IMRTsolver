/*
 * main.cpp
 *
 *  Created on:
 *      Author: leslie
 */

#include <iostream>
#include <iterator>
#include <set>
#include <stack>
#include <unistd.h>
#include <string.h>

#include "Evaluator.h"
#include "EvaluatorGScore.h"
#include "Plan.h"
#include "IntensityILS2.h"
#include "args.hxx"



#include "solver_extra.cpp"

using namespace std;
using namespace imrt;


int main(int argc, char** argv){

    int seed = time(NULL);

    // Budget and execution variables
    int maxtime = 0;
    int maxeval = 0;

    // Aperture and intensity configuration
    StationSetup initial_setup = StationSetup::open_all_min;
    int max_apertures = 5;
    int max_intensity = 28;
    int step_intensity = 1;

    // Type of local search
    vector<NeighborhoodType> neighborhoods;
    neighborhoods.push_back(NeighborhoodType::sequential_a);
    neighborhoods.push_back(NeighborhoodType::sequential_i);


    // Perturbation
    int perturbation_size = 0;

    // Files
    string path = ".";
    string file = "data/testinstance_0_70_140_210_280.txt";
    string file2 = "data/test_instance_coordinates.txt";
    string output_file = "";

    //Acceptation 
    double min_delta_eval = 0.0001;
    double alpha = 1.0;
    int switch_patience=5;

    //evaluator index
    int sf_eval=0;
    int of_eval=0;


    args::ArgumentParser parser("********* IMRT-Light Solver *********",
                                "Example.\n ./AS -s ibo_ls --setup=open_min --ls_sequential=aperture -s ibo_ls --maxeval=15000 --ls=first --perturbation-size=5 --seed=1 --max-intensity=20 --file-coord=data/Equidistantes/equidist-coord.txt --initial-intensity=5");

    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<int>    _seed     (parser, "int", "Seed  (" + 
                                        to_string(seed)+")", {"seed"});

    // Execution parameters
    args::Group budget (parser, "Budget options:");
    args::ValueFlag<int>    _maxeval  (budget, "int",
                                    "Number of evaluations (" +
                                        to_string(maxeval)+ ")", {"maxeval"});


    // Initial collimator setup (initial solution: aperture, intensity)
    args::Group isetup (parser, "Initial collimator setup:");
    args::ValueFlag<int> _max_apertures     (isetup, "int",
                                            "Number of apertures per angle (station) (" +
                                            to_string(max_apertures)+")", {"max-apertures"});

    // Neighborhood parameters
    args::Group neighborhoodsel (parser, "Neighborhood selection:");     
    args::ValueFlag<string> _neighborhoods (neighborhoodsel , "string",
                                "neighborhoods in local search",
                                {"neighborhoods"});


    args::Group accargs (parser, "Acceptation improvement:");
    args::ValueFlag<double> _min_delta_eval (accargs, "int",
                                            "Minimum delta eval for accepting the change",
                                            {"min-delta"});
    args::ValueFlag<double> _alpha (accargs, "double",
                                            "Reduction factor of min-delta after each evaluation",
                                            {"min-delta-red"});
    args::ValueFlag<string> _pr_neigh (accargs, "string",
                                            "Prop. of elements of each neighbourhood",
                                            {"pr-neigh"});


    // Perturbation parameters
    args::Group perargs (parser, "Perturbation:");
    args::ValueFlag<int> _perturbation_size (perargs, "int",
                                            "Perturbation size  (" +
                                            to_string(perturbation_size)+")",
                                            {"perturbation-size"});

    // Objective function
    args::Group objfunct (parser, "Evaluators:");
    args::ValueFlag<string> _evaluators (objfunct, "string",
                                    "Files with function+scores. ",
                                    {"evals"});
    args::ValueFlag<int> _sf_eval (objfunct, "int",
                                    "index of the search function",
                                    {"sf"});
    args::ValueFlag<int> _of_eval (objfunct, "int",
                                    "index of the objective function",
                                    {"of"});


    args::ValueFlag<int> _switch_patience (objfunct, "int",
                                            "Number of movements that does not improve OF  before switching the search function to OF.",
                                            {"switch-patience"});

    


    // Problem file parameters
    args::Group io_opt (parser, "Input output options:");
    args::ValueFlag<string> _file  (io_opt, "string",
                                    "File with the deposition matrix", {"file-dep"});
    args::ValueFlag<string> _file2 (io_opt, "string",
                                    "File with the beam coordinates", {"file-coord"});
    args::ValueFlag<string> _path  (io_opt, "string",
                                    string("Absolute path of the executable ") +
                                    "(if it is executed from other directory)", {"path"});
    args::Flag _verbose               (io_opt, "bool",
                                    "Verbose", {"verbose"});

    // Output file parameters
    args::ValueFlag<string> _output_file (io_opt, "string",
                                    "File to output all indicators for each iteration", 
                                    {"output-file"});
    args::ValueFlag<string> _fm_output (io_opt, "string",
                                    "File to output the fluence map of voxels", 
                                    {"output-fm"});

                                
	try
	{
		parser.ParseCLI(argc, argv);

	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

    if(_maxeval) maxeval = _maxeval.Get();
    if(_seed) seed = _seed.Get();
    if(_max_apertures) max_apertures=_max_apertures.Get();

    //Acceptation improvement
    if(_min_delta_eval) min_delta_eval=_min_delta_eval.Get();
    if(_alpha) alpha=_alpha.Get();

    // Neighborhood
    if(_neighborhoods){
        neighborhoods.clear();
        vector<std::string> neigh = split(_neighborhoods.Get(), ',');
        for(auto nn : neigh){
            if (nn == "intensity") 
                neighborhoods.push_back(NeighborhoodType::intensity);
            else if (nn == "aperture") 
                neighborhoods.push_back(NeighborhoodType::aperture);
            
        }
    }


    // Neighborhood pr
    vector<double> pr_neigh(neighborhoods.size(),1.0);
    if(_pr_neigh){
        pr_neigh.clear();
        vector<std::string> neigh = split(_pr_neigh.Get(), ',');
        for(auto nn : neigh){
            pr_neigh.push_back(stod(nn));
        }
    }

  
    if (_perturbation_size) perturbation_size = _perturbation_size.Get();

    // Archivos del problema
    if (_file) file=_file.Get();
    if (_file2) file2=_file2.Get();
    if (_path) path=_path.Get();

    int aux=chdir(path.c_str());
    string mkdir = "mkdir output";
    system(mkdir.c_str());

    // Iniciar generador de numeros aleatorios
    srand(seed);

    // Create colimator object and volumes
    Collimator collimator (file2);  //, get_angles(file, 5));
    vector<Volume> volumes = createVolumes (file, collimator);
    vector<Evaluator*> evaluators = createEvaluators(collimator, volumes, _evaluators.Get());

    if(_sf_eval) sf_eval=_sf_eval.Get();
    if(_of_eval) of_eval=_of_eval.Get();
    
    if(_switch_patience) switch_patience=_switch_patience.Get();

    //output file

    ofstream output_stream;
    if(_output_file) {
        output_file=_output_file.Get();
        output_stream.open (output_file.c_str(), ios::out);
    }

    cout << endl << "##*********************************** INFO *********************************" << endl;

    cout << "##" << endl << "## Solver: "<< endl;
    cout << "##   Evaluations: " << maxeval << endl;
    cout << "##   Seed: " << seed << endl;

    cout << "##   Apertures: " << max_apertures << endl;
    cout << "##   Max intensity: " << max_intensity << endl;

    cout << "##   Neighborhoods:";
    for(auto neighborhood:neighborhoods){
        if (neighborhood == NeighborhoodType::intensity)
            cout << "intensity - ";
        else if (neighborhood == NeighborhoodType::aperture)
            cout << "aperture - ";
    }

    cout << endl << "##   Perturbation size: " << perturbation_size << endl;

    cout << "##" << endl << "## Colimator configuration: "<< endl;
    cout << "##   Stations: " << collimator.getNbAngles() << endl;
    cout << "##   Angles: ";

    for (int i=0; i<collimator.getNbAngles();i++)
        cout << collimator.getAngle(i) << " ";
    cout << endl;


    cout << "##" << endl << "## Instance information: "<< endl;
    cout << "##   Volumes: " << volumes.size() << endl;

    cout << "##" << endl << "## Output files: " << endl;
    if (_output_file)
    cout << "##   Output file: " << output_file << endl;

    cout << "##********************************** SEARCH ********************************" << endl;

    ifstream _ang_file(file.c_str(), ios::in);
    string angles_str;
    getline(_ang_file, angles_str);
    _ang_file.close();
    list<int> bac = get_angles(angles_str);

    cout << "pr_neigh: " ;
    for(double b:pr_neigh) cout << b << " ";
    cout << endl;


    std::clock_t begin_time = clock();
    int used_evaluations=0;


    double cost = iterated_local_search(collimator, volumes, 
            max_apertures, max_intensity, step_intensity, 
            initial_setup, bac, 
            evaluators, sf_eval, of_eval, 
            maxeval, used_evaluations, neighborhoods, perturbation_size,
            output_stream, begin_time, min_delta_eval, alpha,
            switch_patience, pr_neigh, _verbose );

    

    
    if (output_stream.is_open()){
        for(auto ev:evaluators) ev->incremental_eval() ;
        clock_t time_end = clock();
        double used_time = double(time_end - begin_time) / CLOCKS_PER_SEC;
        output_stream <<  used_evaluations <<";";
        output_stream <<  used_time <<";";
        for (auto score : dynamic_cast<EvaluatorGS*>(evaluators[of_eval])->scores )
            output_stream << score.value << ";";
        for(auto ev:evaluators)
            output_stream << ev->get_evaluation() << ";";
        output_stream << endl;
        output_stream.close();
    }

    ofstream output_fm_stream;
    if(_fm_output) {
        string file=_fm_output.Get();
        output_fm_stream.open (file.c_str(), ios::out);
    }

    for(int o=0; o< evaluators[0]->FM.size(); o++){
        for(int i=0; i<evaluators[0]->FM[o].size(); i++){
            output_fm_stream << evaluators[0]->FM[o][i] << ",";
        }
        output_fm_stream << endl;
    }



	return 0;

}
