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
    LSType ls_type = LSType::first;
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
    double pr_first_neigh=1.0;

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
    args::ValueFlag<double> _pr_first_neigh (accargs, "double",
                                            "Prob. of selecting moves in the first neighbourhood",
                                            {"pr-first-neigh"});


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
    if(_pr_first_neigh) pr_first_neigh=_pr_first_neigh.Get();

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
    vector<double> w={1,1,5};
    vector<double> Zmin={0,0,76};
    vector<double> Zmax={65,65,76};
    Collimator collimator (file2, get_angles(file, 5));
    vector<Volume> volumes = createVolumes (file, collimator);

    FluenceMap fm(volumes, collimator);

    vector<Evaluator*> evaluators;
    
    
    if(_evaluators){
        vector<std::string> files = split(_evaluators.Get(), ',');

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
            
            evaluators.push_back(new EvaluatorGS(fm,w,Zmin,Zmax, scores, t));
        }
    }

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
        if (neighborhood == NeighborhoodType::sequential_i)
            cout << "intensity - ";
        else if (neighborhood == NeighborhoodType::sequential_a)
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

    // Create an initial plan
    Plan P (collimator, volumes, max_apertures,
            max_intensity, 0, step_intensity,
            initial_setup);

    double best_eval = evaluators[sf_eval]->eval(P);
    cout << "## Initial solution: " << best_eval << endl;
 
    //IBO_LS
    ILS* ils = new IntensityILS2(evaluators, sf_eval, of_eval);
    std::clock_t begin_time = clock();
    int used_evaluations=0;
    double cost = ils->iteratedLocalSearch(P, maxeval, neighborhoods, perturbation_size, 
                    output_stream,used_evaluations, begin_time, _verbose, min_delta_eval, alpha, switch_patience, pr_first_neigh);
    
    
  
    EvaluatorF ev(fm,w,Zmin,Zmax);
    cout << "## Best solution found: " << ils->best_evals[of_eval] << endl;
    

    if (output_stream.is_open()){
        evaluators[0]->eval(P);
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




	return 0;

}
