/*
 * main.cpp
 *
 *  Created on:
 *      Author: iaraya, leslie
 */


/*******Sockets stuff*********/
#include <unistd.h>
#include <stdio.h> 
#include <sys/socket.h> 
#include <stdlib.h> 
#include <netinet/in.h> 
/******************/

#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <stack>
#include <unistd.h>
#include <string.h>

#include "IntensityGenerator.h"
#include "EvaluationFunction.h"
#include "Plan.h"
#include "Collimator.h"
#include "Volume.h"
#include "ApertureILS.h"
#include "IntensityILS2.h"
#include "MixedILS.h"
#include "args.hxx"
#include <dirent.h>


using namespace std;
using namespace imrt;


void initialize_socket(int& server_fd, struct sockaddr_in& address, int port){
    int opt = 1; 
    int addrlen = sizeof(address); 

       
    // Creating socket file descriptor 
    if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) == 0) 
    { 
        perror("socket failed"); 
        exit(EXIT_FAILURE); 
    } 
       
    // Forcefully attaching socket to the port 8080 
    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, 
                                                  &opt, sizeof(opt))) 
    { 
        perror("setsockopt"); 
        exit(EXIT_FAILURE); 
    } 
    address.sin_family = AF_INET; 
    address.sin_addr.s_addr = INADDR_ANY; 
    address.sin_port = htons( port ); 
       
    // Forcefully attaching socket to the port 8080 
    if (bind(server_fd, (struct sockaddr *)&address,  
                                 sizeof(address))<0) 
    { 
        perror("bind failed"); 
        exit(EXIT_FAILURE); 
    } 

    if (listen(server_fd, 3) < 0) 
    { 
        perror("listen"); 
        exit(EXIT_FAILURE); 
    }
}

pair<string, int> listen_instruction(int server_fd, struct sockaddr_in address){
	int new_socket, valread; 
    int opt = 1; 
    int addrlen = sizeof(address); 
    char buffer[1024] = {0}; 



    if ((new_socket = accept(server_fd, (struct sockaddr *)&address,  
                       (socklen_t*)&addrlen))<0) 
    { 
        perror("accept"); 
        exit(EXIT_FAILURE); 
    } 
    valread = read( new_socket , buffer, 1024); 
    
	//char *hello = "Hello from server"; 
	//send(new_socket , hello , strlen(hello) , 0 ); 
    //printf("Hello message sent\n"); 
	buffer[strlen(buffer)-1]=0;


    return make_pair(string(buffer),new_socket); 
} 


set<int> get_angles(string file, int n){
  ifstream _file(file.c_str(), ios::in);
  string line;

  if (! _file)
    cerr << "ERROR: unable to open instance file: " << file << ", stopping! \n";

  //cout << "##Reading volume files." << endl;
  getline(_file, line);
  _file.close();

  set<int> angles;
  stack<string> q;
  char delim=' ';
  std::size_t current, previous = 0;
  current = line.find(delim);
  //cout << line << endl;
  while (current != std::string::npos) {
    angles.insert(atoi(line.substr(previous, current - previous).c_str()));
    previous = current + 1;
    current = line.find(delim, previous);
  }
  angles.insert(atoi(line.substr(previous, current - previous).c_str()));

  return angles;
}

list<int> get_angles(string str, string& organ_name) 
{ 
    list<int> angles;
    std::replace(str.begin(), str.end(), '_', ' ');
    std::replace(str.begin(), str.end(), '-', ' ');
    std::replace(str.begin(), str.end(), '/', ' ');
    stringstream ss;     
  
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
            angles.push_back(found); 
            flag=true;
        }else if(flag==true){
          stringstream(temp) >> organ_name;
          //flag=false;
        }
  
        /* To save from space at the end of string */
        temp = ""; 
    } 
    return angles;
} 


vector<Volume> createVolumes (string organ_filename, Collimator& collimator, int max_voxels){
  ifstream organ_file(organ_filename.c_str(), ios::in);
  vector<string> organ_files;
  vector<Volume> volumes;
  string line;

  if (! organ_file)
    cerr << "ERROR: unable to open instance file: " <<
             organ_filename << ", stopping! \n";

  cout << "##Reading volume files." << endl;
  //getline(organ_file, line); //--> angles

  string organ_name = "";
  bool read=false;
  for(int j=0;j<3;j++){
    cout << "new organ" << endl;
    Volume v(collimator, "");
  
    while (true) {
      if(!read && !getline(organ_file, line)) break;
      string tmp;
      list<int> angles = get_angles(line, tmp);
      cout << tmp << endl;
      cout << "##  " << line << endl;

      if(organ_name=="" || tmp == organ_name){
        organ_name=tmp;
        v.set_data(line, max_voxels, angles);
        read=false;
      }else {
        organ_name=tmp;
        read=true;
        break;
      }


    }
    volumes.push_back(v);
  }
  organ_file.close();

  return(volumes);
}


int main(int argc, char** argv){

	int server_fd; struct sockaddr_in address;



  int seed = time(NULL);

  // Budget and execution variables
  int maxiter = 5000;
  int ibo_evals = 2500;
  int maxtime = 0;
  int maxeval = 0; 

  // Aperture and intensity configuration
  StationSetup initial_setup = StationSetup::open_min_k;;
  int max_apertures = 5;
  int I_0 = 5;
  int max_intensity = 28;
  int step_intensity = 1;

  // Type of local search
  string strategy = "dao_ls";
  bool continuous = true;
  double prob_intensity = 0.2;
  int tabu_size = 0;

  // Target beamlet heuristic
  bool targeted_search = false;
  double vsize = 0.002;
  double min_improvement=0.05;


  // Files
  string path = ".";


  string file2 = "data/test_instance_coordinates.txt";
  char* file3 = NULL;
  int max_voxels=100000;
  int port = 8080;


  args::ArgumentParser parser("********* IMRT-Solver (Aperture solver) *********",
                             "Example.\n../AS -s ibo_ls --maxeval=4000 --ls_sequential=intensity --setup=open_min --seed=2 --ls=first  --max-intensity=20");

  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::ValueFlag<int>    _seed     (parser, "int", "Seed  (" +
                                     to_string(seed)+")", {"seed"});
  args::ValueFlag<int>    _port    (parser, "int", "port  (" +
                                     to_string(seed)+")", {"port"});

  // Search strategies
  args::Group strat (parser, "Strategy options:");
  args::ValueFlag<int>    _tabu_size (strat, "int",
                                    "Tabu list size(" +
                                     to_string(tabu_size)+")", {"tabu-size"});


  // Execution parameters
  args::Group budget (parser, "Budget options:");
  args::ValueFlag<int>    _maxtime  (budget, "int",
                                    "Maximum time in seconds (" +
                                     to_string(maxtime)+")", {"maxtime"});
  args::ValueFlag<int>    _maxeval  (budget, "int",
                                    "Number of evaluations (" +
                                     to_string(maxiter)+ ")", {"maxeval"});


  // Initial collimator setup (initial solution: aperture, intensity)
  args::Group isetup (parser, "Initial collimator setup:");
  args::ValueFlag<int> _I_0 (isetup, "int",
                                  "Initial intensities in the fluence map matrices " + to_string(I_0), {"setup"});
  args::ValueFlag<int> _max_apertures     (isetup, "int",
                                          "Number of apertures per angle (station) (" +
                                           to_string(max_apertures)+")", {"max-apertures"});

  // Target beamlet heuristic parameters
  args::Group heur (parser, "Parameters for the ordered neighbourghood (targeted):");
  args::ValueFlag<double>    _vsize    (heur, "double",
                                    "Percentage of considered worst voxels (" +
                                     to_string(vsize)+")", {"vsize"});
 args::ValueFlag<double>    _min_improvement    (heur, "double",
                                   "Minimum beamlet improvement estimation(" +
                                    to_string(min_improvement)+")", {"min_impr"});

  // Problem file parameters
  args::Group io_opt (parser, "Input output options:");
  args::ValueFlag<string> _file  (io_opt, "string",
                                 "Files with the deposition matrices", {"files-dep"});
  args::ValueFlag<string> _file2 (io_opt, "string",
                                 "File with the beam coordinates", {"file-coord"});
  args::ValueFlag<string> _path  (io_opt, "string",
                                 string("Absolute path of the executable ") +
                                 "(if it is executed from other directory)", {"path"}); 
  args::ValueFlag<string> _bac  (io_opt, "string",
                                 "Beam angle configuration ", {"bac"});     
  args::ValueFlag<string> _fm  (io_opt, "string",
                                 "Initial fluence map ", {"fm"});  
  args::ValueFlag<int> _max_voxels  (io_opt, "int",
                                 "Maximum number of voxels per organ ", {"max_voxels"});                        

  // Output file parameters
  args::ValueFlag<string> _convergence_file (io_opt, "string",
                                 "File to output convergence", {"convergence"});

	try { parser.ParseCLI(argc, argv); }
	catch (args::Help&) {  std::cout << parser; return 0; }
	catch (args::ParseError& e)  { std::cerr << e.what() << std::endl;  std::cerr << parser; return 1; }
	catch (args::ValidationError& e) { std::cerr << e.what() << std::endl; std::cerr << parser; return 1; }

  //setting parameters

	if(_maxtime) maxtime = _maxtime.Get();
  if(_maxeval) maxeval = _maxeval.Get();
	if(_seed) seed = _seed.Get();

  // parameters for the ordered neighbourhood (targeted)
	if(_vsize) vsize = _vsize.Get();
  if(_min_improvement) min_improvement=_min_improvement.Get();

  if(_max_apertures) max_apertures=_max_apertures.Get();
  if(_I_0) I_0=_I_0.Get();

  // Local search strategy
  if(_tabu_size) tabu_size = _tabu_size.Get();
  if(_port) port= _port.Get();
 


  if (_file2) file2=_file2.Get();
  if (_path) path=_path.Get();

  int aux=chdir(path.c_str());

  string base_name="output";
  if (_convergence_file) base_name=string("output/") + _convergence_file.Get();
  string convergence_file = base_name + "_" + to_string(seed) + ".conv";

	initialize_socket(server_fd, address, port);
  pair<string, int> message = listen_instruction(server_fd, address);

  // Iniciar generador de numeros aleatorios
  srand(seed); 

  vector<double> w={1,1,5};
  vector<double> Zmin={0,0,76};
  vector<double> Zmax={65,65,76};

  IntensityILS2::vsize=vsize;
  IntensityILS2::min_improvement=min_improvement;

  if(_max_voxels) max_voxels=_max_voxels.Get();

  int angint[]= {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,
                120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,
                215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,
                310,315,320,325,330,335,340,345};
  set<int> angles(angint,angint+70);
  string file = _file.Get();

  Collimator collimator (file2, angles); 
  vector<Volume> volumes = createVolumes (file, collimator, max_voxels);
  //cout << volumes[0].get_n_matrices() << endl;

  cout << "creating evaluation function" << endl;
  EvaluationFunction& ev = EvaluationFunction::getInstance(volumes, collimator);


  /**************** UNTIL HERE: INIT INSTANCE ****************/
  Plan* P= NULL;
  double cost;
  double used_evaluations = 0;
  std::clock_t begin_time = clock();

  std::ostringstream response2;
  response2 << "ready" << endl;
  send(message.second , response2.str().c_str() , response2.str().size() , 0 ); 
	close(message.second);

  server:
  std::ostringstream response;
  response.clear();

  message = listen_instruction(server_fd, address);
  istringstream message_stream(message.first);
  string instruction; message_stream >> instruction;
  

  if(instruction == "get_info"){
    for (auto v : volumes)
      response << v.getNbVoxels() << " ";
    response << endl;

    for (auto angle : collimator.getAngles())
      response << angle << " ";
    response << endl;

    response << "{" ;
    for (auto angle : collimator.getAngles())
      response << angle << ": " << collimator.getNangleBeamlets(angle) << ", ";
    
    response << "}" << endl;

    //matrices shapes

    response << collimator.getXdim() << " " << collimator.getYdim() << endl;
    //valid shapes of matrices
    for (auto angle : collimator.getAngles()){
      for (int i=0; i<collimator.getXdim(); i++) {
        for (int j=0; j<collimator.getYdim(); j++) {
          if (collimator.isActiveBeamAngle(i,j,angle)) {
            response << 1 << " ";
          } else {
            response << 0 << " ";
          }
        }
      }
      response << endl;
    }

  }else if(instruction == "init_fluence_map"){
    cout << instruction << endl;
    /*** init_fluence_map n_angles angle1 angle2 ... angleN fm1 fm2 fm3... fmn ***/
    int n_angles; message_stream >> n_angles;
    vector<int> bac(n_angles);
    for (int i=0; i<n_angles; i++){
       int angle; message_stream >> angle; bac[i] = angle;
    }

    istringstream* fluence_map = &message_stream;
    if(P!= NULL) delete P;
    P = new Plan (w, Zmin, Zmax, collimator, volumes, max_apertures,
            max_intensity, I_0, step_intensity,
            initial_setup, fluence_map, bac);

    response << P->getEvaluation();

  }else if(instruction == "local_search"){
    /** local_search neigh maxeval **/
    cout << "local_search" << endl;
    NeighborhoodType n_type;
    string neigh; message_stream >> neigh;
    if(neigh=="beam_intensity") n_type=NeighborhoodType::aperture;
    else if(neigh=="level_intensity") n_type=NeighborhoodType::intensity;
    else if(neigh=="mixed") n_type=NeighborhoodType::mixed;
    int maxeval; message_stream >> maxeval;

    IntensityILS2 ibo;
    cost = ibo.iteratedLocalSearch(*P, maxtime, maxeval, LSType::first, false /*continuous*/,
              n_type, LSTargetType::target_friends, PerturbationType::none, 
              0 /*perturbation_size*/, tabu_size, convergence_file, 0, begin_time, true /*verbose*/);
    used_evaluations=ibo.used_evaluations;

    response << cost << " " << used_evaluations; 

  }else if(instruction == "perturbation"){
    IntensityILS2 ibo;
    ibo.perturbation(*P, PerturbationType::p_mixed, 3, false /*verbose*/);
    response << "perturbation ready" << endl;
  }else if(instruction == "get_fluence_map"){
    cout << "get_fluence_map" << endl;
    std::streambuf*     oldbuf  = std::cout.rdbuf( response.rdbuf() ); //para retornar cout..
	  for(int i=0;i<P->getNStations();i++){
		  P->printIntensity(i, true);
      cout << endl ; 
    }
    std::cout.rdbuf(oldbuf);
  }else if(instruction == "get_impact_map"){
    cout << "get_impact_map" << endl;
    std::streambuf*     oldbuf  = std::cout.rdbuf( response.rdbuf() ); //para retornar cout..
    multimap < double, pair<int, int>, MagnitudeCompare> beamlets =
    P->getEvaluationFunction()->best_beamlets(*P, vsize);
    map < pair<int,int>, double > beam2impact;
    for( auto b : beamlets)  beam2impact[b.second] = b.first;

	  for(int k=0;k<P->getNStations();k++){
        for (int i=0; i<collimator.getXdim();i++) {
			    for (int j=0; j<collimator.getYdim(); j++) {
            Station* s = P->get_station(k);
            
            if (s->pos2beam.find(make_pair(i,j))!=s->pos2beam.end()){
              int beam = s->pos2beam[make_pair(i,j)];
              if( beam2impact.find(make_pair(k,beam))!= beam2impact.end()  && beam2impact[make_pair(k,beam)] > abs(min_improvement))
                cout << beam2impact[make_pair(k,beam)] << " ";
              else cout << "-1 ";
            }
			    }
		  }
      cout << endl ; 
    }
    std::cout.rdbuf(oldbuf);
  }else if(instruction == "get_dose_vector"){
    vector<vector<double>> Z = ev.get_Z();
    for(vector<double> organZ : Z){
      for(double z : organZ)
        response << z << " ";
      response << endl;
    }
  }else if(instruction == "get_deposition_matrix"){
    int organ;  message_stream >> organ;
    int angle;  message_stream >> angle;
    Matrix& D = volumes[organ].getDepositionMatrix(angle); 
    for(int k=0; k< D.nb_rows(); k++){
      for(int b=0; b< D.nb_cols(); b++)
         response << D(k,b) << " "  ; 
      response << endl; 
    }

  }else if(instruction == "quit"){
    return 0;
  }

  send(message.second , response.str().c_str() , response.str().size() , 0 ); 
	close(message.second);
	cout << "listening..." << endl;

  goto server;


	return 0;

}
