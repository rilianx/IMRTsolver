IMRT Solver
-----------

### Installation


````
mkdir imrt
cd imrt
git clone https://github.com/rilianx/IMRTsolver.git .
git checkout -t origin/develop
cmake .
make
````

Then, download CERR instances from [here](https://drive.google.com/open?id=1C4V0pAilKPJz0L5JIr4QrGsf1cnF2vH0) and uncompress them in the `data` folder.

Now you can run the solver, for example, 
````
./AS -s ibo_ls --setup=open_min --ls_sequential=aperture -s ibo_ls --maxeval=15000 --ls=first --perturbation-size=5 --seed=1 --max-intensity=20 --file-coord=data/Equidistantes/equidist-coord.txt --initial-intensity=5 --obj=mpse  --file-dep=data/Equidistantes/equidist05.tx
````


### Commands

````
./AS {OPTIONS}

    ********* IMRT-Solver (Aperture solver) *********

  OPTIONS:

      -h, --help                        Display this help menu
      --seed=[int]                      Seed (1627334539)
      Strategy options:
        -s[string], --strategy=[string]   Strategy (dao_ls|ibo_ls|mixedILS)
        -l[string], --ls=[string]         Local search strategy (best|first)
        --tabu-size=[int]                 Tabu list size(0)
      Budget options:
        --maxtime=[int]                   Maximum time in seconds (0)
        --maxeval=[int]                   Number of evaluations (0)
      Initial collimator setup:
        -t[string], --setup=[string]      Initial setup
                                          (open_max|open_min|closed_max|closed_min|
                                          random|manual|open_min_min|open_min_k).
                                          * open_min_k initializes intensity of
                                          open apertures with the value of
                                          initial-intensity 2
        --max-apertures=[int]             Number of apertures per angle
                                          (station) (5)
        --initial-intensity=[int]         Initial value aperture intensity (2)
      Intensity options:
        --max-intensity=[int]             Max value aperture intensity (28)
        --step-intensity=[int]            Step size for aperture intensity (1)
      Neighborhood selection:
        --ls_simple=[string]              Simple neighborhood in local search
                                          (aperture|intensity|mixed|imixed)
        --ls_sequential=[string]          Sequential neighborhood in local
                                          search starting by
                                          (intensity|aperture)
        --ls_sequentialp=[double]         Probabilistic sequential neighborhood
                                          in local search [0,1]
      Heuristic options:
        --targeted                        Apply targeted local search
        --bsize=[int]                     Number of considered beamlets for
                                          selection (20)
        --vsize=[double]                  Percentage of considered worst voxels
                                          (0.002000)
        --min_impr=[double]               Minimum beamlet improvement
                                          estimation(0.050000)
      Perturbation:
        --perturbation=[string]           Type of perturbation to be applied
                                          (intensity|aperture|mixed)
        --perturbation-size=[int]         Perturbation size (0)
      Objective function:
        --obj=[string]                    Objective function used in the search
                                          (mpse: mean positive square error|gs:
                                          global score|relu_gs: gs with relu
                                          activation)
        --scores-file=[string]            Scores for the global score objective
                                          function (only if obj in {gs|relu_gs})
        --obj2=[string]                   Secondary objective function (just for
                                          information) (mpse|gs|relu_gs)
        --scores2-file=[string]           Scores for the secondary objective
                                          function (only if sec-obj in
                                          {gs|relu_gs})
      Input output options:
        --file-dep=[string]               File with the deposition matrix
        --file-coord=[string]             File with the beam coordinates
        --path=[string]                   Absolute path of the executable (if it
                                          is executed from other directory)
        --plot                            Generate plot and save in file
        --verbose                         Verbose
        --irace                           To configure with irace
        --convergence=[string]            File to output convergence

    Example.
    ./AS -s ibo_ls --setup=open_min --ls_sequential=aperture -s ibo_ls --maxeval=15000 --ls=first --perturbation-size=5 --seed=1 --max-intensity=20 --file-coord=data/Equidistantes/equidist-coord.txt --initial-intensity=5 --obj=mpse  --file-dep=data/Equidistantes/equidist05.tx
