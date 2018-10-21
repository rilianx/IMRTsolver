IMRT Solver
-----------

**Cómo correr un ejemplo:**

Compilar:

````
cmake .
make
````

Para correr:

    ./AS {OPTIONS}

    ********* IMRT-Solver (Aperture solver) *********

    OPTIONS:

      -h, --help                        Display this help menu
      -s[string], --strategy=[string]   Strategy (dao_ls|ibo_ls)
      --bsize=[int]                     Number of considered beamlets for
                                        selection (20)
      --vsize=[int]                     Number of considered worst voxels (50)
      --int0=[int]                      Initial intensity for beams (4.000000)
      --max_ap=[int]                    Initial intensity for the station (5)
      --maxdelta=[int]                  Max delta (5.000000)
      --maxratio=[int]                  Max ratio (3.000000)
      --alpha=[double]                  Initial temperature for intensities
                                        (1.000000)
      --beta=[double]                   Initial temperature for ratio (1.000000)
      --maxiter=[int]                   Number of iterations (5000)
      --maxtime=[int]                   Maximum time in seconds (0)
      --seed=[int]                      Seed (1533566270)
      Direct aperture local search:
        --ls-aperture                     Apply aperture local search
        --ls-intensity                    Apply intensity local search
        --prob-intensity=[double]         Probability to search over intensity
                                          (0.200000)
      --open-setup                      Initialize apertures as open
      --initial-intensity=[int]         Initial value aperture intensity (2)
      --max-intensity=[int]             Max value aperture intensity (28)
      --step-intensity=[int]            Step size for aperture intensity (2)
      Acceptance criterion:
        --accept-best                     Accept only improvement
        --accept-sa                       Accept as simulated annealing
      --temperature=[double]            Temperature for acceptance criterion
                                        (0.000000)
      --alphaT=[double]                 Reduction rate of the temperature
                                        (0.950000)
      --perturbation-size=[int]         Perturbation size (2)

    Example.
    ./AS -s ibo_ls --maxiter=400 --maxdelta=8 --maxratio=6 --alpha=0.999
    --beta=0.999 --bsize=5 --vsize=20 --max_ap=4 --seed=0 --int0=1 --open-setup

----

## El Algoritmo


    S ← initializeStations (max_apertures, initial_intensity)
    bestF ← eval(S)

    while i < max_iter:

      (b,si,increase) ← select_promising_beamlet(bsize,vsize)
      delta_intensity ← random(1,maxdelta)
      maxdelta ← maxdelta*alpha
      ratio ← random(0,maxratio)
      maxratio ← maxratio*beta

      if increase:
        diff ← increaseIntensity&Repair (b,si,delta_intenisty,ratio)
      else:
        diff ← decreaseIntensity&Repair (b,si,delta_intenisty,ratio)

      F ← incremental_eval (S, diff)

      if F < bestF:
        bestF ← F
        bestS ← S
      else
        revert_changes(s,diff)
        undo_last_eval()

[Más detalles](https://docs.google.com/document/d/1EGoKoLsmik4TSiY_SslWkxddDCrYvUfFKpPLRpQPV_U/edit#)
