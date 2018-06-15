IMRT Solver
-----------

**Cómo correr un ejemplo:**

Compilar:

cmake .
make

Para correr:

´´´
  ./IAS {OPTIONS}

    ********* IMRT-Solver (Intensity-aperture solver) *********

  OPTIONS:

      -h, --help                        Display this help menu
      --bsize=[int]                     Number of considered beamlets for
                                        selection (5)
      --vsize=[int]                     Number of considered worst voxels (20)
      --int0=[int]                      Initial intensity for beams (4.000000)
      --max_ap=[int]                    Initial intensity for the station (4)
      --maxdelta=[int]                  Max delta (5.000000)
      --maxratio=[int]                  Max ratio (3.000000)
      --alpha=[double]                  Initial temperature for intensities
                                        (1.000000)
      --beta=[double]                   Initial temperature for ratio (1.000000)
      --max_iter=[int]                  Number of iterations (100)

    An IMRT Solver.
´´´

Por ejemplo:

     ./IAS --max_iter=400 --maxdelta=8 --maxratio=6 --alpha=0.999 --beta=0.999 --max_ap=4 --help

----



