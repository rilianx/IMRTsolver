IMRT Solver
-----------

**Cómo correr un ejemplo:**

Compilar:

cmake .
make

Para correr:

./IAS {OPTIONS}

    ********* IMRT-Solver *********

  OPTIONS:

      -h, --help                        Display this help menu
      --bsize=[int]                     Number of considered beamlets for
                                        selection (5)
      --vsize=[int]                     Number of considered worst voxels (20)
      --int0=[int]                      Initial intensity for beams (4.000000)
      --maxiter=[int]                   Number of iterations (100)

    An IMRT Solver.

Por ejemplo:
./IAS --vsize=30 --bsize=10 --int0=3 --maxiter=1000 --help

----

**Directories:**


**main:** Para mains

**src:** Aquí val el código!

**data:** Aquí irán las instancias!


----

**Classes:**

**Plan:**
 An IMRT plan
 It consists in a list of stations, each station corresponds to an angle, an aperture and an intensity

**Station:**
 An IMRT station
 A station consists in a beam (matrix of beamlets), an angle, an aperture
 and an intensisty

**EvaluationFunction:**
 The evaluation function based on [this paper](https://drive.google.com/file/d/1YfMNk4GhBK97gSQ0nvpJAnyM6A3EPv61/view).
 Basically it penalizes voxels with doses larger than Zmax (healthy organs), and with doses smaller than Zmin (tumor).

 **Collimator**

 **Volume**

 ----

**Command for generating documentation:**
cldoc generate -std=c++11 -stdlib=libc++ -Isrc/maths -- src/*/* src/*  --basedir src/  --output docs/
