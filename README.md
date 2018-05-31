IMRT Solver
-----------

**Cómo correr un ejemplo:**

Compilar:

cmake .
make

Para correr:

./TEST

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
