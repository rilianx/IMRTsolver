IMRTsolver (Global Score)
==

## Descarga e instalación

```
git clone https://github.com/rilianx/IMRTsolver.git .
git checkout -t origin/gscore
cmake .
make
```

### Descarga de instancias de prueba

Además puedes descargar instancias de prueba de aquí:

- [Instancias TRT](https://drive.google.com/file/d/1b0SSVEScIgdbimNFrW1u8QvY_McogxLg/view?usp=sharing)
- [Instancias CERR](https://drive.google.com/file/d/1C4V0pAilKPJz0L5JIr4QrGsf1cnF2vH0/view?usp=sharing)

Descomprimelas dentro de la carpeta `data`

```
cd data
tar -xzf equidistant-instances.tar.gz
tar -xzf TRT00X-instances.tar.gz
```

## Comando ejecución

### Ejemplo

```bash
#in fwk4exp
#python3 imrt4irace.py instance seed obj epsilon n_bi pert_size n_evals
./AS --maxeval=10000 --path=. --seed=2 \
--neighborhoods=aperture,intensity     \
--epsilon=0.0001 \
--perturbation-size=3 \
--pr-neigh=0.2,1.0 \
--evals=eval_functions/gs76.txt,eval_functions/gs_oar76.txt --sf=0 --of=1 \
--file-coord=data/Equidistantes/equidist-coord.txt \
--file-dep=data/Equidistantes/equidist00.txt \
--output-file=convergence_file.txt \
--output-fm=output/fluence_map_solution.txt
```

### Options

```
-h, --help                        Display this help menu
--seed=[int]                      Seed (1749818971)
**Budget options:**
  --maxeval=[int]                   Number of evaluations (0)
**Initial collimator setup:**
  --max-apertures=[int]             Number of apertures per angle
                                    (station) (5)
**Neighborhood selection:**
  --neighborhoods=[string]          neighborhoods in local search (intensity|aperture)
**Acceptation improvement:**
  --epsilon=[float]                 Minimum delta eval for accepting the
                                    change
  --pr-neigh=[string]               Prop. of elements of each
                                    neighbourhood (1.0 by default)
**Perturbation:**
  --perturbation-size=[int]         Perturbation size (0)
**Evaluators:**
  --evals=[string]                  Files with evaluation functions (global score).
  --sf=[int]                        index of the objective function used 
																	  in ILS (last index=z funct)
--of=[int]                        index of the objective function  (last index=z funct)
**Input/Output options:**
  --file-dep=[string]               File with the deposition matrix
  --file-coord=[string]             File with the beam coordinates
  --path=[string]                   Absolute path of the executable (if it
                                    is executed from other directory)
  --output-file=[string]            File to output all indicators for each
                                    iteration (convergence)
  --output-fm=[string]              File to output the fluence map of
                                    voxels (best solution)
  --verbose                         Verbose
```

### irace wrapper

**Instalación**

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

```

**Ejecución**:

```
# python3 wrapper_irace.py <instance> <seed> <mode> <epsilon> <nbi> <psize> <max_eval>
# Example: 
python3 wrapper_irace.py data/Equidistantes/equidist00.txt 1 gs_ils76 0.0001 1.0 3 10000
```

Solo está implementado el modo gs_ils76:

- **Función objetivo ILS**: gs_ils76
- **Función objetivo final**: gs_oar76
