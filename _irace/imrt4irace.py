import sys
import random
import os
import subprocess

instance,seed,mode,min_delta,pr_first_neigh,pert_size,switch_patience = sys.argv[1:]

if mode=="of76": param=" --evals=scores/of76.txt,scores/of76.txt --sf=0 --of=0 "
if mode=="gs76": param=" --evals=scores/gs76.txt,scores/of76.txt --sf=0 --of=1 "

conv_file= "output/"+instance.split('/')[-1]+'-'+'-'.join(sys.argv[2:])

command = "../AS \
    --maxeval=10000 \
    --neighborhoods=aperture,intensity \
    --min-delta="+min_delta+ \
    " --switch-patience="+switch_patience+ \
    " --perturbation-size="+pert_size+ \
    " --pr-first-neigh="+pr_first_neigh+ \
    param+ \
    "--file-coord=data/Equidistantes/equidist-coord.txt \
    --file-dep="+instance+ \
    " --path=.. --output-file="+conv_file+" --seed="+seed

print(command)

result = subprocess.getoutput(command)

print(float(result.split("\n")[-1].split(":")[1]))
