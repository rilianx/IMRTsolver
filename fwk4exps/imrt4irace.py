import sys
import random
import os
import subprocess

instance,seed,mode = sys.argv[1:4]

if mode=="mse76":
    instance,seed,mode,min_delta,pr_first_neigh,pert_size,switch_patience,w1,w2,w3,z1,z2,z3min,z3max = sys.argv[1:]
else:
    instance,seed,mode,min_delta,pr_first_neigh,pert_size,switch_patience = sys.argv[1:]

if mode=="of76": param=" --evals=scores/of76.txt,scores/of76.txt --sf=0 --of=0 "
if mode=="gs76" or mode=="gs76_4": param=" --evals=scores/gs76.txt,scores/of76.txt --sf=0 --of=1 "
if mode=="gs2_76_4": param=" --evals=scores/gs2-76.txt,scores/of76.txt --sf=0 --of=1 "
if mode=="mse76": param=" --evals=scores/gs2-76.txt,scores/of76.txt --sf=2 --of=1 --w={},{},{} --z={},{},{}.{} ".format(w1,w2,w3,z1,z2,z3min,z3max)

max_iter="10000"
if pert_size=="0": max_iter="0"

pr_neigh_str=" --pr-first-neigh="+pr_first_neigh
if mode=="gs76_4" or mode=="gs2_76_4" or mode=="mse76": 
    pr_neigh_str=" --pr-neigh="+pr_first_neigh+",1.0"


conv_file= "output/"+instance.split('/')[-1]+'-'+'-'.join(sys.argv[2:])
sol_file= "output/"+instance.split('/')[-1]+'-'+'-'.join(sys.argv[2:])+".sol"
file_coord = "--file-coord=data/Equidistantes/equidist-coord.txt"

if "file-coord" in instance:
    file_coord = ""

if "TRT001" in instance:
    file_coord="--file-coord=data/TRT00X/TRT001-coord.txt"

if "TRT002" in instance:
    file_coord="--file-coord=data/TRT00X/TRT002-coord.txt"

if "TRT003" in instance:
    file_coord="--file-coord=data/TRT00X/TRT003-coord.txt"

if "TRT004" in instance:
    file_coord="--file-coord=data/TRT00X/TRT004-coord.txt"

if "TRT005" in instance:
    file_coord="--file-coord=data/TRT00X/TRT005-coord.txt"

command = "../AS \
    --maxeval="+max_iter+ \
    " --neighborhoods=aperture,intensity \
    --min-delta="+min_delta+ \
    " --switch-patience="+switch_patience+ \
    " --perturbation-size="+pert_size+ \
    pr_neigh_str+ \
    param+ \
    file_coord+ \
    " --file-dep="+instance+ \
    " --path=.. --output-file="+conv_file+ " --output-fm=" +sol_file+" --seed=" + seed

#print(command)


result = subprocess.getoutput(command)

print(float(result.split("\n")[-1].split(":")[1]))