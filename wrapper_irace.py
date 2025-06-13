import sys
import random
import os
import subprocess


instance,seed,mode,min_delta,pr_first_neigh,pert_size,max_iter = sys.argv[1:]
#print(instance,seed,mode,min_delta,pr_first_neigh,pert_size,max_iter)

if mode=="gs_ils76": param=" --evals=eval_functions/gs_ils76.txt,eval_functions/gs_oar76.txt --sf=0 --of=1 "

if pert_size=="0": max_iter="0"

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

command = "./AS \
    --maxeval="+max_iter+ \
    " --neighborhoods=aperture,intensity \
    --min-delta="+min_delta+ \
    " --perturbation-size="+pert_size+ \
    pr_neigh_str+ \
    param+ \
    file_coord+ \
    " --file-dep="+instance+ \
    " --path=. --output-file="+conv_file+ " --output-fm=" +sol_file+" --seed=" + seed

#print(command)

result = subprocess.getoutput(command)

print(float(result.split("\n")[-1].split(":")[1]))
