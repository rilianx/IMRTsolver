import sys
import random
import os
import subprocess

filename = "parameterized_scores"+str(random.randint(100000,999999))+".txt"
f = open(filename, "w")

instance,seed,w1,w2,w3,w4,function,maxdmean,type, dose_str, n_sample = sys.argv[1:]

dose = float(dose_str)



f.write("D 95 2 "+ str(round(dose*0.95,2)) +" 0 "+ w1 + "\n")
f.write("D 00 2 0 "+ str(round(dose*1.07,2)) +" "+ w2 + "\n")
f.write("Dmean 00 0 0 "+maxdmean+ " "+ w3 + "\n")
f.write("Dmean 00 1 0 "+maxdmean+ " "+ w4)

f.close()

#with open(filename, 'r') as f:
#    print(f.read())


convergence_file= instance.split('/')[-1]+'-'.join(sys.argv[2:])
score_file=" --scores-file="+filename
if function=="mpse": score_file=""
command = "./AS -s ibo_ls --setup=open_min --ls_sequential=aperture -s ibo_ls " \
         "--maxeval=0 --ls=first --perturbation-size=0 --seed="+seed+" --max-intensity=20 "\
         "--file-coord=data/Equidistantes/equidist-coord.txt --initial-intensity=5 "\
         "--obj="+function+score_file+" --obj2=gs_relu --scores2-file=target_scores"+str(dose_str)+".txt "\
         "--path=/home/iaraya/imrt --file-dep="+instance+" --irace --convergence="+convergence_file + " --n_sample="+n_sample



if type=="scores_command": 
    print (command)
else:
    result = subprocess.getoutput(command)


    try:
        os.remove(filename)
    except FileNotFoundError:
        pass
    #print(result)


    if type=="irace":
        print(float(result.split("\n")[-1].split(",")[-1]))
    elif type=="scores":
        print(result.split("\n")[-3])
    elif type=="write_file":
        sys.argv.pop(2)
        sys.argv.pop(-3)
        file= instance.split('/')[-1][:-4]+'-'+'-'.join(sys.argv[2:])+".res"
        with open(file, 'a') as f:
            f.write(result.split("\n")[-1].split(",")[-1]+"\n")
        print(float(result.split("\n")[-1].split(",")[-1]))

