import sys
import random
import os
import subprocess



filename = "parameterized_scores"+str(random.randint(10000,99999))+".txt"
f = open(filename, "w")

instance,seed,w1,w2,w3,w4,function,maxdmean = sys.argv[1:]

f.write("D 95 2 64.6 0 "+ w1 + "\n")
f.write("D 00 2 0 72.76 "+ w2 + "\n")
f.write("Dmean 00 0 0 "+maxdmean+ " "+ w3 + "\n")
f.write("Dmean 00 1 0 "+maxdmean+ " "+ w4 + "\n")

f.close()

#with open(filename, 'r') as f:
#    print(f.read())


convergence_file= instance.split('/')[-1]+'-'.join(sys.argv[2:])

result = subprocess.getoutput("./AS -s ibo_ls --setup=open_min --ls_sequential=aperture -s ibo_ls " \
         "--maxeval=1000 --ls=first --perturbation-size=5 --seed="+seed+" --max-intensity=20 "\
         "--file-coord=data/Equidistantes/equidist-coord.txt --initial-intensity=5 "\
         "--obj="+function+" --scores-file="+filename+" --obj2=gs_relu --scores2-file=target_scores68.txt "\
         "--path=/home/ignacio/imrt --file-dep="+instance+" --irace --convergence="+convergence_file)


os.remove(filename)
print(result)

print(float(result.split("\n")[-1].split(",")[-1]))

