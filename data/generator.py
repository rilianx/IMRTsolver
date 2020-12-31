import sys 
  
# total arguments 
n = len(sys.argv) 
print("Total arguments passed:", n) 
  
# Arguments passed 
print("\nName of Python script:", sys.argv[0]) 

name = str(sys.argv[1])

organ1 = str(sys.argv[2])
organ2 = str(sys.argv[3])
tumor = str(sys.argv[4])

first_angle = 0

for i in range(14):
    if i<10: f= open("TRT00X/"+name+"_0"+str(i)+".txt","w+")
    else: f= open("TRT00X/"+name+"_"+str(i)+".txt","w+")

    f.write(str(first_angle) + " " +
        str(first_angle+70) + " " +
        str(first_angle+140) + " "+
        str(first_angle+210) + " "+
        str(first_angle+280) + "\n")
    
    strr = str(first_angle) + "-" + str(first_angle+70) + "-" + str(first_angle+140) + "-" + str(first_angle+210) + "-" + str(first_angle+280) + "-"
    
    f.write("data/TRT00X/"+name+"_"+strr+organ1+".dat\n")
    f.write("data/TRT00X/"+name+"_"+strr+organ2+".dat\n")
    f.write("data/TRT00X/"+name+"_"+strr+tumor+".dat\n")
    f.close()

    first_angle += 5