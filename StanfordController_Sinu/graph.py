import numpy as np
import math
import matplotlib.pyplot as plt 


CourbeY=[]
CourbeT=[]
i=1
file = open('/home/xionfire/Bureau/Japan_Project/Controller/StanfordController_Sinu/sinus.txt', "r")
linesY = file.readlines()
file.close()
for line in linesY:
    CourbeY.append(float(line.strip()))
    
file = open('/home/xionfire/Bureau/Japan_Project/Controller/StanfordController_Sinu/time.txt', "r")
linesT = file.readlines()
file.close()
for line in linesT:
    CourbeT.append(float(line.strip()))

for i in range (np.size(CourbeY)):
	if (CourbeY[i]>CourbeY[i-1]):
		Max=CourbeY[i]
	if (CourbeY[i]<CourbeY[i-1]):	
		Min=CourbeY[i]


ymax=max(CourbeY)
ymin=min(CourbeY)
print("Y min is =", ymin)
print("Y max is =", ymax)
plt.plot(CourbeT,CourbeY,"r")
plt.show()
