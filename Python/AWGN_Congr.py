import numpy as np
import matplotlib.pyplot as plt
import math

#Formulae g(t)
emp = []
Thre = np.linspace(-30,120,500)

for xxxx in Thre:
    emp.append((1-(1/(1+math.exp(-(xxxx+8)/5)))))

plt.plot(Thre,emp)
plt.show()
fileID = open("Values_T_TrueAWGN_Congr_-5dB.txt","w")
for i in emp:
    fileID.write(str(i))
    fileID.write("\n")
fileID.close()