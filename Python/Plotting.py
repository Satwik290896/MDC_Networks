import numpy as np
import matplotlib.pyplot as plt
from numpy import matlib as mb

x = range(500)
N = 5
N2 = 10
N3 = 20
Thre = np.linspace(-30,120,500)
f1 = open("Values_T_Model1_-5dB.txt",'r')
f1_lines = f1.readlines()
f2 = open("Values_T_Model1_Rayleigh_-5dB.txt",'r')
f2_lines = f2.readlines()

f11 = open("Values_T_Model2_-5dB.txt",'r')
f11_lines = f11.readlines()
f22 = open("Values_T_Model2_Rayleigh_-5dB.txt",'r')
f22_lines = f22.readlines()

fA = open("Values_T_TrueAWGN_-5dB.txt",'r')
fA_lines = fA.readlines()
fR = open("Values_T_TrueRayleigh_-5dB.txt",'r')
fR_lines = fR.readlines()
  
fA2 = open("Values_T_TrueAWGN_Congr_-5dB.txt",'r')
fA2_lines = fA2.readlines()
fR2 = open("Values_T_TrueRayleigh_Congr_-5dB.txt",'r')
fR2_lines = fR2.readlines()


F1_Arr = []
F2_Arr = []
F11_Arr = []
F22_Arr = []
FA_Arr = []
FR_Arr = []
FA2_Arr = []
FR2_Arr = []

for line in f1_lines:
    F1_Arr.append(float(((line.strip()).split())[2]))

A = np.array([F1_Arr])
np.reshape(A,(1,500))

F1_Arr_rep = mb.repmat(A, 500, 1) 
F1_Arr_rep = F1_Arr_rep.transpose()
R1 = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(F1_Arr_rep[i,:])-np.array(F1_Arr_rep[i+N,:])))
    R1.append(1- (k/500))
R1 = np.array(R1)

R1_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(F1_Arr_rep[i,:])-np.array(F1_Arr_rep[i+N2,:])))
    R1_N2.append(1- (k/500))
R1_N2 = np.array(R1_N2)

R1_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(F1_Arr_rep[i,:])-np.array(F1_Arr_rep[i+N3,:])))
    R1_N3.append(1- (k/500))
R1_N3 = np.array(R1_N3)

x1 = [Thre[i] for i in range(np.size(R1))]
x2 = [Thre[i] for i in range(np.size(R1_N2))]
x3 = [Thre[i] for i in range(np.size(R1_N3))]






for line in f2_lines:
    F2_Arr.append(float(((line.strip()).split())[0]))

A = np.array([F2_Arr])
np.reshape(A,(1,500))

F2_Arr_rep = mb.repmat(A, 500, 1) 
F2_Arr_rep = F2_Arr_rep.transpose()
R2 = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(F2_Arr_rep[i,:])-np.array(F2_Arr_rep[i+N,:])))
    R2.append(1- (k/500))
R2 = np.array(R2)

R2_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(F2_Arr_rep[i,:])-np.array(F2_Arr_rep[i+N2,:])))
    R2_N2.append(1- (k/500))
R2_N2 = np.array(R2_N2)

R2_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(F2_Arr_rep[i,:])-np.array(F2_Arr_rep[i+N3,:])))
    R2_N3.append(1- (k/500))
R2_N3 = np.array(R2_N3)






for line in f11_lines:
    F11_Arr.append(float(((line.strip()).split())[0]))
    
A = np.array([F11_Arr])
np.reshape(A,(1,500))

F11_Arr_rep = mb.repmat(A, 500, 1) 
F11_Arr_rep = F11_Arr_rep.transpose()
R11 = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(F11_Arr_rep[i,:])-np.array(F11_Arr_rep[i+N,:])))
    R11.append(1- (k/500))
R11 = np.array(R11)

R11_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(F11_Arr_rep[i,:])-np.array(F11_Arr_rep[i+N2,:])))
    R11_N2.append(1- (k/500))
R11_N2 = np.array(R11_N2)

R11_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(F11_Arr_rep[i,:])-np.array(F11_Arr_rep[i+N3,:])))
    R11_N3.append(1- (k/500))
R11_N3 = np.array(R11_N3)





for line in f22_lines:
    F22_Arr.append(float(((line.strip()).split())[0]))

A = np.array([F22_Arr])
np.reshape(A,(1,500))

F22_Arr_rep = mb.repmat(A, 500, 1) 
F22_Arr_rep = F22_Arr_rep.transpose()
R22 = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(F22_Arr_rep[i,:])-np.array(F22_Arr_rep[i+N,:])))
    R22.append(1- (k/500))
R22 = np.array(R22)

R22_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(F22_Arr_rep[i,:])-np.array(F22_Arr_rep[i+N2,:])))
    R22_N2.append(1- (k/500))
R22_N2 = np.array(R22_N2)

R22_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(F22_Arr_rep[i,:])-np.array(F22_Arr_rep[i+N3,:])))
    R22_N3.append(1- (k/500))
R22_N3 = np.array(R22_N3)







for line in fA_lines:
    FA_Arr.append(float(((line.strip()).split())[0]))

A = np.array([FA_Arr])
np.reshape(A,(1,500))

FA_Arr_rep = mb.repmat(A, 500, 1) 
FA_Arr_rep = FA_Arr_rep.transpose()
RA = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(FA_Arr_rep[i,:])-np.array(FA_Arr_rep[i+N,:])))
    RA.append(1- (k/500))
RA = np.array(RA)

RA_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(FA_Arr_rep[i,:])-np.array(FA_Arr_rep[i+N2,:])))
    RA_N2.append(1- (k/500))
RA_N2 = np.array(RA_N2)

RA_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(FA_Arr_rep[i,:])-np.array(FA_Arr_rep[i+N3,:])))
    RA_N3.append(1- (k/500))
RA_N3 = np.array(RA_N3)









for line in fR_lines:
    FR_Arr.append(float(((line.strip()).split())[0]))

A = np.array([FR_Arr])
np.reshape(A,(1,500))

FR_Arr_rep = mb.repmat(A, 500, 1) 
FR_Arr_rep = FR_Arr_rep.transpose()
RR = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(FR_Arr_rep[i,:])-np.array(FR_Arr_rep[i+N,:])))
    RR.append(1- (k/500))
RR = np.array(RR)

RR_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(FR_Arr_rep[i,:])-np.array(FR_Arr_rep[i+N2,:])))
    RR_N2.append(1- (k/500))
RR_N2 = np.array(RR_N2)

RR_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(FR_Arr_rep[i,:])-np.array(FR_Arr_rep[i+N3,:])))
    RR_N3.append(1- (k/500))
RR_N3 = np.array(RR_N3)







for line in fA2_lines:
    FA2_Arr.append(float(((line.strip()).split())[0]))

A = np.array([FA2_Arr])
np.reshape(A,(1,500))

FA2_Arr_rep = mb.repmat(A, 500, 1) 
FA2_Arr_rep = FA2_Arr_rep.transpose()
RA2 = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(FA2_Arr_rep[i,:])-np.array(FA2_Arr_rep[i+N,:])))
    RA2.append(1- (k/500))
RA2 = np.array(RA2)

RA2_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(FA2_Arr_rep[i,:])-np.array(FA2_Arr_rep[i+N2,:])))
    RA2_N2.append(1- (k/500))
RA2_N2 = np.array(RA2_N2)

RA2_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(FA2_Arr_rep[i,:])-np.array(FA2_Arr_rep[i+N3,:])))
    RA2_N3.append(1- (k/500))
RA2_N3 = np.array(RA2_N3)









for line in fR2_lines:
    FR2_Arr.append(float(((line.strip()).split())[0]))

A = np.array([FR2_Arr])
np.reshape(A,(1,500))

FR2_Arr_rep = mb.repmat(A, 500, 1) 
FR2_Arr_rep = FR2_Arr_rep.transpose()
RR2 = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(FR2_Arr_rep[i,:])-np.array(FR2_Arr_rep[i+N,:])))
    RR2.append(1- (k/500))
RR2 = np.array(RR2)

RR2_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(FR2_Arr_rep[i,:])-np.array(FR2_Arr_rep[i+N2,:])))
    RR2_N2.append(1- (k/500))
RR2_N2 = np.array(RR2_N2)

RR2_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(FR2_Arr_rep[i,:])-np.array(FR2_Arr_rep[i+N3,:])))
    RR2_N3.append(1- (k/500))
RR2_N3 = np.array(RR2_N3)









A_temp = range(0, np.size(RA), 3)
A_temp_N2 = range(0, np.size(RA_N2), 3)
A_temp_N3 = range(0, np.size(RA_N3), 3)

x1_crisp =[x1[i] for i in A_temp]
x2_crisp =[x2[i] for i in A_temp_N2]
x3_crisp =[x3[i] for i in A_temp_N3]
RA_crisp = [RA[i] for i in A_temp]
RR_crisp = [RR[i] for i in A_temp]
RA_N2_crisp = [RA_N2[i] for i in A_temp_N2]
RR_N2_crisp = [RR_N2[i] for i in A_temp_N2]
RA_N3_crisp = [RA_N3[i] for i in A_temp_N3]
RR_N3_crisp = [RR_N3[i] for i in A_temp_N3]


RA2_crisp = [RA2[i] for i in A_temp]
RR2_crisp = [RR2[i] for i in A_temp]
RA2_N2_crisp = [RA2_N2[i] for i in A_temp_N2]
RR2_N2_crisp = [RR2_N2[i] for i in A_temp_N2]
RA2_N3_crisp = [RA2_N3[i] for i in A_temp_N3]
RR2_N3_crisp = [RR2_N3[i] for i in A_temp_N3]





fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(x1,R1,'k-',label='Idempotent/Congruency (VLDC/AWGN extremity)')
axs[0, 0].plot(x1,R11,'k--',label='Idempotent/Congruency (UDC/Rayleigh extremity)')
axs[0, 0].plot(x1_crisp,RA_crisp,'kx',label='VLDC/AWGN (Ideal fit)')
axs[0, 0].plot(x1,RR,'k:',label='UDC/Rayleigh (Ideal fit)')
axs[0, 0].set_title('Idempotent (extremety vs Ideal fit)\n(a)')
axs[0, 0].set_ylabel('N=5 Correlation')
#axs[0, 0].legend()
axs[0, 0].set_xlim([-30,5])
axs[0, 0].set_ylim([0.4,1.1])

axs[1, 0].plot(x2,R1_N2,'k-')
axs[1, 0].plot(x2,R11_N2,'k--')
axs[1, 0].plot(x2_crisp,RA_N2_crisp,'kx')
axs[1, 0].plot(x2,RR_N2,'k:')
axs[1, 0].set_ylabel('N=10 Correlation')
axs[1, 0].set_title('(b)')
#axs[1, 0].legend()
axs[1, 0].set_xlim([-30,5])
axs[1, 0].set_ylim([0.4,1.1])

axs[2, 0].plot(x3,R1_N3,'k-')
axs[2, 0].plot(x3,R11_N3,'k--')
axs[2, 0].plot(x3_crisp,RA_N3_crisp,'kx')
axs[2, 0].plot(x3,RR_N3,'k:')
axs[2, 0].set_ylabel('N=20 Correlation')
#axs[2, 0].legend()
axs[2, 0].set_xlim([-30,5])
axs[2, 0].set_ylim([0.4,1.1])
axs[2, 0].set_title('(c)')

axs[0, 1].plot(x1,R2,'k-')
axs[0, 1].plot(x1,R22,'k--')
axs[0, 1].plot(x1_crisp,RA2_crisp,'kx')
axs[0, 1].plot(x1,RR2,'k:')
axs[0, 1].set_title('Congruency (extremety vs Ideal fit)\n(d)')
#axs[0, 1].legend()
axs[0, 1].set_xlim([-30,5])
axs[0, 1].set_ylim([0.4,1.1])

axs[1, 1].plot(x2,R2_N2,'k-')
axs[1, 1].plot(x2,R22_N2,'k--')
axs[1, 1].plot(x2_crisp,RA2_N2_crisp,'kx')
axs[1, 1].plot(x2,RR2_N2,'k:')
#axs[1, 1].legend()
axs[1, 1].set_xlim([-30,5])
axs[1, 1].set_ylim([0.4,1.1])
axs[1, 1].set_title('(e)')

axs[2, 1].plot(x3,R2_N3,'k-')
axs[2, 1].plot(x3,R22_N3,'k--')
axs[2, 1].plot(x3_crisp,RA2_N3_crisp,'kx')
axs[2, 1].plot(x3,RR2_N3,'k:')
#axs[2, 1].legend()
axs[2, 1].set_xlim([-30,5])
axs[2, 1].set_ylim([0.4,1.1])
axs[2, 1].set_title('(f)')

fig.legend(loc='lower right')
fig.tight_layout()
plt.savefig("Figure-1.png")
plt.show()