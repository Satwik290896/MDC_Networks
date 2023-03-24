import numpy as np
import matplotlib.pyplot as plt
from numpy import matlib as mb

N = 5
N2 = 10
N3 = 15
N4 = 20
Thre = np.linspace(-30,120,500)
sig = open("Values_T_Sigmoid_N.txt",'r')
sig_lines = sig.readlines()
SIG_plot = []

for line in sig_lines:
    SIG_plot.append(2*(float(((line.strip()).split())[0])))

A = np.array([SIG_plot])
np.reshape(A,(1,500))

SIG_plot_rep = mb.repmat(A, 500, 1) 
SIG_plot_rep = SIG_plot_rep.transpose()
R1 = []
for i in range(0, 500 - N):
    k = sum(abs(np.array(SIG_plot_rep[i,:])-np.array(SIG_plot_rep[i+N,:])))
    R1.append(1- (k/500))
R1 = np.array(R1)

R1_N2 = []
for i in range(0, 500 - N2):
    k = sum(abs(np.array(SIG_plot_rep[i,:])-np.array(SIG_plot_rep[i+N2,:])))
    R1_N2.append(1- (k/500))
R1_N2 = np.array(R1_N2)

R1_N3 = []
for i in range(0, 500 - N3):
    k = sum(abs(np.array(SIG_plot_rep[i,:])-np.array(SIG_plot_rep[i+N3,:])))
    R1_N3.append(1- (k/500))
R1_N3 = np.array(R1_N3)

R1_N4 = []
for i in range(0, 500 - N4):
    k = sum(abs(np.array(SIG_plot_rep[i,:])-np.array(SIG_plot_rep[i+N4,:])))
    R1_N4.append(1- (k/500))
R1_N4 = np.array(R1_N4)

x1 = [Thre[i] for i in range(np.size(R1))]
x2 = [Thre[i] for i in range(np.size(R1_N2))]
x3 = [Thre[i] for i in range(np.size(R1_N3))]
x4 = [Thre[i] for i in range(np.size(R1_N4))]

fig, axs = plt.subplots(1, 2)
axs[0].plot(Thre,SIG_plot,'k-',label='$g(t)=2(1-sigmoid(t)); T=10log_{e}t$')
axs[0].set_title('g(t)')
axs[0].set_ylabel('$P(SINR_{m}>T)$')
axs[0].set_xlabel('T')
axs[0].legend()
axs[0].set_xlim([-30,15])

axs[1].plot(x1,R1,'k-',label='N=5 step')
axs[1].plot(x2,R1_N2,'k--',label='N=10 step')
axs[1].plot(x3,R1_N3,'k-.',label='N=15 step')
axs[1].plot(x4,R1_N4,'k:',label='N=20 step')
axs[1].set_title('N-step Correlation graph')
axs[1].set_ylabel('N-step Correlation')
axs[1].set_xlabel('T')
axs[1].legend()
axs[1].set_xlim([-30,15])

fig.tight_layout()
plt.savefig("GenPractice.png")
plt.show()