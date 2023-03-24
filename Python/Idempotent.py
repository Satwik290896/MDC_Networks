
#First term alone taken and Condition of Model 1
#Simulation technique: 40 matched points were taken

#clear all;

import numpy
import cvxpy as cp
import math
import matplotlib.pyplot as plt
from numpy import matlib as mb
from scipy.optimize import lsq_linear
from scipy.interpolate import interp1d

def Correlating_function(npoints,med_proc,pp,XY_disp,n1):
    pproc = pp - numpy.matlib.repmat(XY_disp,npoints,1)
    A = numpy.zeros((n1,n1))

    for i in range(0,n1-1):
        A[i,i+1] = 1

    C = (numpy.array(med_proc).transpose()).tolist()
    Dt = pproc.transpose()
    b1 = numpy.zeros((n1,1)) + 0.00001
    b1[n1-1,0] = 0
    b2_t = numpy.zeros((n1,1))
    b3_t = numpy.zeros((n1,1)) + 0.00001
    D = []
    Array = numpy.zeros((n1-2,n1))
    
    C1 = []
    C1 = C
    #C1.append(Array)
    #print(i)
    print(Array)
    print(numpy.shape(C), numpy.shape(Array), numpy.shape(C1))

    for i in range(0,npoints):
        j = i % n1
        '''if j==0:
            j = n1'''
        A1 = A
        b2 = b2_t
        b3 = b3_t
        b3[j] = 10
        b2[j] = -10

        if j!=1:
            te = [1]
            te.append(numpy.zeros((1,n1-1)))
            A1[j-1][0:] = te
        print(i)
        #disp(i)
        #Xt = lsqlin(C1,[Dt(:,i);zeros(n1-2,1)],A1,b1,[],[],b2,b3)
        (__,n) = numpy.shape(C1)
        x = cp.Variable(n)
        #b = Dt[:,i].append(numpy.zeros(n1-2).transpose(),2)
        #B = numpy.append(Dt[:,i],(numpy.zeros(n1-2)).transpose(),0)
        #print("hello ", numpy.shape(B))
        #b = B.tolist()
        #objective = cp.Minimize(cp.sum_squares(C1 @ x - b))
        #constraints = [b2[0] <= x[0], b2[1] <= x[1], x[0] <= b3[0], x[1] <= b3[1]]
        #prob = cp.Problem(objective, constraints)
        B = numpy.append(Dt[:,i],(numpy.zeros(n1-2)).transpose(),0)
        print("temp b2 b3 ", b2, b3)
        '''temp = B
        temp[0] = b2[0]
        temp[1] = b2[1]
        b2 = temp'''
        b2 = numpy.ravel(b2)
        print("temp b2 b3 v2 ", b2, b3)
        
        '''temp2 = b2
        temp2[0] = b3[0]
        temp2[1] = b3[1]'''
        #b3 = temp2
        #b3 = b2
        b3 = numpy.ravel(b3)
        print("temp b2 b3 v3 ", b2, b3)
          
        #B[:,0] =  Dt[:,i]
        print(numpy.shape(B), numpy.shape(b2), numpy.shape(b3), b2, b3)
        res = lsq_linear(numpy.array(C1), B, bounds=(numpy.array(b2), numpy.array(b3)), lsmr_tol='auto', verbose=1)
        #res = b2
        print("res unhii v1: ", res)
        print("res unhii v2: ", res.keys())
        print("res unhii v3: ", type(res.items()), res.items())
        print("res unhii v4deqd: ", res['x'], type(res), numpy.shape(res))
        # The optimal objective value is returned by `prob.solve()`.
        res = res['x']
        Xt = res
        print("Xt unhii v5deqd: ", numpy.shape(Xt))
        D.append(Xt)
    
    print("D ", numpy.shape(D))
    return D

def ismember(a_vec, b_vec):
    """ MATLAB equivalent ismember function """

    bool_ind = numpy.isin(a_vec,b_vec)
    common = a_vec[bool_ind]
    common_unique, common_inv  = numpy.unique(common, return_inverse=True)     # common = common_unique[common_inv]
    b_unique, b_ind = numpy.unique(b_vec, return_index=True)  # b_unique = b_vec[b_ind]
    common_ind = b_ind[numpy.isin(b_unique, common_unique, assume_unique=True)]
    return bool_ind, common_ind[common_inv]

def argmax(lst):
  return lst.index(max(lst))

lam = 200
npoints = numpy.random.poisson(lam)
pproc = numpy.random.rand(npoints, 2) - 0.5
n1 = 2
med_proc = [[0.990000,0.990000],[0.000000001,0.000000001]]
XY_dar = 0.5*numpy.random.rand(100,2)
S = []
em_X_cod = 2 - 1
emitter_X = med_proc[em_X_cod][0:]

emp = []

D = numpy.zeros((npoints,n1,100))
P_T = 4000
freq = 1.5*pow(10,9)
wave_l = (3*10^8)/freq
noise = 0.1
Thre = range(-30,10,9)
#AWGN
SNR_avg = 0

for j in range(0,100):
    
    I = 0

    id1 = em_X_cod

    emitter = med_proc[id1][0:]
    d_er = math.sqrt(sum(numpy.multiply(emitter, emitter)))
    P_R = P_T * pow((wave_l/(4*3.14*d_er)), 2)
    
    radd = numpy.sqrt(sum(numpy.multiply(med_proc, med_proc),2))
    for i in range(0,n1):
        rad = radd[i]
        
        if (rad<0.1) and (rad!=d_er):
            alpha = 2
            I = I + P_T*(wave_l/(4*3.14*rad))^2

        if (rad>=0.1) and (rad!=d_er):
            alpha = 4
            I = I + pow(P_T*(wave_l/(4*3.14*0.1)), 2) + pow(P_T*(0.1/rad), alpha)

    SNR = P_R/(I+noise)
    SNR_avg = SNR_avg+ SNR
    
SNR_avg_AWGN = SNR_avg/100
SNR_db_avg = 10*math.log10(SNR_avg_AWGN)


T = -75
T_end = 3000
len = 10000
dif = (T_end - T)/len
Array = []
for i in range(0, 10000):
    Array.append(-75 + i*dif)
#Array = range(T,T_end,dif)
diffs = abs(Array[1]-Array[2])

R_c = 0.1
noise = 0.1

count = numpy.zeros((10000,1))

med_temp = sum((numpy.multiply(med_proc, med_proc)),2)
XPT = mb.repmat(P_T,n1,1)
YPT = numpy.zeros((npoints,100))
Xn = mb.repmat(noise,n1,1)

for i in range(0,100):
    D[:,:,i] = numpy.array(Correlating_function(npoints,med_proc,pproc,XY_dar[i,:],n1))
    S.append(numpy.sqrt(numpy.matmul(numpy.array(numpy.multiply(D[:,:,i], D[:,:,i])), numpy.array(med_temp))))
    K = numpy.matmul(numpy.array(numpy.multiply(D[:,:,i], D[:,:,i])), numpy.array(XPT))
    YPT[:,i] = numpy.matmul(numpy.array(numpy.multiply(D[:,:,i], D[:,:,i])), numpy.array(XPT)).ravel()
    Yn = numpy.matmul(numpy.array(numpy.multiply(D[:,:,i], D[:,:,i])), numpy.array(Xn)).ravel()

Random = numpy.zeros((npoints,1000,100,10))
Avg_Pr_c = 0
mean = YPT/(pow(4*3.14/wave_l, 2))
for i in range(0,100):
    for j in range(0,npoints):
        Random[j,:,i,:] = numpy.random.exponential(mean[j,i],(1000,10))

for xxxx in Thre:
    Avg_Pr_c = 0

    for j in range(0,100):
        idd = []
        d_er = []
        e_d = pproc-numpy.matlib.repmat(XY_dar[j,:],npoints,1)
        e = numpy.multiply(e_d, e_d)
        d = numpy.sqrt(numpy.sum(e,1))
        d1 = d[em_X_cod::n1]
    
        for kk in range(0,40):
            d_er_t = min(d1)
            idd_t = (d1.tolist()).index(d_er_t)
            #[d_er_t,idd_t]= min(d1)
            idd.append(idd_t)
            d_er.append(d_er_t)
            d1[idd_t] = 100000

        id = (numpy.array(idd)-1)+em_X_cod
        YNJ = sum(Yn[id])
        emitter = []
        emitter.append(pproc[id,0]-XY_dar[j,0]) 
        emitter.append(pproc[id,1]-XY_dar[j,1])
        S = numpy.array(S)
        Rad_pr0 = numpy.power(S[j,:], -2)
    
        Rad_pr1 = (R_c*R_c)*(numpy.power(S[j,:], -4))
        print("Rad_pr1: ", numpy.shape(Rad_pr1))
        check = S[j,:]<0.1
        check[id]= 0
        print("check: ", numpy.shape(check))
        s_t0 = numpy.multiply(check, Rad_pr0)
        print("s_t0: ", numpy.shape(s_t0))

        check = S[j,:]>=0.1
        check[id] = 0
        print("check2: ", numpy.shape(check))
        s_t1 = numpy.multiply(check, Rad_pr1)
        print("s_t1: ", numpy.shape(s_t1))

        check = numpy.zeros((npoints,1))
        check[id] = 1
        s_p = check
        for kk in range(0,numpy.size(id)):
            if S[id[kk],j]<0.1:
                s_p[id[kk]] = s_p[id[kk]]*Rad_pr0[id[kk]]
            else:
                s_p[id[kk]] = s_p[id[kk]]*Rad_pr1[id[kk]]
    
        print("s_p: ", numpy.shape(s_p))
        Pr_c = 0
        for c in range(0,10):
            s0 = numpy.multiply(Random[:,:,j,c], (numpy.matlib.repmat(s_t0,1000,1).transpose()))
            s1 = numpy.multiply(Random[:,:,j,c], (numpy.matlib.repmat(s_t1,1000,1)).transpose())
    
            I_RE = sum(s0+s1,1)
    
            P_R_RE = sum(numpy.multiply(Random[:,:,j,c], numpy.matlib.repmat(s_p,1,1000)),1)
        
            SNR_RE = numpy.divide(P_R_RE, (I_RE+YNJ))
            SNR_dB_RE = 10*numpy.log10(SNR_RE)
            #Result = interp1(Array,Array,SNR_dB_RE,'nearest','extrap')
            Result = interp1d(Array,Array)(SNR_dB_RE)
            #Result_uni = unique(Result)
            Result_uni = numpy.unique(Result)
            [y_n,loc] = ismember(Result_uni,Array)

            for i in range(0,numpy.size(loc)):
                num = Result_uni[i]
                count[loc[i]] = sum(Result==num)
        
            Area = sum(count)*diffs
            dist = numpy.divide(count, Area)

            Thresh = xxxx
            lst = [abs(i-Thresh) for i in Array]
            t = min(lst)
            ind = lst.index(t)
            '''[t,ind] = min([abs(i-Thresh) for i in Array])'''
   
            sm = sum(dist[ind:10000])*diffs
    
            Pr_c = Pr_c+sm

        Avg_Pr_c = Avg_Pr_c+Pr_c/10

    print(xxxx)
    #disp(xxxx)
    Avg1 = Avg_Pr_c/100
    emp.append(Avg1)

P_T = 4000
freq = 1.5*pow(10, 9)
wave_l = (3*pow(10, 8))/freq
p11 = 0.4+numpy.random.rand(100,1)*0.2
p22 = 0.4+numpy.random.rand(100,1)*0.2
receiver = []
receiver.append(p11) 
receiver.append(p22)





#AWGN model


print("Thre: ", Thre)
print("emp: ", emp)
#figure(1)
plt.plot(Thre,emp)
plt.show()

X = numpy.array(med_proc).transpose()
fileID = open('Values_T_Model1.txt','a+')
fileID.write('P_c = [%f %f %f %f %f %f %f %f %f],         X = [%f %f;%f %f],      SINR_AWGN = %f\n',emp,X,SNR_db_avg)
fileID.close()


'''Equations using Stochastic Geometry 2D
Mul = wave_l*3.14*((R_c)^2);

func = @(rxxx) exp((-Mul*I1(TT,rxxx,2,2,4))-(TT*noise*(R_c^2)*rxxx));
func1 = @(rxx) exp((-Mul*hyp_geo(-2,TT)*rxx)-(TT*noise*(R_c^2)*(rxx^2)));

func = @(rxxx) exp((-Mul*I(TT,rxxx,2,2,4)));
func1 = @(rxx) exp((-Mul*hyp_geo(-2,TT)*rxx));

Int_func1 = integral(func1,1,Inf);

ls1 = linspace(1,10000,100000);
ls = linspace(0.000001,1,1000);
fs_p1 = [];
fs_p = [];

for i = 1:1000
    fs_p = [fs_p func(ls(i))];
end

for i = 1:100000
    fs_p1 = [fs_p1 func1(ls1(i))];
end

Int_func = sum(fs_p*abs((ls(2)-ls(1))));
Int_func1 = sum(fs_p1*abs((ls1(2)-ls1(1))));
Sim_AvgPr = Mul*(Int_func+Int_func1)'''