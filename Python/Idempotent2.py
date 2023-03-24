import numpy as np
import scipy.stats as stats
import numpy as np
from scipy.optimize import lsq_linear
import cvxpy as cp

def Correlating_function(npoints, med_proc, pp, XY_disp, n1):
    pproc = pp - np.tile(XY_disp, (npoints, 1))
    A = np.zeros((n1, n1))

    for i in range(n1-1):
        A[i, i+1] = 1

    C = np.transpose(med_proc)
    Dt = np.transpose(pproc)
    b1 = np.zeros((n1, 1)) + 0.00001
    b1[n1-1, 0] = 0
    b2_t = np.zeros((n1, 1))
    b3_t = np.zeros((n1, 1)) + 0.00001
    D = []
    Array = np.zeros((n1-2, n1))
    C1 = np.concatenate((C, Array), axis=0)

    for i in range(npoints):
        j = i % n1
        if j == 0:
            j = n1
        A1 = A
        b2 = b2_t
        b3 = b3_t
        b3[j-1] = 10
        b2[j-1] = -10

        if j != 1:
            A1[j-2, :] = np.concatenate(([1], np.zeros(n1-2)), axis=0)

        print(i)
        #constraints=[{"type": "ineq", "fun": lambda x: A1@x-b1.flatten()}]
        B = np.concatenate((Dt[:,i], np.zeros(n1-2)))
        print(np.shape(C1), np.shape(B), np.shape(b2), np.shape(b3))
        
        res = lsq_linear(C1, B, 
                         bounds=(b2.flatten(), b3.flatten()), lsmr_tol='auto', verbose=1)
        Xt = res.x
        D.append(Xt)

    return np.transpose(D)


def ismember(a_vec, b_vec):
    """ MATLAB equivalent ismember function """

    bool_ind = np.isin(a_vec,b_vec)
    common = a_vec[bool_ind]
    common_unique, common_inv  = np.unique(common, return_inverse=True)     # common = common_unique[common_inv]
    b_unique, b_ind = np.unique(b_vec, return_index=True)  # b_unique = b_vec[b_ind]
    common_ind = b_ind[np.isin(b_unique, common_unique, assume_unique=True)]
    return bool_ind, common_ind[common_inv]

lamda = 200
npoints = stats.poisson(lamda).rvs()
pproc = np.random.rand(npoints, 2) - 0.5
n1 = 2
med_proc = np.array([[0.99, 0.99], [0.000000001, 0.000000001]])
XY_dar = 0.5 * np.random.rand(100, 2)
S = []
em_X_cod = 2 - 1
emitter_X = med_proc[em_X_cod]

emp = []

D = np.zeros((npoints, n1, 100))
P_T = 4000
freq = 1.5 * 10 ** 9
wave_l = (3 * 10 ** 8) / freq
noise = 0.1
Thre = np.linspace(-30, 10, 9)
SNR_avg = 0

for j in range(100):
    I = 0
    id1 = em_X_cod
    emitter = med_proc[id1]
    d_er = np.sqrt(np.sum(emitter ** 2))
    P_R = P_T * (wave_l / (4 * 3.14 * d_er)) ** 2

    radd = np.sqrt(np.sum(med_proc ** 2, axis=1))
    for i in range(n1):
        rad = radd[i]

        if rad < 0.1 and rad != d_er:
            alpha = 2
            I = I + P_T * (wave_l / (4 * 3.14 * rad)) ** 2
        if rad >= 0.1 and rad != d_er:
            alpha = 4
            I = I + (P_T * (wave_l / (4 * 3.14 * 0.1)) ** 2) + (P_T * (0.1 / rad) ** alpha)
    SNR = P_R / (I + noise)
    SNR_avg = SNR_avg + SNR

SNR_avg_AWGN = SNR_avg / 100
SNR_db_avg = 10 * np.log10(SNR_avg_AWGN)

T = -75
T_end = 30
len = 10000
Array = np.linspace(T, T_end, len)
diffs = np.abs(Array[0] - Array[1])

R_c = 0.1
noise = 0.1

count = np.zeros((10000, 1))

med_temp = np.sum((med_proc ** 2), axis=1)
XPT = np.repeat(P_T, n1)
YPT = np.zeros((npoints, 100))
Xn = np.repeat(noise, n1)
for i in range(100):
    D[:, :, i] = (Correlating_function(npoints, med_proc, pproc, XY_dar[i], n1)).transpose()
    S.append(np.sqrt((D[:, :, i] ** 2) @ med_temp))
    YPT[:, i] = (D[:, :, i] * D[:, :, i]) @ XPT
    Yn = (D[:, :, i] * D[:, :, i]) @ Xn

S = np.array(S).transpose()
print("S: ", np.shape(S))
Random = np.zeros((npoints, 1000, 100, 10))
Avg_Pr_c = 0

mean = YPT / ((4 * 3.14 / wave_l) ** 2)
Random = np.zeros((npoints, 1000, 100, 10))

for i in range(100):
    for j in range(npoints):
        Random[j, :, i, :] = np.random.exponential(mean[j, i], size=(1000, 10))

Thre = [-30, -25]
for xxxx in Thre:
    Avg_Pr_c = 0

    for j in range(1, 101):
        idd = []
        d_er = []
        e_d = pproc - np.tile(XY_dar[j-1,:], (npoints,1))
        e = e_d * e_d
        d = np.sqrt(np.sum(e, axis=1))
        d1 = d[em_X_cod-1::n1]

        for kk in range(1, 41):
            idd_t = np.argmin(d1)
            d_er_t = d1[idd_t]
            idd.append(idd_t)
            d_er.append(d_er_t)
            d1[idd_t] = 100000

        id = np.array(idd) + em_X_cod - 1
        YNJ = np.sum(Yn[id])
        #emitter = np.concatenate((pproc[id,0]-XY_dar[j-1,0], pproc[id,1]-XY_dar[j-1,1]), axis=1)
        Rad_pr0 = (S[:,j-1]**(-2))
        Rad_pr1 = (R_c**2) * (S[:,j-1]**(-4))
        check = S[:,j-1] < 0.1
        check[id] = 0
        s_t0 = check * Rad_pr0

        check = S[:,j-1] >= 0.1
        check[id] = 0
        s_t1 = check * Rad_pr1

        check = np.zeros((npoints,1))
        check[id] = 1
        s_p = check

        print(type(id))
        for nn in range(0, np.size(id)):
            if S[id[nn],j-1] < 0.1:
                s_p[id[nn]] = s_p[id[nn]] * Rad_pr0[id[nn]]
            else:
                s_p[id[nn]] = s_p[id[nn]] * Rad_pr1[id[nn]]

        Pr_c = 0
        #print("s_t0:", np.shape(s_t0))
        #print("s_t1:", np.shape(s_t1))
        for c in range(1, 11):
            s0 = Random[:,:,j-1,c-1] * np.tile(s_t0, (1000,1)).transpose()
            s1 = Random[:,:,j-1,c-1] * np.tile(s_t1, (1000,1)).transpose()

            I_RE = np.sum(s0 + s1, axis=0)
            #print("Random: ", np.shape(Random[:,:,j-1,c-1]))
            P_R_RE = np.sum(Random[:,:,j-1,c-1] * np.tile(s_p, (1,1000)), axis=0)

            SNR_RE = P_R_RE / (I_RE + YNJ)
            SNR_dB_RE = 10 * np.log10(SNR_RE)
            Result = np.interp(SNR_dB_RE, Array, Array, left=None, right=None, period=None)
            Result_uni = np.unique(Result)
            #loc = np.isin(Array, Result_uni)
            [y_n,loc] = ismember(Result_uni,Array)

            count = np.zeros_like(Array)
            for i in range(0,np.size(loc)):
                #if loc[i]:
                #num = Result_uni[np.argwhere(Result_uni == Array[i])]
                num = Result_uni[i]
                count[loc[i]] = np.sum(Result == num)

            Area = np.sum(count) * diffs
            dist = count / Area

            Thresh = xxxx
            ind = np.argmin(np.abs(Array - Thresh))
            print(xxxx)
            print(ind)
            sm = np.sum(dist[ind:10000]) * diffs
            Pr_c += sm

        Avg_Pr_c += Pr_c / 10

    print(xxxx)
    print(ind)
    Avg1 = Avg_Pr_c/100
    emp.append(Avg1)


print("Thre: ", Thre)
print("emp: ", emp)