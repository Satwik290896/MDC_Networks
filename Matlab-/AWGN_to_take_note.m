clear all;
lambda=200;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2)-0.5;
n1 = 2;
med_proc = [0.89,0.98;0.019,0.04038]
XY_dar = 0.5*rand(100,2);
S = [];
em_X_cod = 2;
emitter_X = med_proc(em_X_cod,:);

emp = [];

D = zeros(npoints,n1,100);
P_T = 4000;
freq = 1.5*10^9;
wave_l = (3*10^8)/freq;
noise = 0.1;
Thre = linspace(-25,-10,4);
%AWGN
SNR_avg = 0;

for j = 1:100
    
    I = 0;
    %e_d = pproc-repmat(XY_dar(j,:),npoints,1);
    %e = e_d.*e_d;
    %d = sqrt(sum(e,2));
    %[d_er,id]= min(d);
    id1 = em_X_cod;
    %if id1==0
    %    id1 = n1;
    %end
    emitter = med_proc(id1,:);
    %d_er=sqrt(sum(med_proc(em_X_cod,:).^2),2);
    d_er = sqrt(sum(emitter.^2));
    P_R = P_T * (wave_l/(4*3.14*d_er))^2;
    
    radd = sqrt(sum(med_proc.^2,2));
    for i = 1:n1
        rad = radd(i);
        
        if rad<0.1 && rad~=d_er
            alpha = 2;
            I = I + P_T*(wave_l/(4*3.14*rad))^2;
        end
        if(rad>=0.1 && rad~=d_er)
            alpha = 4;
            I = I + (P_T*(wave_l/(4*3.14*0.1))^2) + (P_T*(0.1/rad)^alpha);
        end
    end
    SNR = P_R/(I+noise);
    SNR_avg = SNR_avg+ SNR;
end
    
SNR_avg_AWGN = SNR_avg/100;
SNR_db_avg = 10*log10(SNR_avg_AWGN);
fileID = fopen('Values_X.txt','a+');
X = med_proc';
%fprintf(fileID,'X = [%f %f;%f %f],                 SINR = %f\n',X,SNR_db_avg);
fclose(fileID);
