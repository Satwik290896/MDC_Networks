
%First term alone taken and Condition of Model 1
%Simulation technique: 40 matched points were taken

clear all;
lambda=200;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2)-0.5;
n1 = 2;
med_proc = [0.990000 0.990000;0.000000001 0.000000001];
XY_dar = 0.5*rand(100,2);
S = [];
em_X_cod = 2;
emitter_X = med_proc(em_X_cod,:);

emp = [];

D = zeros(npoints,n1,100);
P_T = 4000;
freq = 1.5*10^9;
wave_l = (3*10^8)/freq;
SINR_AWGN = 185; % in dB
noise = 10.^(SINR_AWGN/10);%Conversion from dB

%noise = 0.1;
% Thre = linspace(-30,10,9);
Thre = linspace(-30,15,300)
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


T = -75
T_end = 30
len = 10000;
Array = linspace(T,T_end,len);
diffs = abs(Array(1)-Array(2));

R_c = 0.1;
noise = 10.^(SINR_AWGN/10);%Conversion from dB
% noise = 0.1;
% -5dB
%noise = 0.3162;

count = zeros(10000,1);

med_temp = sum((med_proc.^2),2);
XPT = repmat(P_T,n1,1);
YPT = zeros(npoints,100);
Xn = repmat(noise,n1,1);
for i=1:100
    D(:,:,i) = Correlating_function(npoints,med_proc,pproc,XY_dar(i,:),n1);
    S = [S sqrt((D(:,:,i).^2)*med_temp)];
    YPT(:,i) = (D(:,:,i).*D(:,:,i))*XPT;
    Yn = (D(:,:,i).*D(:,:,i))*Xn;
end

Random = zeros(npoints,1000,100,10);
Avg_Pr_c = 0;
mean = YPT/((4*3.14/wave_l)^2);
for i = 1:100
    for j = 1:npoints
        Random(j,:,i,:) = exprnd(mean(j,i),[1,1000,1,10]);
    end
end
jarvis = 0;
for xxxx = Thre
    jarvis = jarvis + 1;
    disp(jarvis);
    Avg_Pr_c = 0;
for j = 1:10
    idd = [];
    d_er = [];
    e_d = pproc-repmat(XY_dar(j,:),npoints,1);
    e = e_d.*e_d;
    d = sqrt(sum(e,2));
    d1 = d(em_X_cod:n1:end);
    for kk = 1:40
        [d_er_t,idd_t]= min(d1);
        idd = [idd idd_t];
        d_er = [d_er d_er_t];
        d1(idd_t) = 100000;
    end
    id = (idd-1)+em_X_cod;
    %[d_er,idd]= min(d1);
    %id = (idd-1)+em_X_cod;
    YNJ = sum(Yn(id)); 
    emitter = [pproc(id,1)-XY_dar(j,1) pproc(id,2)-XY_dar(j,2)];
    %rad_ar = sqrt(sum((pproc-repmat(receiver(j,:),npoints,1)).^2,2));
    Rad_pr0 = (S(:,j).^(-2));
    Rad_pr1 = (R_c^2)*(S(:,j).^(-4));
    check = S(:,j)<0.1 ;
    %$che_ar = check .* S(:,j);
    %check = che_ar>d_er;
    check(id)= 0;
    s_t0 = check.*Rad_pr0;

    check = S(:,j)>=0.1;
    check(id) = 0;
    s_t1 = check.*Rad_pr1;

    check = zeros(npoints,1);
    check(id) = 1;
    s_p = check;
    for kk = 1:length(id)
        if S(id(kk),j)<0.1
            s_p(id(kk)) = s_p(id(kk))*Rad_pr0(id(kk));
        else
            s_p(id(kk)) = s_p(id(kk))*Rad_pr1(id(kk));
        end
    end
    
    Pr_c = 0;
    for c = 1:10
        s0 = Random(:,:,j,c).*(repmat(s_t0,1,1000));
        s1 = Random(:,:,j,c).*repmat(s_t1,1,1000);
    
        I_RE = sum(s0+s1,1);
    
        %if S(id,j)<0.1
        P_R_RE = sum(Random(:,:,j,c).*repmat(s_p,1,1000),1);
        %else
        %    P_R_RE = sum(Random(:,:,j,c).*repmat(Rad_pr1.*check,1,1000),1);
        %end
        
        SNR_RE = P_R_RE./(I_RE+YNJ);
        SNR_dB_RE = 10*log10(SNR_RE);
        Result = interp1(Array,Array,SNR_dB_RE,'nearest','extrap');
        Result_uni = unique(Result);
        [y_n,loc] = ismember(Result_uni,Array);
        for i = 1:length(loc)
            num = Result_uni(i);
            count(loc(i)) = sum(Result==num);
        end
        Area = sum(count)*diffs;
        dist = count./Area;

        Thresh = xxxx;
        [t,ind] = min(abs(Array-Thresh));
   
        sm = sum(dist(ind:10000))*diffs;
    
        Pr_c = Pr_c+sm;
        %disp(c)
    end
    Avg_Pr_c = Avg_Pr_c+Pr_c/10;
    
end

disp(xxxx);
Avg1 = Avg_Pr_c/10;
disp(Avg1);
emp = [emp Avg1];
end

%TT = 10^(Thresh/10);



P_T = 4000;
freq = 1.5*10^9;
wave_l = (3*10^8)/freq;
p11 = 0.4+rand(100,1)*0.2;
p22 = 0.4+rand(100,1)*0.2;
receiver = [p11 p22];





%AWGN model



figure(1)
plot(Thre,emp);

X = med_proc';
fileID = fopen('Values_T_Model1_S5_185dB.txt','a+');
fprintf(fileID,'%f \n',emp);
fclose(fileID);


%Equations using Stochastic Geometry 2D
%Mul = wave_l*3.14*((R_c)^2);

%func = @(rxxx) exp((-Mul*I1(TT,rxxx,2,2,4))-(TT*noise*(R_c^2)*rxxx));
%func1 = @(rxx) exp((-Mul*hyp_geo(-2,TT)*rxx)-(TT*noise*(R_c^2)*(rxx^2)));

%func = @(rxxx) exp((-Mul*I(TT,rxxx,2,2,4)));
%func1 = @(rxx) exp((-Mul*hyp_geo(-2,TT)*rxx));

%Int_func1 = integral(func1,1,Inf);

%ls1 = linspace(1,10000,100000);
%ls = linspace(0.000001,1,1000);
%fs_p1 = [];
%fs_p = [];

%for i = 1:1000
%    fs_p = [fs_p func(ls(i))];
%end

%for i = 1:100000
%    fs_p1 = [fs_p1 func1(ls1(i))];
%end

%Int_func = sum(fs_p*abs((ls(2)-ls(1))));
%Int_func1 = sum(fs_p1*abs((ls1(2)-ls1(1))));
%Sim_AvgPr = Mul*(Int_func+Int_func1)