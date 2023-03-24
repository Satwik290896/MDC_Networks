%clear all; 
lambda=3000;
%figure
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);


%plot(pproc(:,1)',pproc(:,2),'.g',p1,p2,'.r')
%hold on

%r = 0.1;
%circle(p1,p2,r);


%b = d(d<0.03);

%pi = randi([1,length(b)]);
%where = d == b(pi,1);
%id = find(where);

%plot(pproc(id,1),pproc(id,2),'.b');

%d_er =  sqrt((pproc(id,1)-p1)^2+(pproc(id,2)-p2)^2)

P_T = 4000;
freq = 1.5*10^9;
wave_l = (3*10^8)/freq;







%AWGN model

%I = 0;
%for i = 1:npoints
%    rad = sqrt((pproc(i,1)-p1)^2+(pproc(i,2)-p2)^2)
%    if rad<0.1 && rad > d_er
%        alpha = 2;
%        I = I + P_T*(wave_l/(4*3.14*rad))^2;
%    end
%    if(rad>=0.1)
%        alpha = 4;
%        I = I + (P_T*(wave_l/(4*3.14*0.1))^2) + (P_T*(0.1/rad)^alpha);
%    end
%end
        
noise = 0.1;
%P_R = P_T * (wave_l/(4*3.14*d_er))^2;

%SNR = P_R/(I+noise);

%SNR_dB = 10*log10(SNR);


clear all; 
lambda=3000;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);
P_T = 4000;
freq = 1.5*10^9;
wave_l = (3*10^8)/freq;
noise = 0.1;

%Rayleigh fading

T = -75
T_end = 30
len = 10000;
Array = linspace(T,T_end,len);
diffs = abs(Array(1)-Array(2));
mean = P_T/((4*3.14/wave_l)^2);
R_c = 0.1;

count = zeros(10000,1);

Avg_Pr_c = 0;


%for reps = 1:3000
%    R = exprnd(mean,[npoints,1]);
%    Random = [Random R];
%end

Random = exprnd(mean,[npoints,1000,100]);

for j = 1:100
    p1 = 0.4+rand*0.2;
    p2 = 0.4+rand*0.2;
    receiver = [p1 p2];
    e_d = pproc-repmat(receiver,npoints,1);
    e = e_d.*e_d;
    d = sqrt(sum(e,2));
    [d_er,id]= min(d);
    emitter = [pproc(id,1) pproc(id,2)];
    rad_ar = sqrt(sum((pproc-repmat(receiver,npoints,1)).^2,2));
    Rad_pr0 = (rad_ar.^(-2));
    Rad_pr1 = (R_c^2)*(rad_ar.^(-4));
    check = rad_ar<0.1 ;
    che_ar = check .* rad_ar;
    check = che_ar>d_er;
    s_t0 = check.*Rad_pr0;

    check = rad_ar>=0.1;
    s_t1 = check.*Rad_pr1;

    check = rad_ar == d_er;
    Pr_c = 0;
    for c = 1:100
        s0 = Random(:,:,c).*(repmat(s_t0,1,1000));
        s1 = Random(:,:,c).*repmat(s_t1,1,1000);
    
        I_RE = sum(s0+s1,1);
    
        
        P_R_RE = sum(Random(:,:,c).*repmat(Rad_pr0.*check,1,1000),1);
    
        SNR_RE = P_R_RE./(I_RE+noise);
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

        Thresh = -11;
        [t,ind] = min(abs(Array-Thresh));
   
        sm = sum(dist(ind:10000))*diffs;
    
        Pr_c = Pr_c+sm;
        %disp(c)
    end
    Avg_Pr_c = Avg_Pr_c+Pr_c/100;
    disp(j);
end

Avg1 = Avg_Pr_c/100


TT = 10^(Thresh/10);

%Equations using Stochastic Geometry 2D
Mul = wave_l*3.14*((R_c)^2);

func = @(rxxx) exp((-Mul*I(TT,rxxx,2,2,4))-(TT*noise*(R_c^2)*rxxx));
func1 = @(rxx) exp((-Mul*hyp_geo(-2,TT)*rxx)-(TT*noise*(R_c^2)*(rxx^2)));

%func = @(rxxx) exp((-Mul*I(TT,rxxx,2,2,4)));
%func1 = @(rxx) exp((-Mul*hyp_geo(-2,TT)*rxx));

%Int_func1 = integral(func1,1,Inf);

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
Sim_AvgPr = Mul*(Int_func+Int_func1)


