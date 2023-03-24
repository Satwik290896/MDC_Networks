clear all; 
lambda=3000;
figure
npoints = poissrnd(lambda);
pproc = rand(npoints, 2);
p1 = 0.4+rand*0.2;
p2 = 0.4+rand*0.2;
receiver = [p1 p2];
plot(pproc(:,1)',pproc(:,2),'.g',p1,p2,'.r')
hold on
r = 0.1;
circle(p1,p2,r);

e_d = pproc-repmat(receiver,npoints,1);

e = e_d.*e_d;
d = sqrt(sum(e,2));
%b = d(d<0.03);
[d_er,id]= min(d)
%pi = randi([1,length(b)]);
%where = d == b(pi,1);
%id = find(where);
emitter = [pproc(id,1) pproc(id,2)];
plot(pproc(id,1),pproc(id,2),'.b');

%d_er =  sqrt((pproc(id,1)-p1)^2+(pproc(id,2)-p2)^2)

P_T = 4000;
freq = 1.5*10^9;
wave_l = (3*10^8)/freq;

%AWGN model

I = 0;
for i = 1:npoints
    rad = sqrt((pproc(i,1)-p1)^2+(pproc(i,2)-p2)^2)
    if rad<0.1 && rad > d_er
        alpha = 2;
        I = I + P_T*(wave_l/(4*3.14*rad))^2;
    end
    if(rad>=0.1)
        alpha = 4;
        I = I + (P_T*(wave_l/(4*3.14*0.1))^2) + (P_T*(0.1/rad)^alpha);
    end
end
        
noise = 0.1;
P_R = P_T * (wave_l/(4*3.14*d_er))^2;

SNR = P_R/(I+noise);

SNR_dB = 10*log10(SNR);






%Rayleigh fading

T = -75
T_end = 30
len = 10000;
Array = linspace(T,T_end,len);
mean = P_T/((4*3.14/wave_l)^2);
R_c = 0.1;

count = zeros(10000,1);
Pr_c = 0;

for c = 1:10
    for reps = 1:3000
        Random = exprnd(mean,[npoints,1]);
    
        rad_ar = sqrt(sum((pproc-repmat(receiver,npoints,1)).^2,2));
        Rad_pr0 = (rad_ar.^(-2));
        Rad_pr1 = (R_c^2)*(rad_ar.^(-4));
    
        check = rad_ar<0.1 ;
        che_ar = check .* rad_ar;
        check = che_ar>d_er;
        s_t0 = check.*Rad_pr0;
        s0 = Random.*(s_t0);
    
        check = rad_ar>=0.1;
        s_t1 = check.*Rad_pr1;
        s1 = Random.*s_t1;
    
        I_RE = sum(s0+s1);
    
        check = rad_ar == d_er;
        P_R_RE = sum(Random.*Rad_pr0.*check);
    
        SNR_RE = P_R_RE/(I_RE+noise);
        SNR_dB_RE = 10*log10(SNR_RE);
        [t,index] =  min(abs(Array-SNR_dB_RE));
        count(index) = count(index)+1;
    
    end
    
    diffs = abs(Array(1)-Array(2));
    Area = sum(count)*diffs;
    dist = count./Area;
    
    %figure(2)
    %plot(1:10000,dist)

    Thresh = -11;
    [t,ind] = min(abs(Array-Thresh));
   
    sm = sum(dist(ind:10000))*diffs;
    
    Pr_c = Pr_c+sm;
    disp(c)
end

Avg_Pr_c = Pr_c/10

