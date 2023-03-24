clear all;
lambda=200;
npoints = poissrnd(lambda);
pproc = rand(npoints, 2)-0.5;

figure(1)
plot(pproc(:,1),pproc(:,2),'.');

%very rarely plotted
if npoints>1000
    med_proc_1 = rand(1000,2);
    med_proc_2 = (1/1000)*rand(npoints-1000,2);
    med_proc = [med_proc_1;med_proc_2];
else
    med_proc = rand(npoints-160,2);
end

A = zeros(npoints-160,npoints-160);

for i=1:npoints-161
    A(i,i+1) = 1;
end

C = med_proc';
Dt = pproc';
b1 = zeros(npoints-160,1)+0.00001;
b1(npoints-160,1) = 0 ;
b2_t = zeros(npoints-160,1);
b3_t = zeros(npoints-160,1)+0.00001;
D = [];
Array = zeros(npoints-162,npoints-160);
C1 = [C;Array]; 

for i = 1:npoints
    j = rem(i,npoints-160);
    if j==0
        j = npoints-160;
    end
    A1 = A;
    b2 = b2_t;
    b3 = b3_t;
    b3(j) = 10;
    b2(j) = -10;

    if j~=1 
        A1(j-1,:) = [1 zeros(1,npoints-161)];
    end
    disp(i);
    Xt = lsqlin(C1,[Dt(:,i);zeros(npoints-162,1)],A1,b1,[],[],b2,b3);
    D = [D;Xt'];
end

%D = pproc/med_proc;

%summation_of_distances

                                